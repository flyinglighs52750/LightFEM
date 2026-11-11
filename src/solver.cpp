#include "lightfem/solver.hpp"
#include "lightfem/bench_mode.hpp"
#include "lightfem/precond.hpp"
#include "lightfem/reordering.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#ifdef LIGHTFEM_HAS_OPENMP
#include <omp.h>
#endif

namespace lightfem {

namespace {

struct AssembledSystem {
    CsrMatrix stiffness;
    std::vector<double> load;
    std::size_t threads_used{1};
    std::size_t coloring_count{0};
    double pattern_ms{0.0};
    double element_stiffness_ms{0.0};
    double assembly_ms{0.0};
};

double dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) throw std::invalid_argument("dot product vector size mismatch");
    double result = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) result += left[i] * right[i];
    return result;
}

double norm(const std::vector<double>& vector) { return std::sqrt(std::max(0.0, dot(vector, vector))); }

double elapsed_ms(const std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

std::array<Vec2, 4> q4_coordinates(const Model& model, const Q4Element& element) {
    std::array<Vec2, 4> coordinates{};
    for (std::size_t i = 0; i < 4; ++i) {
        const Node* node = model.find_node(element.node_ids[i]);
        if (node == nullptr) throw std::logic_error("validated element references an unknown node");
        coordinates[i] = node->position;
    }
    return coordinates;
}

// Q4 elements come first, then bars. Every consumer below iterates in this order.
std::vector<std::vector<int>> element_node_lists(const Model& model) {
    std::vector<std::vector<int>> lists;
    lists.reserve(model.elements.size() + model.trusses.size());
    for (const Q4Element& element : model.elements) {
        lists.push_back({element.node_ids[0], element.node_ids[1], element.node_ids[2], element.node_ids[3]});
    }
    for (const TrussElement& truss : model.trusses) {
        lists.push_back({truss.node_ids[0], truss.node_ids[1]});
    }
    return lists;
}

ElementDofList element_dof_lists(const Model& model) {
    ElementDofList element_dofs;
    element_dofs.reserve(model.elements.size() + model.trusses.size());
    for (const Q4Element& element : model.elements) {
        const auto dofs = q4_dof_indices(model, element);
        element_dofs.emplace_back(dofs.begin(), dofs.end());
    }
    for (const TrussElement& truss : model.trusses) {
        const auto dofs = truss_dof_indices(model, truss);
        element_dofs.emplace_back(dofs.begin(), dofs.end());
    }
    return element_dofs;
}

std::size_t configure_openmp(const Model& model) {
#ifdef LIGHTFEM_HAS_OPENMP
    if (model.settings.threads > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("threads exceeds OpenMP supported range");
    }
    if (model.settings.threads != 0) omp_set_num_threads(static_cast<int>(model.settings.threads));
    return static_cast<std::size_t>(omp_get_max_threads());
#else
    if (model.settings.assembly != "serial") {
        throw std::runtime_error("parallel assembly was requested but LightFEM was built without OpenMP");
    }
    return 1;
#endif
}

std::vector<std::vector<std::size_t>> element_colors(const Model& model) {
    const auto nodes_of = element_node_lists(model);
    const std::size_t total = nodes_of.size();
    std::unordered_map<int, std::vector<std::size_t>> node_elements;
    node_elements.reserve(model.nodes.size());
    for (std::size_t element = 0; element < total; ++element) {
        for (const int node_id : nodes_of[element]) node_elements[node_id].push_back(element);
    }

    std::vector<int> colors(total, -1);
    int maximum_color = -1;
    std::vector<std::size_t> unavailable_stamp(total, 0);
    std::size_t stamp = 0;
    for (std::size_t element = 0; element < total; ++element) {
        ++stamp;
        for (const int node_id : nodes_of[element]) {
            for (const std::size_t neighbour : node_elements.at(node_id)) {
                if (colors[neighbour] >= 0) {
                    unavailable_stamp[static_cast<std::size_t>(colors[neighbour])] = stamp;
                }
            }
        }
        int color = 0;
        while (unavailable_stamp[static_cast<std::size_t>(color)] == stamp) ++color;
        colors[element] = color;
        maximum_color = std::max(maximum_color, color);
    }

    std::vector<std::vector<std::size_t>> groups(static_cast<std::size_t>(maximum_color + 1));
    for (std::size_t element = 0; element < colors.size(); ++element) {
        groups[static_cast<std::size_t>(colors[element])].push_back(element);
    }
    return groups;
}

void serial_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                     const std::vector<std::vector<double>>& local_stiffness) {
    for (std::size_t element = 0; element < local_stiffness.size(); ++element) {
        const auto& mapping = pattern.element_value_indices[element];
        const auto& local = local_stiffness[element];
        for (std::size_t index = 0; index < local.size(); ++index) stiffness.values()[mapping[index]] += local[index];
    }
}

void atomic_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                     const std::vector<std::vector<double>>& local_stiffness) {
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
    for (long long element = 0; element < static_cast<long long>(local_stiffness.size()); ++element) {
        const auto& mapping = pattern.element_value_indices[static_cast<std::size_t>(element)];
        const auto& local = local_stiffness[static_cast<std::size_t>(element)];
        for (std::size_t index = 0; index < local.size(); ++index) {
#pragma omp atomic
            stiffness.values()[mapping[index]] += local[index];
        }
    }
#else
    (void)stiffness;
    (void)pattern;
    (void)local_stiffness;
#endif
}

void coloring_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                       const std::vector<std::vector<double>>& local_stiffness,
                       const std::vector<std::vector<std::size_t>>& colors) {
    for (const auto& group : colors) {
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long long local = 0; local < static_cast<long long>(group.size()); ++local) {
            const std::size_t element = group[static_cast<std::size_t>(local)];
            const auto& mapping = pattern.element_value_indices[element];
            const auto& local_values = local_stiffness[element];
            for (std::size_t index = 0; index < local_values.size(); ++index) {
                stiffness.values()[mapping[index]] += local_values[index];
            }
        }
    }
}

void thread_local_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                           const std::vector<std::vector<double>>& local_stiffness,
                           const std::size_t threads_used) {
    const std::size_t nnz = stiffness.nonzeros();
    if (nnz != 0 && threads_used > std::numeric_limits<std::size_t>::max() / nnz) {
        throw std::overflow_error("thread-local CSR assembly buffer size overflows");
    }
    std::vector<std::vector<double>> buffers(threads_used, std::vector<double>(nnz, 0.0));
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel
    {
        const std::size_t thread_id = static_cast<std::size_t>(omp_get_thread_num());
        std::vector<double>& local_values = buffers[thread_id];
#pragma omp for schedule(static)
        for (long long element = 0; element < static_cast<long long>(local_stiffness.size()); ++element) {
            const std::size_t element_index = static_cast<std::size_t>(element);
            const auto& mapping = pattern.element_value_indices[element_index];
            const auto& local = local_stiffness[element_index];
            for (std::size_t index = 0; index < local.size(); ++index) {
                local_values[mapping[index]] += local[index];
            }
        }
    }
#else
    (void)pattern;
    (void)local_stiffness;
#endif

#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long index = 0; index < static_cast<long long>(nnz); ++index) {
        double sum = 0.0;
        for (const auto& buffer : buffers) sum += buffer[static_cast<std::size_t>(index)];
        stiffness.values()[static_cast<std::size_t>(index)] = sum;
    }
}

AssembledSystem assemble(const Model& model) {
    const std::size_t dofs = model.summary().dofs;
    // The "scan" strategy is a measurement control: it assembles with the same
    // CSR pattern but resolves every entry by searching the row, so the gap to
    // the mapped path is exactly what the precomputed map buys. Its setup skips
    // the map, which makes the build cost of the map visible as the difference
    // between the two pattern timings.
    const bool use_precomputed_map = model.settings.assembly != "scan";
    const auto pattern_start = std::chrono::steady_clock::now();
    const ElementDofList element_dofs = element_dof_lists(model);
    CsrAssemblyPattern pattern;
    if (use_precomputed_map) {
        pattern = build_csr_assembly_pattern(dofs, element_dofs);
    } else {
        pattern.matrix = CsrMatrix::from_element_dofs(dofs, element_dofs);
    }
    const double pattern_ms = elapsed_ms(pattern_start);

    AssembledSystem system;
    system.pattern_ms = pattern_ms;
    system.stiffness = std::move(pattern.matrix);
    system.load.assign(dofs, 0.0);
    system.threads_used = configure_openmp(model);
    const std::size_t total_elements = model.elements.size() + model.trusses.size();
    std::vector<std::vector<double>> local_stiffness(total_elements);
    const auto stiffness_start = std::chrono::steady_clock::now();
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long element_index = 0; element_index < static_cast<long long>(total_elements); ++element_index) {
        const std::size_t e = static_cast<std::size_t>(element_index);
        if (e < model.elements.size()) {
            const Q4Element& element = model.elements[e];
            const Material* material = model.find_material(element.material_id);
            if (material == nullptr) throw std::logic_error("validated element references an unknown material");
            const ElementMatrix matrix =
                q4_stiffness(q4_coordinates(model, element), *material, model.settings.plane_mode);
            local_stiffness[e].assign(matrix.begin(), matrix.end());
        } else {
            const TrussElement& truss = model.trusses[e - model.elements.size()];
            const Material* material = model.find_material(truss.material_id);
            if (material == nullptr) throw std::logic_error("validated truss references an unknown material");
            std::array<Vec2, 2> coordinates{};
            for (std::size_t i = 0; i < 2; ++i) {
                const Node* node = model.find_node(truss.node_ids[i]);
                if (node == nullptr) throw std::logic_error("validated truss references an unknown node");
                coordinates[i] = node->position;
            }
            const auto matrix = truss_stiffness(coordinates, truss.area, material->young_modulus);
            local_stiffness[e].assign(matrix.begin(), matrix.end());
        }
    }
    const auto stiffness_end = std::chrono::steady_clock::now();
    system.element_stiffness_ms =
        std::chrono::duration<double, std::milli>(stiffness_end - stiffness_start).count();

    const auto assembly_start = std::chrono::steady_clock::now();
    if (model.settings.assembly == "serial") {
        serial_assembly(system.stiffness, pattern, local_stiffness);
    } else if (model.settings.assembly == "scan") {
        for (std::size_t e = 0; e < total_elements; ++e) {
            const std::vector<std::size_t>& element_dof_list = element_dofs[e];
            const std::size_t count = element_dof_list.size();
            for (std::size_t a = 0; a < count; ++a) {
                for (std::size_t b = 0; b < count; ++b) {
                    system.stiffness.add(element_dof_list[a], element_dof_list[b],
                                         local_stiffness[e][a * count + b]);
                }
            }
        }
    } else if (model.settings.assembly == "atomic") {
        atomic_assembly(system.stiffness, pattern, local_stiffness);
    } else if (model.settings.assembly == "coloring") {
        const auto colors = element_colors(model);
        system.coloring_count = colors.size();
        coloring_assembly(system.stiffness, pattern, local_stiffness, colors);
    } else if (model.settings.assembly == "thread_local") {
        thread_local_assembly(system.stiffness, pattern, local_stiffness, system.threads_used);
    } else {
        throw std::logic_error("validated model has an unknown assembly strategy");
    }
    const auto assembly_end = std::chrono::steady_clock::now();
    system.assembly_ms = std::chrono::duration<double, std::milli>(assembly_end - assembly_start).count();

    for (const NodalLoad& load : model.loads) {
        const Node* node = model.find_node(load.node_id);
        if (node == nullptr) throw std::logic_error("validated load references an unknown node");
        system.load[node->dof_x] += load.fx;
        system.load[node->dof_y] += load.fy;
    }
    const auto accumulate_element_load = [&model, &system](const Q4Element& element,
                                                            const std::array<double, 8>& local_load) {
        const auto dofs = q4_dof_indices(model, element);
        for (std::size_t i = 0; i < 8; ++i) system.load[dofs[i]] += local_load[i];
    };
    for (const BodyForce& force : model.body_forces) {
        const Q4Element* element = model.find_element(force.element_id);
        if (element == nullptr) throw std::logic_error("validated body force references an unknown element");
        const Material* material = model.find_material(element->material_id);
        if (material == nullptr) throw std::logic_error("validated element references an unknown material");
        accumulate_element_load(*element, q4_body_force(q4_coordinates(model, *element), *material,
                                                         force.bx, force.by));
    }
    for (const EdgeLoad& load : model.edge_loads) {
        const Q4Element* element = model.find_element(load.element_id);
        if (element == nullptr) throw std::logic_error("validated edge load references an unknown element");
        const Material* material = model.find_material(element->material_id);
        if (material == nullptr) throw std::logic_error("validated element references an unknown material");
        const auto coordinates = q4_coordinates(model, *element);
        const auto local_load = load.type == EdgeLoadType::Traction
                                    ? q4_edge_traction(coordinates, *material, load.edge, load.first, load.second)
                                    : q4_edge_pressure(coordinates, *material, load.edge, load.first);
        accumulate_element_load(*element, local_load);
    }
    return system;
}

CsrProfile profile_for_model(const Model& model) {
    return csr_profile(build_csr_assembly_pattern(model.summary().dofs, element_dof_lists(model)).matrix);
}

std::vector<std::pair<std::size_t, double>> collect_constraints(const Model& model) {
    std::vector<std::pair<std::size_t, double>> constraints;
    for (const DirichletBC& boundary : model.boundaries) {
        const Node* node = model.find_node(boundary.node_id);
        if (node == nullptr) throw std::logic_error("validated boundary references an unknown node");
        if (boundary.ux) constraints.emplace_back(node->dof_x, *boundary.ux);
        if (boundary.uy) constraints.emplace_back(node->dof_y, *boundary.uy);
    }
    std::sort(constraints.begin(), constraints.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    return constraints;
}

// Symmetric elimination: Kfc*uc is moved to the right-hand side, then both
// the constrained row and column are zeroed and the diagonal is set to one.
void apply_dirichlet(CsrMatrix& stiffness, std::vector<double>& load,
                     const std::vector<std::pair<std::size_t, double>>& constraints) {
    if (load.size() != stiffness.rows()) throw std::invalid_argument("Dirichlet RHS size mismatch");
    if (constraints.empty()) return;

    if (legacy_stage("LIGHTFEM_LEGACY_BOUNDARY")) {
        for (const auto& constraint : constraints) {
            const std::size_t constrained_dof = constraint.first;
            const double prescribed_value = constraint.second;
            for (std::size_t row = 0; row < stiffness.rows(); ++row) {
                if (row == constrained_dof) continue;
                if (!stiffness.contains(row, constrained_dof)) continue;
                const double column_value = stiffness.at(row, constrained_dof);
                if (column_value != 0.0) {
                    load[row] -= column_value * prescribed_value;
                    stiffness.set(row, constrained_dof, 0.0);
                }
            }
            const auto begin = stiffness.row_offsets()[constrained_dof];
            const auto end = stiffness.row_offsets()[constrained_dof + 1];
            for (std::size_t index = begin; index < end; ++index) stiffness.values()[index] = 0.0;
            stiffness.set(constrained_dof, constrained_dof, 1.0);
            load[constrained_dof] = prescribed_value;
        }
        return;
    }

    // The loop above visits every row once per constrained DOF, so it costs
    // O(constraints * dofs): 2.6e8 row lookups on a 400x400 mesh. Marking the
    // constrained columns first turns the elimination into one sweep over the
    // nonzero values, which is O(nnz).
    const std::size_t dofs = stiffness.rows();
    std::vector<char> constrained(dofs, 0);
    std::vector<double> prescribed(dofs, 0.0);
    for (const auto& constraint : constraints) {
        constrained[constraint.first] = 1;
        prescribed[constraint.first] = constraint.second;
    }

    auto& values = stiffness.values();
    const auto& columns = stiffness.column_indices();
    const auto& offsets = stiffness.row_offsets();
    for (std::size_t row = 0; row < dofs; ++row) {
        const std::size_t begin = offsets[row];
        const std::size_t end = offsets[row + 1];
        if (constrained[row] != 0) {
            for (std::size_t index = begin; index < end; ++index) values[index] = 0.0;
            stiffness.set(row, row, 1.0);
            load[row] = prescribed[row];
            continue;
        }
        for (std::size_t index = begin; index < end; ++index) {
            const std::size_t column = columns[index];
            if (constrained[column] == 0 || values[index] == 0.0) continue;
            load[row] -= values[index] * prescribed[column];
            values[index] = 0.0;
        }
    }
}

// Pre-optimization CG loop, kept only for the scaling study: it allocates a
// fresh product vector per iteration and sweeps the vectors seven times per
// step. The optimized path below fuses those sweeps.
LinearSolverReport pcg_jacobi_legacy(const CsrMatrix& matrix, const std::vector<double>& rhs,
                                     const double tolerance, const std::size_t max_iterations,
                                     std::vector<double>& solution) {
    if (rhs.size() != matrix.rows()) throw std::invalid_argument("PCG RHS size mismatch");
    solution.assign(matrix.rows(), 0.0);
    std::vector<double> residual = rhs;
    const std::vector<double> diagonal = matrix.diagonal();
    std::vector<double> inverse_diagonal(matrix.rows(), 0.0);
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        inverse_diagonal[i] = std::abs(diagonal[i]) > 1.0e-30 ? 1.0 / diagonal[i] : 1.0;
    }

    LinearSolverReport report;
    const auto clock_start = std::chrono::steady_clock::now();
    double matvec_total_ms = 0.0;
    std::size_t matvec_calls = 0;
    const auto finish = [&]() {
        report.solve_ms = elapsed_ms(clock_start);
        report.matvec_ms = matvec_total_ms;
        report.matvec_calls = matvec_calls;
        return report;
    };

    report.initial_residual = norm(residual);
    const double target = tolerance * std::max(1.0, report.initial_residual);
    if (report.initial_residual <= target) {
        report.converged = true;
        report.final_residual = report.initial_residual;
        report.message = "initial residual satisfies tolerance";
        return finish();
    }

    std::vector<double> preconditioned(matrix.rows(), 0.0);
    for (std::size_t i = 0; i < matrix.rows(); ++i) preconditioned[i] = inverse_diagonal[i] * residual[i];
    std::vector<double> direction = preconditioned;
    double rho = dot(residual, preconditioned);
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("PCG failed to form a positive initial preconditioned residual");
    }

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const auto matvec_start = std::chrono::steady_clock::now();
        const std::vector<double> matrix_direction = matrix.multiply(direction);
        matvec_total_ms += elapsed_ms(matvec_start);
        ++matvec_calls;
        const double curvature = dot(direction, matrix_direction);
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            throw std::runtime_error("PCG encountered non-positive curvature; system is not SPD");
        }
        const double alpha = rho / curvature;
        for (std::size_t i = 0; i < matrix.rows(); ++i) {
            solution[i] += alpha * direction[i];
            residual[i] -= alpha * matrix_direction[i];
        }
        report.iterations = iteration;
        report.final_residual = norm(residual);
        if (report.final_residual <= target) {
            report.converged = true;
            report.message = "converged";
            return finish();
        }
        for (std::size_t i = 0; i < matrix.rows(); ++i) preconditioned[i] = inverse_diagonal[i] * residual[i];
        const double next_rho = dot(residual, preconditioned);
        if (!(next_rho > 0.0) || !std::isfinite(next_rho)) {
            throw std::runtime_error("PCG preconditioned residual became non-positive or non-finite");
        }
        const double beta = next_rho / rho;
        for (std::size_t i = 0; i < matrix.rows(); ++i) direction[i] = preconditioned[i] + beta * direction[i];
        rho = next_rho;
    }

    report.final_residual = norm(residual);
    report.message = "maximum iterations reached without convergence";
    return finish();
}

// Diagonal scaling. Kept as a functor so the Krylov loop below is written once
// and can be instantiated against either storage scheme and either
// preconditioner.
class JacobiPreconditioner {
public:
    explicit JacobiPreconditioner(std::vector<double> inverse_diagonal)
        : inverse_diagonal_(std::move(inverse_diagonal)) {}

    void apply(const std::vector<double>& input, std::vector<double>& output) const {
        const std::size_t n = input.size();
        for (std::size_t i = 0; i < n; ++i) output[i] = inverse_diagonal_[i] * input[i];
    }

private:
    std::vector<double> inverse_diagonal_;
};

// Adds the timing a heavy preconditioner reports on its own behalf, so the
// Krylov loop stays free of accounting code.
template <typename Preconditioner>
class TimedPreconditioner {
public:
    explicit TimedPreconditioner(const Preconditioner& preconditioner)
        : preconditioner_(preconditioner) {}

    void apply(const std::vector<double>& input, std::vector<double>& output) const {
        const auto start = std::chrono::steady_clock::now();
        preconditioner_.apply(input, output);
        preconditioner_.account_apply(elapsed_ms(start));
    }

private:
    const Preconditioner& preconditioner_;
};

// Preconditioned CG. Both storage schemes expose rows() and multiply_into(), and
// both preconditioners expose apply(), so one template covers all four
// combinations and every combination runs identical arithmetic.
template <typename Matrix, typename Preconditioner>
LinearSolverReport pcg_solve(const Matrix& matrix, const Preconditioner& preconditioner,
                             const std::vector<double>& rhs, const double tolerance,
                             const std::size_t max_iterations, std::vector<double>& solution) {
    if (rhs.size() != matrix.rows()) throw std::invalid_argument("PCG RHS size mismatch");
    const std::size_t n = matrix.rows();
    solution.assign(n, 0.0);
    std::vector<double> residual = rhs;

    LinearSolverReport report;
    const auto clock_start = std::chrono::steady_clock::now();
    double matvec_total_ms = 0.0;
    std::size_t matvec_calls = 0;
    const auto finish = [&]() {
        report.solve_ms = elapsed_ms(clock_start);
        report.matvec_ms = matvec_total_ms;
        report.matvec_calls = matvec_calls;
        return report;
    };

    // Workspace is allocated once. The earliest loop built a fresh vector for the
    // matrix product on every iteration and swept the vectors seven times per
    // step; both show up as pure memory traffic at this problem size, so the
    // products write into a reused buffer and the level-1 operations are fused.
    std::vector<double> preconditioned(n, 0.0);
    std::vector<double> direction(n, 0.0);
    std::vector<double> matrix_direction(n, 0.0);

    // x0 = 0 makes r0 = b, so the initial norm, the preconditioned residual and
    // the first search direction all come out of one preconditioner call.
    preconditioner.apply(residual, preconditioned);
    double residual_norm_squared = 0.0;
    double rho = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double value = residual[i];
        residual_norm_squared += value * value;
        direction[i] = preconditioned[i];
        rho += value * preconditioned[i];
    }
    report.initial_residual = std::sqrt(std::max(0.0, residual_norm_squared));
    const double target = tolerance * std::max(1.0, report.initial_residual);
    if (report.initial_residual <= target) {
        report.converged = true;
        report.final_residual = report.initial_residual;
        report.message = "initial residual satisfies tolerance";
        return finish();
    }
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("PCG failed to form a positive initial preconditioned residual");
    }

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const auto matvec_start = std::chrono::steady_clock::now();
        matrix.multiply_into(direction, matrix_direction);
        matvec_total_ms += elapsed_ms(matvec_start);
        ++matvec_calls;

        double curvature = 0.0;
        for (std::size_t i = 0; i < n; ++i) curvature += direction[i] * matrix_direction[i];
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            throw std::runtime_error("PCG encountered non-positive curvature; system is not SPD");
        }
        const double alpha = rho / curvature;

        // A single pass advances the solution and the residual, then forms the
        // next residual norm, which is otherwise a separate sweep.
        residual_norm_squared = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            solution[i] += alpha * direction[i];
            const double value = residual[i] - alpha * matrix_direction[i];
            residual[i] = value;
            residual_norm_squared += value * value;
        }
        report.iterations = iteration;
        report.final_residual = std::sqrt(std::max(0.0, residual_norm_squared));
        if (report.final_residual <= target) {
            report.converged = true;
            report.message = "converged";
            return finish();
        }

        preconditioner.apply(residual, preconditioned);
        double next_rho = 0.0;
        for (std::size_t i = 0; i < n; ++i) next_rho += residual[i] * preconditioned[i];
        if (!(next_rho > 0.0) || !std::isfinite(next_rho)) {
            throw std::runtime_error("PCG preconditioned residual became non-positive or non-finite");
        }
        const double beta = next_rho / rho;
        for (std::size_t i = 0; i < n; ++i) direction[i] = preconditioned[i] + beta * direction[i];
        rho = next_rho;
    }

    report.final_residual = std::sqrt(std::max(0.0, residual_norm_squared));
    report.message = "maximum iterations reached without convergence";
    return finish();
}

// Maps every scalar unknown to its subdomain. Subdomains are unions of complete
// nodes, so a 2x2 storage block never straddles two subdomains.
std::vector<std::uint32_t> build_dof_subdomain(const Model& model, const std::string& split,
                                               const std::size_t nodes_per_subdomain,
                                               std::size_t& subdomain_count) {
    const std::size_t node_count = model.nodes.size();
    std::vector<double> x(node_count, 0.0);
    std::vector<double> y(node_count, 0.0);
    for (const Node& node : model.nodes) {
        const std::size_t index = node.dof_x / 2;
        if (index >= node_count) continue;
        x[index] = node.position.x;
        y[index] = node.position.y;
    }
    SubdomainSplit split_result = split == "chunk"
                                      ? split_by_chunks(node_count, nodes_per_subdomain)
                                      : split_by_coordinates(x, y, nodes_per_subdomain);
    subdomain_count = split_result.count;
    std::vector<std::uint32_t> dof_subdomain(2 * node_count, 0);
    for (std::size_t node = 0; node < node_count; ++node) {
        const std::uint32_t owner = split_result.node_subdomain[node];
        dof_subdomain[2 * node] = owner;
        dof_subdomain[2 * node + 1] = owner;
    }
    return dof_subdomain;
}

std::vector<double> inverse_diagonal_of(const std::vector<double>& diagonal) {
    std::vector<double> inverse(diagonal.size(), 0.0);
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
        if (!std::isfinite(diagonal[i])) {
            throw std::runtime_error("Jacobi preconditioner encountered a non-finite diagonal");
        }
        // A truss end node has no transverse stiffness, so its diagonal is exactly
        // zero. Such a DOF is decoupled from every other DOF; with a zero residual
        // it never moves, so a unit scaling keeps PCG well defined instead of
        // aborting the whole analysis.
        inverse[i] = std::abs(diagonal[i]) > 1.0e-30 ? 1.0 / diagonal[i] : 1.0;
    }
    return inverse;
}

// The four combinations of storage and preconditioner. Each is a separate
// instantiation of the same Krylov loop, so the only difference between runs is
// how many bytes a product moves and how many iterations the preconditioner
// needs.
LinearSolverReport solve_iterative(const Model& model, const CsrMatrix& matrix,
                                   const std::vector<double>& rhs, std::vector<double>& solution) {
    const SolverSettings& settings = model.settings;
    if (legacy_stage("LIGHTFEM_LEGACY_PCG")) {
        LinearSolverReport report = pcg_jacobi_legacy(matrix, rhs, settings.tolerance,
                                                      settings.max_iterations, solution);
        report.storage = "csr";
        report.preconditioner = "jacobi";
        report.matrix_bytes_per_matvec =
            matrix.nonzeros() * (sizeof(double) + sizeof(std::size_t));
        return report;
    }

    const bool use_bsr = settings.storage == "bsr";
    BsrMatrix bsr;
    double convert_ms = 0.0;
    if (use_bsr) {
        const auto convert_start = std::chrono::steady_clock::now();
        bsr = BsrMatrix::from_csr(matrix);
        convert_ms = elapsed_ms(convert_start);
    }

    SubdomainBlockJacobi block_preconditioner;
    if (settings.preconditioner == "block") {
        std::size_t subdomain_count = 0;
        const std::vector<std::uint32_t> dof_subdomain =
            build_dof_subdomain(model, settings.subdomain_split, settings.subdomain_nodes,
                                subdomain_count);
        block_preconditioner = SubdomainBlockJacobi(matrix, dof_subdomain, subdomain_count);
    }
    std::optional<IncompleteCholesky0> incomplete_cholesky;
    if (settings.preconditioner == "ic0") {
        incomplete_cholesky.emplace(matrix);
    }

    LinearSolverReport report;
    if (settings.preconditioner == "ic0") {
        const TimedPreconditioner<IncompleteCholesky0> preconditioner(*incomplete_cholesky);
        report = use_bsr ? pcg_solve(bsr, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution)
                         : pcg_solve(matrix, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution);
        const IncompleteCholeskyStats stats = incomplete_cholesky->stats();
        report.preconditioner_setup_ms = stats.setup_ms;
        report.preconditioner_apply_ms = stats.apply_ms;
    } else if (settings.preconditioner == "block") {
        const TimedPreconditioner<SubdomainBlockJacobi> preconditioner(block_preconditioner);
        report = use_bsr ? pcg_solve(bsr, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution)
                         : pcg_solve(matrix, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution);
        const BlockJacobiStats stats = block_preconditioner.stats();
        report.subdomains = stats.subdomains;
        report.capture_ratio = stats.capture_ratio;
        report.preconditioner_setup_ms = stats.setup_ms;
        report.preconditioner_apply_ms = stats.apply_ms;
    } else {
        const std::vector<double> diagonal = matrix.diagonal();
        const JacobiPreconditioner preconditioner(inverse_diagonal_of(diagonal));
        report = use_bsr ? pcg_solve(bsr, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution)
                         : pcg_solve(matrix, preconditioner, rhs, settings.tolerance,
                                     settings.max_iterations, solution);
    }

    report.storage = settings.storage;
    report.preconditioner = settings.preconditioner;
    report.storage_convert_ms = convert_ms;
    report.matrix_bytes_per_matvec =
        use_bsr ? bsr.block_nonzeros() * BsrMatrix::block_bytes()
                : matrix.nonzeros() * (sizeof(double) + sizeof(std::size_t));
    return report;
}

// A bar carries a single axial stress. To feed the existing nodal-stress
// smoother (which expects a plane stress tensor per element) the axial stress is
// written back as the rank-one tensor sigma * n x n.
std::array<double, 3> axial_stress_tensor(double sigma, double c, double s) {
    return {sigma * c * c, sigma * s * s, sigma * c * s};
}

std::vector<ElementResult> postprocess(const Model& model, const std::vector<double>& displacement) {
    std::vector<ElementResult> results;
    results.reserve(model.elements.size() + model.trusses.size());
    for (const Q4Element& element : model.elements) {
        const auto dofs = q4_dof_indices(model, element);
        std::array<double, 8> local_displacement{};
        for (std::size_t i = 0; i < 8; ++i) local_displacement[i] = displacement[dofs[i]];
        const Material* material = model.find_material(element.material_id);
        if (material == nullptr) throw std::logic_error("validated element references an unknown material");
        const auto coordinates = q4_coordinates(model, element);
        results.push_back(ElementResult{element.id,
                                        q4_center_field(coordinates, *material, model.settings.plane_mode,
                                                        local_displacement),
                                        q4_gauss_point_fields(coordinates, *material,
                                                              model.settings.plane_mode, local_displacement)});
    }

    for (const TrussElement& truss : model.trusses) {
        const auto dofs = truss_dof_indices(model, truss);
        std::array<double, 4> local_displacement{};
        for (std::size_t i = 0; i < 4; ++i) local_displacement[i] = displacement[dofs[i]];
        const Material* material = model.find_material(truss.material_id);
        if (material == nullptr) throw std::logic_error("validated truss references an unknown material");
        std::array<Vec2, 2> coordinates{};
        for (std::size_t i = 0; i < 2; ++i) {
            const Node* node = model.find_node(truss.node_ids[i]);
            if (node == nullptr) throw std::logic_error("validated truss references an unknown node");
            coordinates[i] = node->position;
        }
        const double dx = coordinates[1].x - coordinates[0].x;
        const double dy = coordinates[1].y - coordinates[0].y;
        const double length = std::sqrt(dx * dx + dy * dy);
        const double c = length > 0.0 ? dx / length : 1.0;
        const double s = length > 0.0 ? dy / length : 0.0;
        ElementField field{};
        field.strain = {truss_axial_strain(coordinates, local_displacement), 0.0, 0.0};
        field.stress = axial_stress_tensor(
            truss_axial_stress(coordinates, material->young_modulus, local_displacement), c, s);
        field.von_mises = std::abs(truss_axial_stress(coordinates, material->young_modulus,
                                                      local_displacement));
        ElementResult result;
        result.element_id = truss.id;
        result.field = field;
        for (auto& gp : result.gauss_points) {
            gp.field = field;  // a bar is uniaxial: every Gauss point carries the same state
        }
        results.push_back(result);
    }
    return results;
}

std::vector<NodeResult> smooth_nodal_stresses(const Model& model,
                                              const std::vector<ElementResult>& element_results) {
    const std::size_t total_elements = model.elements.size() + model.trusses.size();
    if (element_results.size() != total_elements) {
        throw std::logic_error("element result count does not match mesh during nodal stress smoothing");
    }
    constexpr double root_three = 1.7320508075688772935;
    const double alpha = 0.5 * (1.0 + root_three);
    const double beta = 0.5 * (1.0 - root_three);
    const std::array<std::array<double, 4>, 4> extrapolation{{
        {{alpha * alpha, alpha * beta, beta * beta, alpha * beta}},
        {{alpha * beta, alpha * alpha, alpha * beta, beta * beta}},
        {{beta * beta, alpha * beta, alpha * alpha, alpha * beta}},
        {{alpha * beta, beta * beta, alpha * beta, alpha * alpha}},
    }};
    std::vector<std::array<double, 3>> sums(model.nodes.size());
    std::vector<std::size_t> counts(model.nodes.size(), 0);
    for (std::size_t e = 0; e < model.elements.size(); ++e) {
        const Q4Element& element = model.elements[e];
        for (std::size_t local_node = 0; local_node < 4; ++local_node) {
            const Node* node = model.find_node(element.node_ids[local_node]);
            if (node == nullptr) throw std::logic_error("validated element references unknown node");
            const std::size_t node_index = node->dof_x / 2;
            for (std::size_t component = 0; component < 3; ++component) {
                for (std::size_t gp = 0; gp < 4; ++gp) {
                    sums[node_index][component] += extrapolation[local_node][gp] *
                                                   element_results[e].gauss_points[gp].field.stress[component];
                }
            }
            ++counts[node_index];
        }
    }
    for (std::size_t t = 0; t < model.trusses.size(); ++t) {
        const TrussElement& truss = model.trusses[t];
        const ElementField& field = element_results[model.elements.size() + t].field;
        for (const int node_id : truss.node_ids) {
            const Node* node = model.find_node(node_id);
            if (node == nullptr) throw std::logic_error("validated truss references an unknown node");
            const std::size_t node_index = node->dof_x / 2;
            for (std::size_t component = 0; component < 3; ++component) {
                sums[node_index][component] += field.stress[component];
            }
            ++counts[node_index];
        }
    }
    std::vector<NodeResult> results(model.nodes.size());
    for (std::size_t node_index = 0; node_index < model.nodes.size(); ++node_index) {
        if (counts[node_index] == 0) {
            // A node that no element touches (e.g. a bar end left floating in space)
            // carries no stress; keep the export usable instead of aborting.
            results[node_index].node_id = model.nodes[node_index].id;
            continue;
        }
        NodeResult& result = results[node_index];
        result.node_id = model.nodes[node_index].id;
        for (std::size_t component = 0; component < 3; ++component) {
            result.stress[component] = sums[node_index][component] / static_cast<double>(counts[node_index]);
        }
        const double sx = result.stress[0];
        const double sy = result.stress[1];
        const double txy = result.stress[2];
        result.von_mises = std::sqrt(std::max(0.0, sx * sx - sx * sy + sy * sy + 3.0 * txy * txy));
    }
    return results;
}

}  // namespace

AnalysisResult LinearElasticSolver::solve(Model& model) const {
    if (!model.is_indexed()) {
        throw std::logic_error("model must pass validation and sort_and_index before solving");
    }
    const auto solve_start = std::chrono::steady_clock::now();
    AnalysisResult result;
    result.performance.reordering = model.settings.reordering;
    // The "before" profile only means something when the mesh is about to be
    // reordered. With reordering=none it would rebuild the whole CSR pattern
    // for nothing, which on a 400x400 mesh is over a second of pure overhead,
    // so the field is simply mirrored from the post-assembly profile instead.
    const bool reorder_requested = model.settings.reordering == "rcm";
    const auto profile_start = std::chrono::steady_clock::now();
    if (reorder_requested) result.performance.before = profile_for_model(model);
    result.performance.profile_ms = elapsed_ms(profile_start);
    if (reorder_requested) {
        model.index_with_node_order(reverse_cuthill_mckee(model));
    }

    AssembledSystem assembled = assemble(model);
    const auto after_profile_start = std::chrono::steady_clock::now();
    result.performance.after = csr_profile(assembled.stiffness);
    result.performance.profile_ms += elapsed_ms(after_profile_start);
    if (!reorder_requested) result.performance.before = result.performance.after;
    result.performance.assembly = model.settings.assembly;
    result.performance.requested_threads = model.settings.threads;
    result.performance.threads_used = assembled.threads_used;
    result.performance.coloring_count = assembled.coloring_count;
    result.performance.pattern_ms = assembled.pattern_ms;
    result.performance.element_stiffness_ms = assembled.element_stiffness_ms;
    result.performance.assembly_ms = assembled.assembly_ms;
    const auto retain_start = std::chrono::steady_clock::now();
    result.stiffness = assembled.stiffness;
    result.load = assembled.load;
    result.performance.retain_matrix_ms = elapsed_ms(retain_start);

    CsrMatrix constrained_stiffness = std::move(assembled.stiffness);
    std::vector<double> constrained_load = assembled.load;
    const auto boundary_start = std::chrono::steady_clock::now();
    apply_dirichlet(constrained_stiffness, constrained_load, collect_constraints(model));

    // A bar has no stiffness across its axis, so a bar end that no other element
    // touches leaves a zero pivot. If that DOF carries no load it is simply a
    // dangling direction: pin it so the system stays solvable. If it does carry
    // load, the model really is a mechanism and the user must be told.
    {
        const std::vector<double> diagonal = constrained_stiffness.diagonal();
        double load_scale = 0.0;
        for (const double value : constrained_load) load_scale = std::max(load_scale, std::abs(value));
        if (load_scale <= 0.0) load_scale = 1.0;
        for (std::size_t i = 0; i < diagonal.size(); ++i) {
            if (std::abs(diagonal[i]) > 1.0e-30) continue;
            if (std::abs(constrained_load[i]) > 1.0e-12 * load_scale) {
                throw std::runtime_error(
                    "the structure is a mechanism: a bar end has no transverse stiffness but "
                    "carries a transverse load; add a support there or load the bar along its axis [dof "
                    + std::to_string(i) + " load " + std::to_string(constrained_load[i])
                    + " diag " + std::to_string(diagonal[i]) + " node "
                    + std::to_string(model.nodes[i / 2].id) + " xy "
                    + std::to_string(model.nodes[i / 2].position.x) + ","
                    + std::to_string(model.nodes[i / 2].position.y) + "]");
            }
            constrained_stiffness.set(i, i, 1.0);
            constrained_load[i] = 0.0;
        }
    }
    // Only the Dirichlet reduction and the zero-pivot sweep belong to "boundary".
    // Storage conversion and preconditioner setup happen inside the next call and
    // are accounted for by their own fields in the solver report, so stopping the
    // clock here keeps the budget additive instead of double counting.
    result.performance.boundary_ms = elapsed_ms(boundary_start);
    result.solver = solve_iterative(model, constrained_stiffness, constrained_load,
                                    result.displacement);
    if (!result.solver.converged) {
        throw std::runtime_error("PCG did not converge: " + result.solver.message);
    }

    const auto postprocess_start = std::chrono::steady_clock::now();
    result.reaction = result.stiffness.multiply(result.displacement);
    for (std::size_t i = 0; i < result.reaction.size(); ++i) result.reaction[i] -= result.load[i];
    result.element_results = postprocess(model, result.displacement);
    result.node_results = smooth_nodal_stresses(model, result.element_results);
    result.performance.postprocess_ms = elapsed_ms(postprocess_start);
    result.performance.total_ms = elapsed_ms(solve_start);
    return result;
}

}  // namespace lightfem
