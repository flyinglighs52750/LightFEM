#include "lightfem/solver.hpp"
#include "lightfem/reordering.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
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

std::array<Vec2, 4> q4_coordinates(const Model& model, const Q4Element& element) {
    std::array<Vec2, 4> coordinates{};
    for (std::size_t i = 0; i < 4; ++i) {
        const Node* node = model.find_node(element.node_ids[i]);
        if (node == nullptr) throw std::logic_error("validated element references an unknown node");
        coordinates[i] = node->position;
    }
    return coordinates;
}

std::vector<std::array<std::size_t, 8>> element_dof_lists(const Model& model) {
    std::vector<std::array<std::size_t, 8>> element_dofs;
    element_dofs.reserve(model.elements.size());
    for (const Q4Element& element : model.elements) element_dofs.push_back(q4_dof_indices(model, element));
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
    std::unordered_map<int, std::vector<std::size_t>> node_elements;
    node_elements.reserve(model.nodes.size());
    for (std::size_t element = 0; element < model.elements.size(); ++element) {
        for (const int node_id : model.elements[element].node_ids) node_elements[node_id].push_back(element);
    }

    std::vector<int> colors(model.elements.size(), -1);
    int maximum_color = -1;
    std::vector<std::size_t> unavailable_stamp(model.elements.size(), 0);
    std::size_t stamp = 0;
    for (std::size_t element = 0; element < model.elements.size(); ++element) {
        ++stamp;
        for (const int node_id : model.elements[element].node_ids) {
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
                     const std::vector<ElementMatrix>& local_stiffness) {
    for (std::size_t element = 0; element < local_stiffness.size(); ++element) {
        const auto& mapping = pattern.element_value_indices[element];
        for (std::size_t index = 0; index < 64; ++index) stiffness.values()[mapping[index]] += local_stiffness[element][index];
    }
}

void atomic_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                     const std::vector<ElementMatrix>& local_stiffness) {
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
    for (long long element = 0; element < static_cast<long long>(local_stiffness.size()); ++element) {
        const auto& mapping = pattern.element_value_indices[static_cast<std::size_t>(element)];
        for (std::size_t index = 0; index < 64; ++index) {
#pragma omp atomic
            stiffness.values()[mapping[index]] += local_stiffness[static_cast<std::size_t>(element)][index];
        }
    }
#else
    (void)stiffness;
    (void)pattern;
    (void)local_stiffness;
#endif
}

void coloring_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                       const std::vector<ElementMatrix>& local_stiffness,
                       const std::vector<std::vector<std::size_t>>& colors) {
    for (const auto& group : colors) {
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long long local = 0; local < static_cast<long long>(group.size()); ++local) {
            const std::size_t element = group[static_cast<std::size_t>(local)];
            const auto& mapping = pattern.element_value_indices[element];
            for (std::size_t index = 0; index < 64; ++index) {
                stiffness.values()[mapping[index]] += local_stiffness[element][index];
            }
        }
    }
}

void thread_local_assembly(CsrMatrix& stiffness, const CsrAssemblyPattern& pattern,
                           const std::vector<ElementMatrix>& local_stiffness,
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
            for (std::size_t index = 0; index < 64; ++index) {
                local_values[mapping[index]] += local_stiffness[element_index][index];
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
    const std::vector<std::array<std::size_t, 8>> element_dofs = element_dof_lists(model);
    CsrAssemblyPattern pattern = build_csr_assembly_pattern(dofs, element_dofs);

    AssembledSystem system;
    system.stiffness = std::move(pattern.matrix);
    system.load.assign(dofs, 0.0);
    system.threads_used = configure_openmp(model);
    std::vector<ElementMatrix> local_stiffness(model.elements.size());
    const auto stiffness_start = std::chrono::steady_clock::now();
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long element_index = 0; element_index < static_cast<long long>(model.elements.size()); ++element_index) {
        const std::size_t e = static_cast<std::size_t>(element_index);
        const Q4Element& element = model.elements[e];
        const Material* material = model.find_material(element.material_id);
        if (material == nullptr) throw std::logic_error("validated element references an unknown material");
        local_stiffness[e] = q4_stiffness(q4_coordinates(model, element), *material, model.settings.plane_mode);
    }
    const auto stiffness_end = std::chrono::steady_clock::now();
    system.element_stiffness_ms =
        std::chrono::duration<double, std::milli>(stiffness_end - stiffness_start).count();

    const auto assembly_start = std::chrono::steady_clock::now();
    if (model.settings.assembly == "serial") {
        serial_assembly(system.stiffness, pattern, local_stiffness);
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
}

LinearSolverReport pcg_jacobi(const CsrMatrix& matrix, const std::vector<double>& rhs,
                               const double tolerance, const std::size_t max_iterations,
                               std::vector<double>& solution) {
    if (rhs.size() != matrix.rows()) throw std::invalid_argument("PCG RHS size mismatch");
    solution.assign(matrix.rows(), 0.0);
    std::vector<double> residual = rhs;
    const std::vector<double> diagonal = matrix.diagonal();
    std::vector<double> inverse_diagonal(matrix.rows(), 0.0);
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        if (!std::isfinite(diagonal[i]) || std::abs(diagonal[i]) <= 1.0e-30) {
            throw std::runtime_error("Jacobi preconditioner encountered a zero or non-finite diagonal");
        }
        inverse_diagonal[i] = 1.0 / diagonal[i];
    }

    LinearSolverReport report;
    report.initial_residual = norm(residual);
    const double target = tolerance * std::max(1.0, report.initial_residual);
    if (report.initial_residual <= target) {
        report.converged = true;
        report.final_residual = report.initial_residual;
        report.message = "initial residual satisfies tolerance";
        return report;
    }

    std::vector<double> preconditioned(matrix.rows(), 0.0);
    for (std::size_t i = 0; i < matrix.rows(); ++i) preconditioned[i] = inverse_diagonal[i] * residual[i];
    std::vector<double> direction = preconditioned;
    double rho = dot(residual, preconditioned);
    if (!(rho > 0.0) || !std::isfinite(rho)) {
        throw std::runtime_error("PCG failed to form a positive initial preconditioned residual");
    }

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const std::vector<double> matrix_direction = matrix.multiply(direction);
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
            return report;
        }

        for (std::size_t i = 0; i < matrix.rows(); ++i) {
            preconditioned[i] = inverse_diagonal[i] * residual[i];
        }
        const double next_rho = dot(residual, preconditioned);
        if (!(next_rho > 0.0) || !std::isfinite(next_rho)) {
            throw std::runtime_error("PCG preconditioned residual became non-positive or non-finite");
        }
        const double beta = next_rho / rho;
        for (std::size_t i = 0; i < matrix.rows(); ++i) {
            direction[i] = preconditioned[i] + beta * direction[i];
        }
        rho = next_rho;
    }

    report.final_residual = norm(residual);
    report.message = "maximum iterations reached without convergence";
    return report;
}

std::vector<ElementResult> postprocess(const Model& model, const std::vector<double>& displacement) {
    std::vector<ElementResult> results;
    results.reserve(model.elements.size());
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
    return results;
}

std::vector<NodeResult> smooth_nodal_stresses(const Model& model,
                                              const std::vector<ElementResult>& element_results) {
    if (element_results.size() != model.elements.size()) {
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
    std::vector<NodeResult> results(model.nodes.size());
    for (std::size_t node_index = 0; node_index < model.nodes.size(); ++node_index) {
        if (counts[node_index] == 0) throw std::logic_error("node is not connected during stress smoothing");
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
    AnalysisResult result;
    result.performance.reordering = model.settings.reordering;
    result.performance.before = profile_for_model(model);
    if (model.settings.reordering == "rcm") {
        model.index_with_node_order(reverse_cuthill_mckee(model));
    }

    AssembledSystem assembled = assemble(model);
    result.performance.after = csr_profile(assembled.stiffness);
    result.performance.assembly = model.settings.assembly;
    result.performance.requested_threads = model.settings.threads;
    result.performance.threads_used = assembled.threads_used;
    result.performance.coloring_count = assembled.coloring_count;
    result.performance.element_stiffness_ms = assembled.element_stiffness_ms;
    result.performance.assembly_ms = assembled.assembly_ms;
    result.stiffness = assembled.stiffness;
    result.load = assembled.load;

    CsrMatrix constrained_stiffness = std::move(assembled.stiffness);
    std::vector<double> constrained_load = assembled.load;
    apply_dirichlet(constrained_stiffness, constrained_load, collect_constraints(model));
    result.solver = pcg_jacobi(constrained_stiffness, constrained_load, model.settings.tolerance,
                               model.settings.max_iterations, result.displacement);
    if (!result.solver.converged) {
        throw std::runtime_error("PCG did not converge: " + result.solver.message);
    }

    result.reaction = result.stiffness.multiply(result.displacement);
    for (std::size_t i = 0; i < result.reaction.size(); ++i) result.reaction[i] -= result.load[i];
    result.element_results = postprocess(model, result.displacement);
    result.node_results = smooth_nodal_stresses(model, result.element_results);
    return result;
}

}  // namespace lightfem
