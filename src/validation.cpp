#include "lightfem/validation.hpp"

#include "lightfem/bench_mode.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lightfem {

namespace {

template <typename T>
void validate_ids(const std::vector<T>& entities, const char* entity_name, const char* duplicate_code,
                  ValidationReport& report) {
    std::unordered_set<int> ids;
    for (const auto& entity : entities) {
        if (entity.id <= 0) {
            report.error("E_NONPOSITIVE_ID",
                         std::string(entity_name) + " ID must be positive: " + std::to_string(entity.id),
                         entity.source_line);
        }
        if (!ids.insert(entity.id).second) {
            report.error(duplicate_code,
                         std::string("duplicate ") + entity_name + " ID: " + std::to_string(entity.id),
                         entity.source_line);
        }
    }
}

double q4_det_jacobian(const std::array<const Node*, 4>& nodes, const double xi, const double eta) {
    const std::array<double, 4> dxi{{-(1.0 - eta), +(1.0 - eta), +(1.0 + eta), -(1.0 + eta)}};
    const std::array<double, 4> deta{{-(1.0 - xi), -(1.0 + xi), +(1.0 + xi), +(1.0 - xi)}};
    double dx_dxi = 0.0;
    double dy_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_deta = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        dx_dxi += 0.25 * dxi[i] * nodes[i]->position.x;
        dy_dxi += 0.25 * dxi[i] * nodes[i]->position.y;
        dx_deta += 0.25 * deta[i] * nodes[i]->position.x;
        dy_deta += 0.25 * deta[i] * nodes[i]->position.y;
    }
    return dx_dxi * dy_deta - dx_deta * dy_dxi;
}

int rigid_constraint_rank(std::vector<std::array<double, 3>> rows) {
    if (rows.empty()) return 0;
    int rank = 0;
    constexpr double tolerance = 1.0e-12;
    for (int col = 0; col < 3 && rank < static_cast<int>(rows.size()); ++col) {
        int pivot = rank;
        for (int row = rank + 1; row < static_cast<int>(rows.size()); ++row) {
            if (std::abs(rows[row][col]) > std::abs(rows[pivot][col])) pivot = row;
        }
        if (std::abs(rows[pivot][col]) <= tolerance) continue;
        std::swap(rows[pivot], rows[rank]);
        const double scale = rows[rank][col];
        for (int j = col; j < 3; ++j) rows[rank][j] /= scale;
        for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
            if (row == rank) continue;
            const double factor = rows[row][col];
            for (int j = col; j < 3; ++j) rows[row][j] -= factor * rows[rank][j];
        }
        ++rank;
    }
    return rank;
}

std::string element_nodes_text(const Q4Element& element) {
    std::ostringstream out;
    out << '[' << element.node_ids[0] << ", " << element.node_ids[1] << ", "
        << element.node_ids[2] << ", " << element.node_ids[3] << ']';
    return out.str();
}

}  // namespace

void ValidationReport::error(std::string code, std::string message, const std::size_t source_line) {
    diagnostics_.push_back(Diagnostic{Severity::Error, std::move(code), std::move(message), source_line});
}

void ValidationReport::warning(std::string code, std::string message, const std::size_t source_line) {
    diagnostics_.push_back(Diagnostic{Severity::Warning, std::move(code), std::move(message), source_line});
}

bool ValidationReport::ok() const noexcept { return error_count() == 0; }

std::size_t ValidationReport::error_count() const noexcept {
    return static_cast<std::size_t>(std::count_if(diagnostics_.begin(), diagnostics_.end(),
                                                  [](const Diagnostic& diagnostic) {
                                                      return diagnostic.severity == Severity::Error;
                                                  }));
}

std::size_t ValidationReport::warning_count() const noexcept {
    return static_cast<std::size_t>(std::count_if(diagnostics_.begin(), diagnostics_.end(),
                                                  [](const Diagnostic& diagnostic) {
                                                      return diagnostic.severity == Severity::Warning;
                                                  }));
}

ValidationReport ModelValidator::validate(const Model& model) const {
    ValidationReport report;

    if (model.nodes.size() < 2) report.error("E_TOO_FEW_NODES", "at least two nodes are required");
    if (model.materials.empty()) report.error("E_NO_MATERIAL", "at least one material is required");
    if (model.elements.empty() && model.trusses.empty()) {
        report.error("E_NO_ELEMENT", "at least one Q4 element or truss member is required");
    }

    validate_ids(model.nodes, "node", "E_DUPLICATE_NODE_ID", report);
    validate_ids(model.materials, "material", "E_DUPLICATE_MATERIAL_ID", report);
    validate_ids(model.elements, "element", "E_DUPLICATE_ELEMENT_ID", report);
    validate_ids(model.trusses, "truss", "E_DUPLICATE_TRUSS_ID", report);

    std::unordered_map<int, const Node*> nodes;
    for (const Node& node : model.nodes) {
        nodes.emplace(node.id, &node);
        if (!std::isfinite(node.position.x) || !std::isfinite(node.position.y)) {
            report.error("E_INVALID_COORDINATE", "node coordinates must be finite", node.source_line);
        }
    }

    // Coincident-node detection. Comparing every pair is O(n^2) and dominated the
    // whole run time once meshes reached tens of thousands of nodes (4.5 s of a
    // 10 s run on 40k nodes). Sorting by x first turns it into an O(n log n)
    // sweep: only nodes whose x separation is inside the tolerance can be
    // coincident, so the inner loop advances a short window.
    if (legacy_stage("LIGHTFEM_LEGACY_VALIDATE")) {
        for (std::size_t i = 0; i < model.nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < model.nodes.size(); ++j) {
                if (std::abs(model.nodes[i].position.x - model.nodes[j].position.x) <= 1.0e-14 &&
                    std::abs(model.nodes[i].position.y - model.nodes[j].position.y) <= 1.0e-14) {
                    report.error("E_COINCIDENT_NODES",
                                 "nodes " + std::to_string(model.nodes[i].id) + " and " +
                                     std::to_string(model.nodes[j].id) + " have coincident coordinates",
                                 model.nodes[j].source_line);
                }
            }
        }
    } else {
        std::vector<std::size_t> by_x(model.nodes.size());
        for (std::size_t i = 0; i < by_x.size(); ++i) by_x[i] = i;
        std::sort(by_x.begin(), by_x.end(), [&model](const std::size_t a, const std::size_t b) {
            return model.nodes[a].position.x < model.nodes[b].position.x;
        });
        for (std::size_t a = 0; a < by_x.size(); ++a) {
            const Node& left = model.nodes[by_x[a]];
            for (std::size_t b = a + 1; b < by_x.size(); ++b) {
                const Node& right = model.nodes[by_x[b]];
                if (right.position.x - left.position.x > 1.0e-14) break;
                if (std::abs(left.position.y - right.position.y) <= 1.0e-14) {
                    report.error("E_COINCIDENT_NODES",
                                 "nodes " + std::to_string(left.id) + " and " +
                                     std::to_string(right.id) + " have coincident coordinates",
                                 right.source_line);
                }
            }
        }
    }

    std::unordered_map<int, const Material*> materials;
    for (const Material& material : model.materials) {
        materials.emplace(material.id, &material);
        if (!std::isfinite(material.young_modulus) || material.young_modulus <= 0.0) {
            report.error("E_INVALID_YOUNG_MODULUS", "Young's modulus must be finite and positive",
                         material.source_line);
        }
        if (!std::isfinite(material.poisson_ratio) || material.poisson_ratio <= -1.0 ||
            material.poisson_ratio >= 0.5) {
            report.error("E_INVALID_POISSON_RATIO", "Poisson's ratio must satisfy -1 < nu < 0.5",
                         material.source_line);
        }
        if (!std::isfinite(material.thickness) || material.thickness <= 0.0) {
            report.error("E_INVALID_THICKNESS", "material thickness must be finite and positive",
                         material.source_line);
        }
    }

    std::unordered_set<int> referenced_nodes;
    std::unordered_set<int> used_materials;
    std::unordered_map<int, const Q4Element*> elements;
    std::unordered_map<int, std::vector<int>> adjacency;
    constexpr double gauss = 0.57735026918962576451;
    const std::array<std::array<double, 2>, 5> points{{
        {{0.0, 0.0}}, {{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}},
    }};

    for (const Q4Element& element : model.elements) {
        elements.emplace(element.id, &element);
        const std::set<int> distinct(element.node_ids.begin(), element.node_ids.end());
        if (distinct.size() != 4) {
            report.error("E_REPEATED_ELEMENT_NODE",
                         "Q4 element " + std::to_string(element.id) +
                             " must reference four distinct nodes: " + element_nodes_text(element),
                         element.source_line);
        }

        std::array<const Node*, 4> element_nodes{{nullptr, nullptr, nullptr, nullptr}};
        bool references_ok = true;
        for (std::size_t i = 0; i < 4; ++i) {
            const auto it = nodes.find(element.node_ids[i]);
            if (it == nodes.end()) {
                report.error("E_UNKNOWN_ELEMENT_NODE",
                             "element " + std::to_string(element.id) + " references unknown node " +
                                 std::to_string(element.node_ids[i]),
                             element.source_line);
                references_ok = false;
            } else {
                element_nodes[i] = it->second;
                referenced_nodes.insert(element.node_ids[i]);
            }
        }
        if (materials.find(element.material_id) == materials.end()) {
            report.error("E_UNKNOWN_MATERIAL",
                         "element " + std::to_string(element.id) + " references unknown material " +
                             std::to_string(element.material_id),
                         element.source_line);
        } else {
            used_materials.insert(element.material_id);
        }

        if (references_ok && distinct.size() == 4) {
            for (const auto& point : points) {
                const double det = q4_det_jacobian(element_nodes, point[0], point[1]);
                if (!(det > 1.0e-14)) {
                    report.error("E_NONPOSITIVE_JACOBIAN",
                                 "element " + std::to_string(element.id) +
                                     " has non-positive or near-zero Jacobian; use counter-clockwise, "
                                     "non-self-intersecting Q4 connectivity",
                                 element.source_line);
                    break;
                }
            }
            for (std::size_t i = 0; i < 4; ++i) {
                const int a = element.node_ids[i];
                const int b = element.node_ids[(i + 1) % 4];
                adjacency[a].push_back(b);
                adjacency[b].push_back(a);
            }
        }
    }

    for (const TrussElement& truss : model.trusses) {
        if (truss.node_ids[0] == truss.node_ids[1]) {
            report.error("E_REPEATED_TRUSS_NODE",
                         "truss " + std::to_string(truss.id) +
                             " must reference two distinct nodes: [" +
                             std::to_string(truss.node_ids[0]) + ", " +
                             std::to_string(truss.node_ids[1]) + "]",
                         truss.source_line);
        }
        std::array<const Node*, 2> truss_nodes{{nullptr, nullptr}};
        bool references_ok = true;
        for (std::size_t i = 0; i < 2; ++i) {
            const auto it = nodes.find(truss.node_ids[i]);
            if (it == nodes.end()) {
                report.error("E_UNKNOWN_TRUSS_NODE",
                             "truss " + std::to_string(truss.id) + " references unknown node " +
                                 std::to_string(truss.node_ids[i]),
                             truss.source_line);
                references_ok = false;
            } else {
                truss_nodes[i] = it->second;
                referenced_nodes.insert(truss.node_ids[i]);
            }
        }
        if (materials.find(truss.material_id) == materials.end()) {
            report.error("E_UNKNOWN_MATERIAL",
                         "truss " + std::to_string(truss.id) + " references unknown material " +
                             std::to_string(truss.material_id),
                         truss.source_line);
        } else {
            used_materials.insert(truss.material_id);
        }
        if (!std::isfinite(truss.area) || truss.area <= 0.0) {
            report.error("E_INVALID_TRUSS_AREA", "truss cross-section area must be finite and positive",
                         truss.source_line);
        }
        if (references_ok && truss.node_ids[0] != truss.node_ids[1]) {
            const double dx = truss_nodes[1]->position.x - truss_nodes[0]->position.x;
            const double dy = truss_nodes[1]->position.y - truss_nodes[0]->position.y;
            if (!(std::sqrt(dx * dx + dy * dy) > 1.0e-12)) {
                report.error("E_ZERO_LENGTH_TRUSS",
                             "truss " + std::to_string(truss.id) + " has zero length", truss.source_line);
            } else {
                adjacency[truss.node_ids[0]].push_back(truss.node_ids[1]);
                adjacency[truss.node_ids[1]].push_back(truss.node_ids[0]);
            }
        }
    }

    for (const Node& node : model.nodes) {
        if (referenced_nodes.find(node.id) == referenced_nodes.end()) {
            report.warning("W_UNUSED_NODE", "node " + std::to_string(node.id) + " is not used by any element",
                           node.source_line);
        }
    }
    for (const Material& material : model.materials) {
        if (used_materials.find(material.id) == used_materials.end()) {
            report.warning("W_UNUSED_MATERIAL",
                           "material " + std::to_string(material.id) + " is not used by any element",
                           material.source_line);
        }
    }

    if (!referenced_nodes.empty()) {
        std::unordered_set<int> visited;
        std::size_t components = 0;
        for (const int start : referenced_nodes) {
            if (visited.find(start) != visited.end()) continue;
            ++components;
            std::queue<int> pending;
            pending.push(start);
            visited.insert(start);
            while (!pending.empty()) {
                const int current = pending.front();
                pending.pop();
                for (const int next : adjacency[current]) {
                    if (visited.insert(next).second) pending.push(next);
                }
            }
        }
        if (components > 1) {
            report.error("E_DISCONNECTED_MESH",
                         "the current solver scope requires one connected mesh; found " +
                             std::to_string(components) + " components");
        }
    }

    std::set<std::pair<int, char>> constrained_dofs;
    std::vector<std::array<double, 3>> rigid_rows;
    bool has_nonzero_prescribed_displacement = false;
    for (const DirichletBC& boundary : model.boundaries) {
        const auto node_it = nodes.find(boundary.node_id);
        if (node_it == nodes.end()) {
            report.error("E_UNKNOWN_BOUNDARY_NODE",
                         "boundary condition references unknown node " + std::to_string(boundary.node_id),
                         boundary.source_line);
            continue;
        }
        if (!boundary.ux && !boundary.uy) {
            report.warning("W_EMPTY_BOUNDARY", "boundary record constrains neither ux nor uy",
                           boundary.source_line);
        }
        const Node& node = *node_it->second;
        if (boundary.ux) {
            if (!constrained_dofs.emplace(boundary.node_id, 'x').second) {
                report.error("E_DUPLICATE_CONSTRAINT",
                             "ux is constrained more than once at node " + std::to_string(boundary.node_id),
                             boundary.source_line);
            }
            rigid_rows.push_back({{1.0, 0.0, -node.position.y}});
            has_nonzero_prescribed_displacement |= std::abs(*boundary.ux) > 0.0;
        }
        if (boundary.uy) {
            if (!constrained_dofs.emplace(boundary.node_id, 'y').second) {
                report.error("E_DUPLICATE_CONSTRAINT",
                             "uy is constrained more than once at node " + std::to_string(boundary.node_id),
                             boundary.source_line);
            }
            rigid_rows.push_back({{0.0, 1.0, node.position.x}});
            has_nonzero_prescribed_displacement |= std::abs(*boundary.uy) > 0.0;
        }
    }
    if (rigid_constraint_rank(rigid_rows) < 3) {
        report.error("E_RIGID_BODY_MODES",
                     "displacement constraints do not remove all three 2D rigid-body modes");
    }

    std::unordered_map<int, std::size_t> load_count;
    bool has_nonzero_load = false;
    for (const NodalLoad& load : model.loads) {
        if (nodes.find(load.node_id) == nodes.end()) {
            report.error("E_UNKNOWN_LOAD_NODE", "load references unknown node " + std::to_string(load.node_id),
                         load.source_line);
            continue;
        }
        if (++load_count[load.node_id] > 1) {
            report.warning("W_REPEATED_LOAD_NODE",
                           "multiple load records at node " + std::to_string(load.node_id) +
                               " will be accumulated by the solver",
                           load.source_line);
        }
        has_nonzero_load |= std::abs(load.fx) > 0.0 || std::abs(load.fy) > 0.0;
        if (std::abs(load.fx) > 0.0 && constrained_dofs.count({load.node_id, 'x'}) != 0) {
            report.warning("W_LOAD_ON_CONSTRAINED_DOF",
                           "nonzero fx is applied to constrained ux at node " + std::to_string(load.node_id),
                           load.source_line);
        }
        if (std::abs(load.fy) > 0.0 && constrained_dofs.count({load.node_id, 'y'}) != 0) {
            report.warning("W_LOAD_ON_CONSTRAINED_DOF",
                           "nonzero fy is applied to constrained uy at node " + std::to_string(load.node_id),
                           load.source_line);
        }
    }
    for (const BodyForce& force : model.body_forces) {
        if (elements.find(force.element_id) == elements.end()) {
            report.error("E_UNKNOWN_BODY_FORCE_ELEMENT",
                         "body force references unknown element " + std::to_string(force.element_id),
                         force.source_line);
            continue;
        }
        has_nonzero_load |= std::abs(force.bx) > 0.0 || std::abs(force.by) > 0.0;
    }
    for (const EdgeLoad& load : model.edge_loads) {
        if (elements.find(load.element_id) == elements.end()) {
            report.error("E_UNKNOWN_EDGE_LOAD_ELEMENT",
                         "edge load references unknown element " + std::to_string(load.element_id),
                         load.source_line);
            continue;
        }
        if (load.edge < 1 || load.edge > 4) {
            report.error("E_INVALID_EDGE_LOAD_EDGE", "Q4 edge load edge must be in [1, 4]", load.source_line);
        }
        has_nonzero_load |= std::abs(load.first) > 0.0 || std::abs(load.second) > 0.0;
    }
    if (!has_nonzero_load && !has_nonzero_prescribed_displacement) {
        report.warning("W_NO_EXCITATION", "model has no nonzero load or prescribed displacement");
    }

    if (!std::isfinite(model.settings.tolerance) || model.settings.tolerance <= 0.0 ||
        model.settings.tolerance >= 1.0) {
        report.error("E_INVALID_TOLERANCE", "solver tolerance must satisfy 0 < tolerance < 1");
    }
    if (model.settings.max_iterations == 0) {
        report.error("E_INVALID_MAX_ITERATIONS", "max_iterations must be greater than zero");
    }
    if (!std::isfinite(model.settings.deformation_scale) || model.settings.deformation_scale < 0.0) {
        report.error("E_INVALID_DEFORMATION_SCALE", "deformation_scale must be finite and non-negative");
    }
    if (model.settings.solver != "pcg") {
        report.error("E_UNSUPPORTED_SOLVER", "only solver=pcg is supported in the current project scope");
    }
    if (model.settings.reordering != "none" && model.settings.reordering != "rcm") {
        report.error("E_UNSUPPORTED_REORDERING", "reordering must be none or rcm");
    }
    if (model.settings.assembly != "serial" && model.settings.assembly != "atomic" &&
        model.settings.assembly != "coloring" && model.settings.assembly != "thread_local" &&
        model.settings.assembly != "scan") {
        report.error("E_UNSUPPORTED_ASSEMBLY",
                     "assembly must be serial, atomic, coloring, thread_local, or scan");
    }
    if (model.settings.storage != "csr" && model.settings.storage != "bsr") {
        report.error("E_UNSUPPORTED_STORAGE", "storage must be csr or bsr");
    }
    if (model.settings.preconditioner != "jacobi" && model.settings.preconditioner != "block" &&
        model.settings.preconditioner != "ic0") {
        report.error("E_UNSUPPORTED_PRECONDITIONER", "preconditioner must be jacobi, block, or ic0");
    }
    if (model.settings.subdomain_split != "rcb" && model.settings.subdomain_split != "chunk") {
        report.error("E_UNSUPPORTED_SUBDOMAIN_SPLIT", "subdomain_split must be rcb or chunk");
    }
    if (model.settings.subdomain_nodes == 0) {
        report.error("E_INVALID_SUBDOMAIN_SIZE", "subdomain_nodes must be at least 1");
    }
    if (model.settings.output_file.empty()) {
        report.error("E_EMPTY_OUTPUT", "output file name must not be empty");
    } else if (std::filesystem::path(model.settings.output_file).extension() != ".vtk") {
        report.warning("W_OUTPUT_EXTENSION", "VTK output file normally uses the .vtk extension");
    }
    if (model.settings.csv_prefix.empty() && !model.settings.output_file.empty()) {
        // An empty prefix is intentional: the solver derives it from output.
    }

    return report;
}

const char* to_string(const Severity severity) noexcept {
    return severity == Severity::Error ? "ERROR" : "WARNING";
}

}  // namespace lightfem
