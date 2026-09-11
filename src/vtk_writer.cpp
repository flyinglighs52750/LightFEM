#include "lightfem/vtk_writer.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lightfem {

namespace {

void require_size(const std::vector<double>& values, const std::size_t expected, const char* label) {
    if (values.size() != expected) {
        throw std::invalid_argument(std::string("VTK ") + label + " size does not match model DOFs");
    }
}

}  // namespace

void write_vtk(const Model& model, const AnalysisResult& result, const std::string& path) {
    if (!model.is_indexed()) throw std::logic_error("model must be indexed before VTK export");
    const std::size_t dofs = model.summary().dofs;
    require_size(result.displacement, dofs, "displacement");
    require_size(result.reaction, dofs, "reaction");
    if (result.element_results.size() != model.elements.size()) {
        throw std::invalid_argument("VTK element result count does not match model elements");
    }
    if (result.node_results.size() != model.nodes.size()) {
        throw std::invalid_argument("VTK nodal result count does not match model nodes");
    }

    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open VTK output file: " + path);
    output << std::setprecision(17);
    output << "# vtk DataFile Version 3.0\n"
           << "LightFEM linear static Q4 result\n"
           << "ASCII\n"
           << "DATASET UNSTRUCTURED_GRID\n"
           << "POINTS " << model.nodes.size() << " double\n";
    for (const Node& node : model.nodes) {
        output << node.position.x + model.settings.deformation_scale * result.displacement[node.dof_x] << ' '
               << node.position.y + model.settings.deformation_scale * result.displacement[node.dof_y] << " 0\n";
    }

    output << "CELLS " << model.elements.size() << ' ' << 5 * model.elements.size() << "\n";
    for (const Q4Element& element : model.elements) {
        output << "4";
        for (const int node_id : element.node_ids) {
            const Node* node = model.find_node(node_id);
            if (node == nullptr) throw std::logic_error("validated element references an unknown node");
            output << ' ' << node->dof_x / 2;
        }
        output << '\n';
    }
    output << "CELL_TYPES " << model.elements.size() << "\n";
    for (std::size_t i = 0; i < model.elements.size(); ++i) output << "9\n";  // VTK_QUAD

    output << "POINT_DATA " << model.nodes.size() << "\n"
           << "VECTORS displacement double\n";
    for (const Node& node : model.nodes) {
        output << result.displacement[node.dof_x] << ' ' << result.displacement[node.dof_y] << " 0\n";
    }
    const std::array<const char*, 3> nodal_stress_names{{"nodal_stress_xx", "nodal_stress_yy", "nodal_stress_xy"}};
    for (std::size_t component = 0; component < 3; ++component) {
        output << "SCALARS " << nodal_stress_names[component] << " double 1\nLOOKUP_TABLE default\n";
        for (const NodeResult& node : result.node_results) output << node.stress[component] << '\n';
    }
    output << "SCALARS nodal_von_mises_in_plane double 1\nLOOKUP_TABLE default\n";
    for (const NodeResult& node : result.node_results) output << node.von_mises << '\n';
    output << "VECTORS reaction double\n";
    for (const Node& node : model.nodes) {
        output << result.reaction[node.dof_x] << ' ' << result.reaction[node.dof_y] << " 0\n";
    }
    output << "SCALARS displacement_magnitude double 1\nLOOKUP_TABLE default\n";
    for (const Node& node : model.nodes) {
        output << std::hypot(result.displacement[node.dof_x], result.displacement[node.dof_y]) << '\n';
    }

    output << "CELL_DATA " << model.elements.size() << "\n"
           << "SCALARS element_id int 1\nLOOKUP_TABLE default\n";
    for (const ElementResult& element : result.element_results) output << element.element_id << '\n';

    const std::array<const char*, 3> strain_names{{"strain_xx", "strain_yy", "strain_xy"}};
    const std::array<const char*, 3> stress_names{{"stress_xx", "stress_yy", "stress_xy"}};
    for (std::size_t component = 0; component < 3; ++component) {
        output << "SCALARS " << strain_names[component] << " double 1\nLOOKUP_TABLE default\n";
        for (const ElementResult& element : result.element_results) {
            output << element.field.strain[component] << '\n';
        }
    }
    for (std::size_t component = 0; component < 3; ++component) {
        output << "SCALARS " << stress_names[component] << " double 1\nLOOKUP_TABLE default\n";
        for (const ElementResult& element : result.element_results) {
            output << element.field.stress[component] << '\n';
        }
    }
    output << "SCALARS von_mises_in_plane double 1\nLOOKUP_TABLE default\n";
    for (const ElementResult& element : result.element_results) output << element.field.von_mises << '\n';

    const std::array<const char*, 3> gauss_stress_names{{"stress_xx", "stress_yy", "stress_xy"}};
    for (std::size_t gp = 0; gp < 4; ++gp) {
        for (std::size_t component = 0; component < 3; ++component) {
            output << "SCALARS " << gauss_stress_names[component] << "_gp" << (gp + 1)
                   << " double 1\nLOOKUP_TABLE default\n";
            for (const ElementResult& element : result.element_results) {
                output << element.gauss_points[gp].field.stress[component] << '\n';
            }
        }
        output << "SCALARS von_mises_gp" << (gp + 1) << " double 1\nLOOKUP_TABLE default\n";
        for (const ElementResult& element : result.element_results) {
            output << element.gauss_points[gp].field.von_mises << '\n';
        }
    }

    if (!output) throw std::runtime_error("failed while writing VTK output file: " + path);
}

CsvExportPaths write_csv_results(const Model& model, const AnalysisResult& result, const std::string& prefix) {
    if (!model.is_indexed()) throw std::logic_error("model must be indexed before CSV export");
    if (prefix.empty()) throw std::invalid_argument("CSV prefix must not be empty");
    const std::size_t dofs = model.summary().dofs;
    require_size(result.displacement, dofs, "displacement");
    require_size(result.reaction, dofs, "reaction");
    if (result.node_results.size() != model.nodes.size() || result.element_results.size() != model.elements.size()) {
        throw std::invalid_argument("CSV result count does not match model");
    }
    CsvExportPaths paths{prefix + "_nodes.csv", prefix + "_gauss_points.csv"};
    std::ofstream nodes(paths.nodes);
    if (!nodes) throw std::runtime_error("cannot open node CSV output: " + paths.nodes);
    nodes << std::setprecision(17)
          << "node_id,x,y,deformed_x,deformed_y,ux,uy,reaction_x,reaction_y,stress_xx,stress_yy,stress_xy,von_mises\n";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const Node& node = model.nodes[i];
        const NodeResult& field = result.node_results[i];
        nodes << node.id << ',' << node.position.x << ',' << node.position.y << ','
              << node.position.x + model.settings.deformation_scale * result.displacement[node.dof_x] << ','
              << node.position.y + model.settings.deformation_scale * result.displacement[node.dof_y] << ','
              << result.displacement[node.dof_x] << ',' << result.displacement[node.dof_y] << ','
              << result.reaction[node.dof_x] << ',' << result.reaction[node.dof_y] << ','
              << field.stress[0] << ',' << field.stress[1] << ',' << field.stress[2] << ',' << field.von_mises
              << '\n';
    }
    std::ofstream gauss_points(paths.gauss_points);
    if (!gauss_points) throw std::runtime_error("cannot open Gauss-point CSV output: " + paths.gauss_points);
    gauss_points << std::setprecision(17)
                 << "element_id,gauss_point,xi,eta,strain_xx,strain_yy,strain_xy,stress_xx,stress_yy,stress_xy,von_mises\n";
    for (const ElementResult& element : result.element_results) {
        for (std::size_t gp = 0; gp < 4; ++gp) {
            const GaussPointField& point = element.gauss_points[gp];
            gauss_points << element.element_id << ',' << (gp + 1) << ',' << point.xi << ',' << point.eta << ','
                         << point.field.strain[0] << ',' << point.field.strain[1] << ',' << point.field.strain[2] << ','
                         << point.field.stress[0] << ',' << point.field.stress[1] << ',' << point.field.stress[2] << ','
                         << point.field.von_mises << '\n';
        }
    }
    if (!nodes || !gauss_points) throw std::runtime_error("failed while writing CSV results");
    return paths;
}

}  // namespace lightfem
