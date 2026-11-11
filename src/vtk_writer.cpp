#include "lightfem/vtk_writer.hpp"

#include "lightfem/bench_mode.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <system_error>

namespace lightfem {

namespace {

void require_size(const std::vector<double>& values, const std::size_t expected, const char* label) {
    if (values.size() != expected) {
        throw std::invalid_argument(std::string("VTK ") + label + " size does not match model DOFs");
    }
}

// Result export writes a few million numbers, so the formatting path dominates.
// Two changes keep it cheap: numbers go through std::to_chars (shortest exact
// round-trip form, roughly ten times faster than std::ostream with
// setprecision(17)), and the text is collected in a one-megabyte buffer that is
// flushed in bulk instead of being pushed through the stream one operator<< at
// a time.
class TextFile {
public:
    explicit TextFile(const std::string& path)
        : path_(path),
          output_(path, std::ios::binary),
          legacy_(legacy_stage("LIGHTFEM_LEGACY_EXPORT")) {
        if (!output_) throw std::runtime_error("cannot open output file: " + path);
        if (legacy_) output_ << std::setprecision(17);
        buffer_.reserve(kFlushSize + 4096);
    }

    void text(const char* value) {
        if (legacy_) {
            output_ << value;
            return;
        }
        buffer_ += value;
        flush_if_full();
    }
    void character(const char value) {
        if (legacy_) {
            output_ << value;
            return;
        }
        buffer_ += value;
        flush_if_full();
    }
    void integer(const std::size_t value) { write_integer(static_cast<long long>(value)); }
    void integer(const int value) { write_integer(value); }
    void integer(const long long value) { write_integer(value); }
    void real(const double value) {
        if (legacy_) {
            output_ << value;
            return;
        }
        char text[40];
        const auto result = std::to_chars(text, text + sizeof(text), value);
        if (result.ec != std::errc{}) throw std::runtime_error("failed to format a number for " + path_);
        buffer_.append(text, result.ptr);
        flush_if_full();
    }

    void finish() {
        flush();
        if (!output_) throw std::runtime_error("failed while writing " + path_);
        output_.close();
        if (!output_) throw std::runtime_error("failed while closing " + path_);
    }

private:
    static constexpr std::size_t kFlushSize = std::size_t{1} << 20;

    void write_integer(const long long value) {
        if (legacy_) {
            output_ << value;
            return;
        }
        char text[24];
        const auto result = std::to_chars(text, text + sizeof(text), value);
        if (result.ec != std::errc{}) throw std::runtime_error("failed to format an integer for " + path_);
        buffer_.append(text, result.ptr);
        flush_if_full();
    }
    void flush_if_full() {
        if (buffer_.size() >= kFlushSize) flush();
    }
    void flush() {
        if (buffer_.empty()) return;
        output_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        buffer_.clear();
    }

    std::string path_;
    std::ofstream output_;
    bool legacy_;
    std::string buffer_;
};

}  // namespace

void write_vtk(const Model& model, const AnalysisResult& result, const std::string& path) {
    if (!model.is_indexed()) throw std::logic_error("model must be indexed before VTK export");
    const std::size_t dofs = model.summary().dofs;
    require_size(result.displacement, dofs, "displacement");
    require_size(result.reaction, dofs, "reaction");
    const std::size_t total_cells = model.elements.size() + model.trusses.size();
    if (result.element_results.size() != total_cells) {
        throw std::invalid_argument("VTK element result count does not match model elements");
    }
    if (result.node_results.size() != model.nodes.size()) {
        throw std::invalid_argument("VTK nodal result count does not match model nodes");
    }

    TextFile output(path);
    output.text("# vtk DataFile Version 3.0\n");
    output.text("LightFEM linear static Q4 result\n");
    output.text("ASCII\n");
    output.text("DATASET UNSTRUCTURED_GRID\n");
    output.text("POINTS ");
    output.integer(model.nodes.size());
    output.text(" double\n");
    for (const Node& node : model.nodes) {
        output.real(node.position.x + model.settings.deformation_scale * result.displacement[node.dof_x]);
        output.character(' ');
        output.real(node.position.y + model.settings.deformation_scale * result.displacement[node.dof_y]);
        output.text(" 0\n");
    }

    output.text("CELLS ");
    output.integer(total_cells);
    output.character(' ');
    output.integer(5 * model.elements.size() + 3 * model.trusses.size());
    output.character('\n');
    for (const Q4Element& element : model.elements) {
        output.character('4');
        for (const int node_id : element.node_ids) {
            const Node* node = model.find_node(node_id);
            if (node == nullptr) throw std::logic_error("validated element references an unknown node");
            output.character(' ');
            output.integer(static_cast<std::size_t>(node->dof_x / 2));
        }
        output.character('\n');
    }
    for (const TrussElement& truss : model.trusses) {
        output.character('2');
        for (const int node_id : truss.node_ids) {
            const Node* node = model.find_node(node_id);
            if (node == nullptr) throw std::logic_error("validated truss references an unknown node");
            output.character(' ');
            output.integer(static_cast<std::size_t>(node->dof_x / 2));
        }
        output.character('\n');
    }
    output.text("CELL_TYPES ");
    output.integer(total_cells);
    output.character('\n');
    for (std::size_t i = 0; i < model.elements.size(); ++i) output.text("9\n");  // VTK_QUAD
    for (std::size_t i = 0; i < model.trusses.size(); ++i) output.text("3\n");    // VTK_LINE

    output.text("POINT_DATA ");
    output.integer(model.nodes.size());
    output.text("\nVECTORS displacement double\n");
    for (const Node& node : model.nodes) {
        output.real(result.displacement[node.dof_x]);
        output.character(' ');
        output.real(result.displacement[node.dof_y]);
        output.text(" 0\n");
    }
    const std::array<const char*, 3> nodal_stress_names{{"nodal_stress_xx", "nodal_stress_yy", "nodal_stress_xy"}};
    for (std::size_t component = 0; component < 3; ++component) {
        output.text("SCALARS ");
        output.text(nodal_stress_names[component]);
        output.text(" double 1\nLOOKUP_TABLE default\n");
        for (const NodeResult& node : result.node_results) {
            output.real(node.stress[component]);
            output.character('\n');
        }
    }
    output.text("SCALARS nodal_von_mises_in_plane double 1\nLOOKUP_TABLE default\n");
    for (const NodeResult& node : result.node_results) {
        output.real(node.von_mises);
        output.character('\n');
    }
    output.text("VECTORS reaction double\n");
    for (const Node& node : model.nodes) {
        output.real(result.reaction[node.dof_x]);
        output.character(' ');
        output.real(result.reaction[node.dof_y]);
        output.text(" 0\n");
    }
    output.text("SCALARS displacement_magnitude double 1\nLOOKUP_TABLE default\n");
    for (const Node& node : model.nodes) {
        output.real(std::hypot(result.displacement[node.dof_x], result.displacement[node.dof_y]));
        output.character('\n');
    }

    output.text("CELL_DATA ");
    output.integer(total_cells);
    output.text("\nSCALARS element_id int 1\nLOOKUP_TABLE default\n");
    for (const ElementResult& element : result.element_results) {
        output.integer(element.element_id);
        output.character('\n');
    }

    const std::array<const char*, 3> strain_names{{"strain_xx", "strain_yy", "strain_xy"}};
    const std::array<const char*, 3> stress_names{{"stress_xx", "stress_yy", "stress_xy"}};
    for (std::size_t component = 0; component < 3; ++component) {
        output.text("SCALARS ");
        output.text(strain_names[component]);
        output.text(" double 1\nLOOKUP_TABLE default\n");
        for (const ElementResult& element : result.element_results) {
            output.real(element.field.strain[component]);
            output.character('\n');
        }
    }
    for (std::size_t component = 0; component < 3; ++component) {
        output.text("SCALARS ");
        output.text(stress_names[component]);
        output.text(" double 1\nLOOKUP_TABLE default\n");
        for (const ElementResult& element : result.element_results) {
            output.real(element.field.stress[component]);
            output.character('\n');
        }
    }
    output.text("SCALARS von_mises_in_plane double 1\nLOOKUP_TABLE default\n");
    for (const ElementResult& element : result.element_results) {
        output.real(element.field.von_mises);
        output.character('\n');
    }

    for (std::size_t gp = 0; gp < 4; ++gp) {
        for (std::size_t component = 0; component < 3; ++component) {
            output.text("SCALARS ");
            output.text(stress_names[component]);
            output.text("_gp");
            output.integer(gp + 1);
            output.text(" double 1\nLOOKUP_TABLE default\n");
            for (const ElementResult& element : result.element_results) {
                output.real(element.gauss_points[gp].field.stress[component]);
                output.character('\n');
            }
        }
        output.text("SCALARS von_mises_gp");
        output.integer(gp + 1);
        output.text(" double 1\nLOOKUP_TABLE default\n");
        for (const ElementResult& element : result.element_results) {
            output.real(element.gauss_points[gp].field.von_mises);
            output.character('\n');
        }
    }

    output.finish();
}

CsvExportPaths write_csv_results(const Model& model, const AnalysisResult& result, const std::string& prefix) {
    if (!model.is_indexed()) throw std::logic_error("model must be indexed before CSV export");
    if (prefix.empty()) throw std::invalid_argument("CSV prefix must not be empty");
    const std::size_t dofs = model.summary().dofs;
    require_size(result.displacement, dofs, "displacement");
    require_size(result.reaction, dofs, "reaction");
    if (result.node_results.size() != model.nodes.size() ||
        result.element_results.size() != model.elements.size() + model.trusses.size()) {
        throw std::invalid_argument("CSV result count does not match model");
    }
    CsvExportPaths paths{prefix + "_nodes.csv", prefix + "_gauss_points.csv"};

    TextFile nodes(paths.nodes);
    nodes.text(
        "node_id,x,y,deformed_x,deformed_y,ux,uy,reaction_x,reaction_y,stress_xx,stress_yy,stress_xy,von_mises\n");
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        const Node& node = model.nodes[i];
        const NodeResult& field = result.node_results[i];
        nodes.integer(node.id);
        nodes.character(',');
        nodes.real(node.position.x);
        nodes.character(',');
        nodes.real(node.position.y);
        nodes.character(',');
        nodes.real(node.position.x + model.settings.deformation_scale * result.displacement[node.dof_x]);
        nodes.character(',');
        nodes.real(node.position.y + model.settings.deformation_scale * result.displacement[node.dof_y]);
        nodes.character(',');
        nodes.real(result.displacement[node.dof_x]);
        nodes.character(',');
        nodes.real(result.displacement[node.dof_y]);
        nodes.character(',');
        nodes.real(result.reaction[node.dof_x]);
        nodes.character(',');
        nodes.real(result.reaction[node.dof_y]);
        for (std::size_t component = 0; component < 3; ++component) {
            nodes.character(',');
            nodes.real(field.stress[component]);
        }
        nodes.character(',');
        nodes.real(field.von_mises);
        nodes.character('\n');
    }
    nodes.finish();

    TextFile gauss_points(paths.gauss_points);
    gauss_points.text(
        "element_id,gauss_point,xi,eta,strain_xx,strain_yy,strain_xy,stress_xx,stress_yy,stress_xy,von_mises\n");
    for (const ElementResult& element : result.element_results) {
        for (std::size_t gp = 0; gp < 4; ++gp) {
            const GaussPointField& point = element.gauss_points[gp];
            gauss_points.integer(element.element_id);
            gauss_points.character(',');
            gauss_points.integer(static_cast<std::size_t>(gp + 1));
            gauss_points.character(',');
            gauss_points.real(point.xi);
            gauss_points.character(',');
            gauss_points.real(point.eta);
            for (std::size_t component = 0; component < 3; ++component) {
                gauss_points.character(',');
                gauss_points.real(point.field.strain[component]);
            }
            for (std::size_t component = 0; component < 3; ++component) {
                gauss_points.character(',');
                gauss_points.real(point.field.stress[component]);
            }
            gauss_points.character(',');
            gauss_points.real(point.field.von_mises);
            gauss_points.character('\n');
        }
    }
    gauss_points.finish();
    return paths;
}

}  // namespace lightfem
