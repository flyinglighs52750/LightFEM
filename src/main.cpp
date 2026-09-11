#include "lightfem/input_parser.hpp"
#include "lightfem/model.hpp"
#include "lightfem/solver.hpp"
#include "lightfem/validation.hpp"
#include "lightfem/vtk_writer.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

void print_usage() {
    std::cout << "LightFEM 0.6.0 - 2D Q4 linear elastic solver\n\n"
              << "Usage:\n"
              << "  lightfem check <model.fem>\n"
              << "  lightfem solve <model.fem> [result.vtk]\n"
              << "  lightfem --help\n"
              << "  lightfem --version\n\n"
              << "solve uses 2x2 Gauss Q4 stiffness, preallocated CSR assembly,\n"
              << "symmetric Dirichlet elimination and Jacobi-preconditioned CG.\n";
}

void print_model_summary(const lightfem::Model& model) {
    const lightfem::ModelSummary summary = model.summary();
    std::cout << "Model: " << model.source_name << '\n'
              << "  nodes=" << summary.nodes << ", elements=" << summary.elements
              << ", materials=" << summary.materials << ", boundary_records="
              << summary.boundary_records << ", nodal_loads=" << summary.loads
              << ", body_forces=" << summary.body_forces << ", edge_loads=" << summary.edge_loads << '\n'
              << "  plane_mode=" << lightfem::to_string(model.settings.plane_mode)
              << ", solver=" << model.settings.solver
              << ", tolerance=" << model.settings.tolerance
              << ", max_iterations=" << model.settings.max_iterations
              << ", reordering=" << model.settings.reordering
              << ", assembly=" << model.settings.assembly
              << ", threads=" << model.settings.threads
              << ", deformation_scale=" << model.settings.deformation_scale << '\n';
}

void print_diagnostics(const lightfem::ValidationReport& report) {
    for (const lightfem::Diagnostic& diagnostic : report.diagnostics()) {
        std::cout << lightfem::to_string(diagnostic.severity) << " [" << diagnostic.code << "]";
        if (diagnostic.source_line != 0) std::cout << " line " << diagnostic.source_line;
        std::cout << ": " << diagnostic.message << '\n';
    }
}

}  // namespace

int main(const int argc, char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::cout << "LightFEM 0.6.0\n";
        return 0;
    }
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        print_usage();
        return 0;
    }
    if (argc < 3 || argc > 4) {
        print_usage();
        return 1;
    }

    const std::string command = argv[1];
    if (command != "check" && command != "solve") {
        print_usage();
        return 1;
    }
    if (command == "check" && argc != 3) {
        print_usage();
        return 1;
    }

    try {
        lightfem::Model model = lightfem::InputParser{}.parse_file(argv[2]);
        const lightfem::ValidationReport report = lightfem::ModelValidator{}.validate(model);
        print_model_summary(model);
        print_diagnostics(report);
        if (!report.ok()) {
            std::cout << "INVALID: " << report.error_count() << " error(s), "
                      << report.warning_count() << " warning(s)\n";
            return 3;
        }

        model.sort_and_index();
        if (command == "check") {
            std::cout << "VALID: " << model.summary().dofs << " displacement DOFs, "
                      << report.warning_count() << " warning(s)\n";
            return 0;
        }

        const std::string output_path = argc == 4 ? argv[3] : model.settings.output_file;
        const lightfem::AnalysisResult result = lightfem::LinearElasticSolver{}.solve(model);
        lightfem::write_vtk(model, result, output_path);
        const std::string csv_prefix = model.settings.csv_prefix.empty()
                                           ? std::filesystem::path(output_path).replace_extension().string()
                                           : model.settings.csv_prefix;
        const lightfem::CsvExportPaths csv_paths = lightfem::write_csv_results(model, result, csv_prefix);
        std::cout << "SOLVED: dofs=" << result.displacement.size()
                  << ", nnz=" << result.stiffness.nonzeros()
                  << ", pcg_iterations=" << result.solver.iterations
                  << ", initial_residual=" << result.solver.initial_residual
                  << ", final_residual=" << result.solver.final_residual << '\n'
                  << "PROFILE: ordering=" << result.performance.reordering
                  << ", before(nnz=" << result.performance.before.nonzeros
                  << ", bandwidth=" << result.performance.before.bandwidth
                  << ", profile=" << result.performance.before.profile << ')'
                  << ", after(nnz=" << result.performance.after.nonzeros
                  << ", bandwidth=" << result.performance.after.bandwidth
                  << ", profile=" << result.performance.after.profile << ")\n"
                  << "PARALLEL: assembly=" << result.performance.assembly
                  << ", requested_threads=" << result.performance.requested_threads
                  << ", threads_used=" << result.performance.threads_used
                  << ", colors=" << result.performance.coloring_count
                  << ", element_stiffness_ms=" << result.performance.element_stiffness_ms
                  << ", assembly_ms=" << result.performance.assembly_ms << '\n'
                  << "VTK: " << output_path << " (deformation_scale=" << model.settings.deformation_scale << ")\n"
                  << "CSV: " << csv_paths.nodes << ", " << csv_paths.gauss_points << '\n';
        return 0;
    } catch (const lightfem::ParseError& error) {
        std::cerr << "PARSE ERROR: " << error.what() << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "FATAL ERROR: " << error.what() << '\n';
        return 4;
    }
}
