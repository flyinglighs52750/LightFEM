#include "lightfem/element.hpp"
#include "lightfem/input_parser.hpp"
#include "lightfem/solver.hpp"
#include "lightfem/validation.hpp"
#include "lightfem/vtk_writer.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

int failures = 0;

void expect(const bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void expect_near(const double actual, const double expected, const double tolerance,
                 const std::string& message) {
    expect(std::abs(actual - expected) <= tolerance,
           message + " (actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected) + ")");
}

lightfem::Model parse_and_index(const std::string& text) {
    std::istringstream input(text);
    lightfem::Model model = lightfem::InputParser{}.parse(input, "solver_test.fem");
    const lightfem::ValidationReport report = lightfem::ModelValidator{}.validate(model);
    if (!report.ok()) throw std::runtime_error("solver test fixture did not pass validation");
    model.sort_and_index();
    return model;
}

const char* cantilever = R"FEM(
[settings]
plane_mode plane_stress
solver pcg
tolerance 1e-12
max_iterations 100
output ignored.vtk
[materials]
1 1000000 0.25 1
[nodes]
1 0 0
2 1 0
3 1 1
4 0 1
[elements]
1 1 2 3 4 1
[boundary]
1 0 0
4 0 0
[loads]
2 0 -1
)FEM";

void test_q4_stiffness_and_constant_strain() {
    const lightfem::Material material{1, 1000.0, 0.25, 1.0, 0};
    const std::array<lightfem::Vec2, 4> coordinates{{{0.0, 0.0}, {2.0, 0.0},
                                                       {2.0, 1.0}, {0.0, 1.0}}};
    const lightfem::ElementMatrix stiffness =
        lightfem::q4_stiffness(coordinates, material, lightfem::PlaneMode::PlaneStress);
    for (std::size_t row = 0; row < 8; ++row) {
        for (std::size_t column = 0; column < 8; ++column) {
            expect_near(stiffness[8 * row + column], stiffness[8 * column + row], 1.0e-10,
                        "Q4 stiffness must be symmetric");
        }
    }

    const std::array<double, 8> displacement{{0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0}};
    const lightfem::ElementField field = lightfem::q4_center_field(
        coordinates, material, lightfem::PlaneMode::PlaneStress, displacement);
    expect_near(field.strain[0], 1.0, 1.0e-12, "linear ux=x field must produce exx=1");
    expect_near(field.strain[1], 0.0, 1.0e-12, "linear ux=x field must produce eyy=0");
    expect_near(field.strain[2], 0.0, 1.0e-12, "linear ux=x field must produce gamma_xy=0");
    expect_near(field.stress[0], 1000.0 / (1.0 - 0.25 * 0.25), 1.0e-9,
                "plane-stress constitutive response must be correct");
}

void test_consistent_body_and_edge_loads() {
    const lightfem::Material material{1, 1000.0, 0.25, 2.0, 0};
    const std::array<lightfem::Vec2, 4> coordinates{{{0.0, 0.0}, {1.0, 0.0},
                                                       {1.0, 1.0}, {0.0, 1.0}}};
    const auto body = lightfem::q4_body_force(coordinates, material, 3.0, -4.0);
    for (std::size_t node = 0; node < 4; ++node) {
        expect_near(body[2 * node], 1.5, 1.0e-12, "unit-square body force must distribute fx consistently");
        expect_near(body[2 * node + 1], -2.0, 1.0e-12, "unit-square body force must distribute fy consistently");
    }
    const auto traction = lightfem::q4_edge_traction(coordinates, material, 2, 5.0, -7.0);
    expect_near(traction[2], 5.0, 1.0e-12, "edge traction must load first edge endpoint");
    expect_near(traction[3], -7.0, 1.0e-12, "edge traction must load first endpoint y");
    expect_near(traction[4], 5.0, 1.0e-12, "edge traction must load second edge endpoint");
    expect_near(traction[5], -7.0, 1.0e-12, "edge traction must load second endpoint y");
    const auto pressure = lightfem::q4_edge_pressure(coordinates, material, 2, 6.0);
    expect_near(pressure[2], -6.0, 1.0e-12, "positive right-edge pressure must act inward at node 2");
    expect_near(pressure[4], -6.0, 1.0e-12, "positive right-edge pressure must act inward at node 3");
}

void test_csr_pcg_boundary_and_vtk() {
    lightfem::Model model = parse_and_index(cantilever);
    const lightfem::AnalysisResult result = lightfem::LinearElasticSolver{}.solve(model);
    expect(result.solver.converged, "Jacobi-PCG must converge for constrained Q4 cantilever");
    expect(result.stiffness.rows() == 8, "single Q4 model must create 8 global DOFs");
    expect(result.stiffness.nonzeros() == 64, "single Q4 must preallocate a full 8x8 CSR block");
    expect(result.displacement.size() == 8 && result.reaction.size() == 8,
           "solver output vectors must have one value per global DOF");

    const lightfem::Node* fixed_bottom = model.find_node(1);
    const lightfem::Node* fixed_top = model.find_node(4);
    const lightfem::Node* loaded = model.find_node(2);
    expect_near(result.displacement[fixed_bottom->dof_x], 0.0, 1.0e-14,
                "fixed ux at node 1 must be enforced");
    expect_near(result.displacement[fixed_bottom->dof_y], 0.0, 1.0e-14,
                "fixed uy at node 1 must be enforced");
    expect_near(result.displacement[fixed_top->dof_x], 0.0, 1.0e-14,
                "fixed ux at node 4 must be enforced");
    expect_near(result.displacement[fixed_top->dof_y], 0.0, 1.0e-14,
                "fixed uy at node 4 must be enforced");
    expect(result.displacement[loaded->dof_y] < 0.0, "downward nodal load must produce downward displacement");

    double reaction_y = 0.0;
    for (std::size_t i = 1; i < result.reaction.size(); i += 2) reaction_y += result.reaction[i];
    expect_near(reaction_y, 1.0, 1.0e-10, "reactions must balance the total vertical load");

    const std::filesystem::path output = std::filesystem::current_path() / "lightfem_solver_test.vtk";
    lightfem::write_vtk(model, result, output.string());
    std::ifstream vtk(output);
    std::string header;
    std::getline(vtk, header);
    expect(header == "# vtk DataFile Version 3.0", "VTK writer must produce legacy VTK header");
    expect(std::filesystem::file_size(output) > 0, "VTK writer must produce a non-empty result file");
    vtk.close();
    std::filesystem::remove(output);
}

const char* scrambled_mesh = R"FEM(
[settings]
plane_mode plane_stress
solver pcg
tolerance 1e-12
max_iterations 200
reordering none
output ignored.vtk
[materials]
1 1000000 0.25 1
[nodes]
1 0 0
9 1 0
3 2 0
8 0 1
2 1 1
7 2 1
4 0 2
6 1 2
5 2 2
[elements]
1 1 9 2 8 1
2 9 3 7 2 1
3 8 2 6 4 1
4 2 7 5 6 1
[boundary]
1 0 0
4 0 0
[loads]
3 0 -1
5 0 -1
)FEM";

void test_rcm_preserves_solution_and_improves_scrambled_pattern() {
    lightfem::Model baseline = parse_and_index(scrambled_mesh);
    lightfem::Model reordered = parse_and_index(scrambled_mesh);
    reordered.settings.reordering = "rcm";

    const lightfem::AnalysisResult baseline_result = lightfem::LinearElasticSolver{}.solve(baseline);
    const lightfem::AnalysisResult reordered_result = lightfem::LinearElasticSolver{}.solve(reordered);
    expect(reordered_result.performance.before.nonzeros == reordered_result.performance.after.nonzeros,
           "RCM must not change CSR nonzero count");
    expect(reordered_result.performance.after.bandwidth <= reordered_result.performance.before.bandwidth,
           "RCM must not increase bandwidth for the scrambled test mesh");
    expect(reordered_result.performance.after.profile <= reordered_result.performance.before.profile,
           "RCM must not increase profile for the scrambled test mesh");
    expect(reordered_result.performance.reordering == "rcm", "performance report must retain RCM mode");

    for (int node_id = 1; node_id <= 9; ++node_id) {
        const lightfem::Node* original_node = baseline.find_node(node_id);
        const lightfem::Node* reordered_node = reordered.find_node(node_id);
        expect_near(baseline_result.displacement[original_node->dof_x],
                    reordered_result.displacement[reordered_node->dof_x], 1.0e-13,
                    "RCM must preserve ux at each physical node");
        expect_near(baseline_result.displacement[original_node->dof_y],
                    reordered_result.displacement[reordered_node->dof_y], 1.0e-13,
                    "RCM must preserve uy at each physical node");
    }
}

void test_parallel_assembly_strategies_match_serial() {
    lightfem::Model baseline = parse_and_index(scrambled_mesh);
    baseline.settings.threads = 2;
    baseline.settings.assembly = "serial";
    const lightfem::AnalysisResult expected = lightfem::LinearElasticSolver{}.solve(baseline);

#ifdef LIGHTFEM_HAS_OPENMP
    for (const std::string strategy : {"atomic", "coloring", "thread_local"}) {
        lightfem::Model candidate = parse_and_index(scrambled_mesh);
        candidate.settings.threads = 2;
        candidate.settings.assembly = strategy;
        const lightfem::AnalysisResult actual = lightfem::LinearElasticSolver{}.solve(candidate);
        expect(actual.solver.converged, strategy + " assembly must converge");
        expect(actual.performance.assembly == strategy, strategy + " must be recorded in performance report");
        expect(actual.performance.threads_used >= 1, strategy + " must report OpenMP thread capacity");
        if (strategy == "coloring") {
            expect(actual.performance.coloring_count > 0, "coloring must report at least one color");
        }
        for (int node_id = 1; node_id <= 9; ++node_id) {
            const lightfem::Node* expected_node = baseline.find_node(node_id);
            const lightfem::Node* actual_node = candidate.find_node(node_id);
            expect_near(actual.displacement[actual_node->dof_x], expected.displacement[expected_node->dof_x],
                        1.0e-11, strategy + " must preserve ux");
            expect_near(actual.displacement[actual_node->dof_y], expected.displacement[expected_node->dof_y],
                        1.0e-11, strategy + " must preserve uy");
        }
    }
#endif
}

const char* extended_load_mesh = R"FEM(
[settings]
plane_mode plane_stress
solver pcg
tolerance 1e-12
max_iterations 100
output ignored.vtk
[materials]
1 1000000 0.25 1
[nodes]
1 0 0
2 1 0
3 1 1
4 0 1
[elements]
1 1 2 3 4 1
[boundary]
1 fixed
4 symmetry_x
4 roller_y
[body_forces]
1 0 -2
[edge_loads]
1 2 traction 0 -3
)FEM";

void test_extended_load_input_and_equilibrium() {
    lightfem::Model model = parse_and_index(extended_load_mesh);
    expect(model.body_forces.size() == 1 && model.edge_loads.size() == 1,
           "extended load sections must be parsed");
    const lightfem::AnalysisResult result = lightfem::LinearElasticSolver{}.solve(model);
    double reaction_y = 0.0;
    for (std::size_t i = 1; i < result.reaction.size(); i += 2) reaction_y += result.reaction[i];
    expect_near(reaction_y, 5.0, 1.0e-10,
                "body and uniform edge traction must be assembled into global equilibrium");
}

}  // namespace

int main() {
    test_q4_stiffness_and_constant_strain();
    test_consistent_body_and_edge_loads();
    test_csr_pcg_boundary_and_vtk();
    test_rcm_preserves_solution_and_improves_scrambled_pattern();
    test_parallel_assembly_strategies_match_serial();
    test_extended_load_input_and_equilibrium();
    if (failures != 0) {
        std::cerr << failures << " solver test assertion(s) failed\n";
        return 1;
    }
    std::cout << "All LightFEM solver tests passed.\n";
    return 0;
}
