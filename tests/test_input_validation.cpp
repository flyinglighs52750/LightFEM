#include "lightfem/input_parser.hpp"
#include "lightfem/validation.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

int failures = 0;

void expect(const bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

bool has_code(const lightfem::ValidationReport& report, const std::string& code) {
    for (const auto& diagnostic : report.diagnostics()) {
        if (diagnostic.code == code) return true;
    }
    return false;
}

lightfem::Model parse_text(const std::string& text) {
    std::istringstream input(text);
    return lightfem::InputParser{}.parse(input, "test.fem");
}

const char* valid_model = R"FEM(
[settings]
plane_mode plane_stress
solver pcg
tolerance 1e-10
max_iterations 200
output test.vtk
[materials]
1 210e9 0.3 0.01
[nodes]
4 0 1
2 1 0
1 0 0
3 1 1
[elements]
1 1 2 3 4 1
[boundary]
1 0 0
4 0 0
[loads]
2 0 -1
)FEM";

void test_valid_model_and_indexing() {
    auto model = parse_text(valid_model);
    const auto report = lightfem::ModelValidator{}.validate(model);
    expect(report.ok(), "valid model must pass semantic validation");
    model.sort_and_index();
    expect(model.is_indexed(), "model must report indexed state");
    expect(model.find_node(1) != nullptr, "node lookup must find ID 1");
    expect(model.find_node(1)->dof_x == 0 && model.find_node(1)->dof_y == 1,
           "DOFs must be assigned by ascending node ID");
    expect(model.find_node(4)->dof_x == 6 && model.find_node(4)->dof_y == 7,
           "last node DOFs must be deterministic");
    expect(model.summary().dofs == 8, "four nodes must create eight displacement DOFs");
}

void test_clockwise_element_is_rejected() {
    const std::string text = R"FEM(
[materials]
1 1e6 0.25 1
[nodes]
1 0 0
2 1 0
3 1 1
4 0 1
[elements]
1 1 4 3 2 1
[boundary]
1 0 0
4 0 0
[loads]
2 0 -1
)FEM";
    const auto report = lightfem::ModelValidator{}.validate(parse_text(text));
    expect(!report.ok(), "clockwise element must be invalid");
    expect(has_code(report, "E_NONPOSITIVE_JACOBIAN"),
           "clockwise element must report E_NONPOSITIVE_JACOBIAN");
}

void test_unknown_node_is_rejected() {
    const std::string text = R"FEM(
[materials]
1 1e6 0.25 1
[nodes]
1 0 0
2 1 0
3 1 1
4 0 1
[elements]
1 1 2 99 4 1
[boundary]
1 0 0
4 0 0
[loads]
2 0 -1
)FEM";
    const auto report = lightfem::ModelValidator{}.validate(parse_text(text));
    expect(has_code(report, "E_UNKNOWN_ELEMENT_NODE"),
           "unknown element node must be diagnosed");
}

void test_duplicate_node_id_is_rejected() {
    const std::string text = R"FEM(
[materials]
1 1e6 0.25 1
[nodes]
1 0 0
1 1 0
3 1 1
4 0 1
[elements]
1 1 3 4 1 1
[boundary]
1 0 0
[loads]
3 0 -1
)FEM";
    const auto report = lightfem::ModelValidator{}.validate(parse_text(text));
    expect(has_code(report, "E_DUPLICATE_NODE_ID"), "duplicate node ID must be diagnosed");
}

void test_parse_error_has_line_number() {
    const std::string text = "[nodes]\n1 not_a_number 0\n";
    try {
        (void)parse_text(text);
        expect(false, "invalid numeric token must throw ParseError");
    } catch (const lightfem::ParseError& error) {
        expect(error.line() == 2, "ParseError must preserve source line number");
    }
}

void test_constraints_remove_rigid_modes() {
    std::string text = valid_model;
    const std::string old_constraints = "1 0 0\n4 0 0";
    const std::string weak_constraints = "1 0 free\n4 0 free";
    text.replace(text.find(old_constraints), old_constraints.size(), weak_constraints);
    const auto report = lightfem::ModelValidator{}.validate(parse_text(text));
    expect(has_code(report, "E_RIGID_BODY_MODES"),
           "insufficient constraints must report remaining rigid-body modes");
}

}  // namespace

int main() {
    test_valid_model_and_indexing();
    test_clockwise_element_is_rejected();
    test_unknown_node_is_rejected();
    test_duplicate_node_id_is_rejected();
    test_parse_error_has_line_number();
    test_constraints_remove_rigid_modes();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }
    std::cout << "All LightFEM input and validation tests passed.\n";
    return 0;
}

