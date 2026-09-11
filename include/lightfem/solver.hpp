#pragma once

#include "lightfem/element.hpp"
#include "lightfem/sparse.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace lightfem {

struct LinearSolverReport {
    bool converged{false};
    std::size_t iterations{0};
    double initial_residual{0.0};
    double final_residual{0.0};
    std::string message;
};

struct ElementResult {
    int element_id{0};
    ElementField field{};
    std::array<GaussPointField, 4> gauss_points{};
};

struct NodeResult {
    int node_id{0};
    std::array<double, 3> stress{{0.0, 0.0, 0.0}};
    double von_mises{0.0};
};

struct PerformanceReport {
    std::string reordering{"none"};
    CsrProfile before;
    CsrProfile after;
    std::string assembly{"serial"};
    std::size_t requested_threads{0};
    std::size_t threads_used{1};
    std::size_t coloring_count{0};
    double element_stiffness_ms{0.0};
    double assembly_ms{0.0};
};

struct AnalysisResult {
    // Unmodified global system, retained for reaction recovery and inspection.
    CsrMatrix stiffness;
    std::vector<double> load;
    std::vector<double> displacement;
    std::vector<double> reaction;
    std::vector<ElementResult> element_results;
    std::vector<NodeResult> node_results;
    LinearSolverReport solver;
    PerformanceReport performance;
};

// Linear static solver for validated, indexed 2D Q4 models.
class LinearElasticSolver {
public:
    [[nodiscard]] AnalysisResult solve(Model& model) const;
};

}  // namespace lightfem
