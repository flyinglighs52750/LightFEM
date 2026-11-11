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
    // Iterative-phase breakdown. End-to-end studies need the split between the
    // solve loop and everything else, because assembly is a small share of it.
    std::string preconditioner{"jacobi"};
    double solve_ms{0.0};
    double matvec_ms{0.0};
    std::size_t matvec_calls{0};
    std::string storage{"csr"};
    // Bytes the matrix alone moves per product, plus the setup cost paid once.
    std::size_t matrix_bytes_per_matvec{0};
    double storage_convert_ms{0.0};
    double preconditioner_setup_ms{0.0};
    double preconditioner_apply_ms{0.0};
    std::size_t subdomains{0};
    // Share of nonzeros whose row and column share a subdomain. Values near one
    // mean the block preconditioner sees most of the coupling.
    double capture_ratio{0.0};
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
    // Measuring the bandwidth/profile of the system walks every nonzero. When a
    // reordering is requested that has to be done twice (before and after), so
    // the cost is charged to its own line instead of hiding inside the solve.
    double profile_ms{0.0};
    // Building the CSR pattern and the precomputed element-to-value map is setup
    // work that the numeric assembly cannot start without, so it is reported
    // separately from the assembly loop itself.
    double pattern_ms{0.0};
    double element_stiffness_ms{0.0};
    double assembly_ms{0.0};
    // Post-assembly phases, so the end-to-end budget can be closed exactly.
    double boundary_ms{0.0};
    double postprocess_ms{0.0};
    // The result object keeps its own copy of the unconstrained matrix so the
    // caller can inspect it after the solve. That copy is pure memory traffic
    // and is charged to no other stage.
    double retain_matrix_ms{0.0};
    double total_ms{0.0};
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
