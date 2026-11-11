#pragma once

#include "lightfem/sparse.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace lightfem {

// A block preconditioner needs the unknowns split into subdomains whose
// diagonal blocks are inverted densely. Two splits are provided because they
// answer different questions:
//
//   Chunk - consecutive node indices. Whether a chunk is spatially compact
//           depends entirely on the reordering, so this split measures what a
//           bandwidth-reducing permutation does for the preconditioner.
//   Rcb   - recursive coordinate bisection on the node coordinates. Compact
//           regardless of how the unknowns are numbered, so it is the
//           reordering-independent reference.
//
// Splitting is defined on nodes, not on scalar unknowns, so every subdomain is
// a union of complete (u, v) pairs.
struct SubdomainSplit {
    std::vector<std::uint32_t> node_subdomain;
    std::size_t count{0};
};

[[nodiscard]] SubdomainSplit split_by_chunks(std::size_t node_count, std::size_t nodes_per_subdomain);

[[nodiscard]] SubdomainSplit split_by_coordinates(const std::vector<double>& x,
                                                  const std::vector<double>& y,
                                                  std::size_t nodes_per_subdomain);

// Interior fraction: the share of stored nonzeros whose row and column fall in
// the same subdomain. It is computable from the sparsity pattern alone, before
// any iteration, and it is the quantity a block preconditioner actually
// exploits. Bandwidth and profile measure a direct method's fill, which a
// preconditioned Krylov solver never produces.
[[nodiscard]] double subdomain_capture_ratio(const CsrMatrix& matrix,
                                             const std::vector<std::uint32_t>& dof_subdomain);

struct BlockJacobiStats {
    std::size_t subdomains{0};
    std::size_t dense_dofs{0};       // sum of subdomain sizes
    std::size_t factor_nonzeros{0};
    double capture_ratio{0.0};
    double setup_ms{0.0};
    double apply_ms{0.0};
    std::size_t apply_calls{0};
    std::size_t degenerate_subdomains{0};
};

// Additive Schwarz without overlap: one dense Cholesky factorization per
// subdomain, applied as a block-diagonal inverse.
class SubdomainBlockJacobi {
public:
    SubdomainBlockJacobi() = default;
    SubdomainBlockJacobi(const CsrMatrix& matrix, std::vector<std::uint32_t> dof_subdomain,
                         std::size_t subdomain_count);

    [[nodiscard]] bool empty() const noexcept { return subdomain_count_ == 0; }
    void apply(const std::vector<double>& input, std::vector<double>& output) const;
    [[nodiscard]] BlockJacobiStats stats() const { return stats_; }
    void account_apply(double milliseconds) const;

private:
    std::size_t dof_count_{0};
    std::size_t subdomain_count_{0};
    std::vector<std::size_t> dof_begin_{0};
    std::vector<std::uint32_t> dofs_;
    std::vector<std::size_t> factor_begin_{0};
    std::vector<double> factor_;
    std::size_t largest_block_{0};
    mutable BlockJacobiStats stats_{};
};

struct IncompleteCholeskyStats {
    double setup_ms{0.0};
    double apply_ms{0.0};
    std::size_t apply_calls{0};
    std::size_t factor_entries{0};
    std::size_t breakdowns{0};
};

// Incomplete Cholesky with zero fill: L keeps exactly the pattern of the lower
// triangle of K, so both triangular sweeps cost the same order as one
// matrix-vector product instead of growing with a subdomain size. This is the
// preconditioner that spends the memory budget on iterations rather than on
// dense blocks.
class IncompleteCholesky0 {
public:
    explicit IncompleteCholesky0(const CsrMatrix& matrix);

    void apply(const std::vector<double>& input, std::vector<double>& output) const;
    [[nodiscard]] IncompleteCholeskyStats stats() const { return stats_; }
    void account_apply(double milliseconds) const;

private:
    std::size_t rows_{0};
    std::vector<std::size_t> row_offsets_{0};
    std::vector<std::uint32_t> columns_;
    std::vector<double> values_;
    std::vector<double> inverse_diagonal_;
    mutable IncompleteCholeskyStats stats_{};
};

}  // namespace lightfem
