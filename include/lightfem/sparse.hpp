#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace lightfem {

struct CsrAssemblyPattern;
// Element DOF lists may have different lengths (8 for Q4, 4 for a two-node bar),
// so each entry is a plain vector rather than a fixed-size array.
using ElementDofList = std::vector<std::vector<std::size_t>>;

[[nodiscard]] CsrAssemblyPattern build_csr_assembly_pattern(std::size_t dimension,
                                                            const ElementDofList& element_dofs);

// Symmetric finite-element matrices are stored as a complete CSR pattern so
// Dirichlet elimination can zero rows and columns without changing the pattern.
class CsrMatrix {
public:
    CsrMatrix() = default;

    [[nodiscard]] static CsrMatrix from_element_dofs(std::size_t dimension,
                                                     const ElementDofList& element_dofs);

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return values_.size(); }
    [[nodiscard]] const std::vector<std::size_t>& row_offsets() const noexcept { return row_offsets_; }
    [[nodiscard]] const std::vector<std::size_t>& column_indices() const noexcept {
        return column_indices_;
    }
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }
    [[nodiscard]] std::vector<double>& values() noexcept { return values_; }

    void add(std::size_t row, std::size_t column, double value);
    void set(std::size_t row, std::size_t column, double value);
    [[nodiscard]] double at(std::size_t row, std::size_t column) const;
    [[nodiscard]] bool contains(std::size_t row, std::size_t column) const;
    [[nodiscard]] std::size_t value_index(std::size_t row, std::size_t column) const;
    void zero_values() noexcept;
    [[nodiscard]] std::vector<double> multiply(const std::vector<double>& vector) const;
    // Same product into a caller-owned buffer. The iterative solver calls this
    // hundreds of times per analysis, so it must not allocate or zero the
    // result: rows are overwritten outright.
    void multiply_into(const std::vector<double>& vector, std::vector<double>& result) const;
    [[nodiscard]] std::vector<double> diagonal() const;

private:
    friend CsrAssemblyPattern build_csr_assembly_pattern(std::size_t dimension,
                                                         const ElementDofList& element_dofs);

    std::size_t rows_{0};
    std::vector<std::size_t> row_offsets_{0};
    std::vector<std::size_t> column_indices_;
    std::vector<double> values_;
};

// Two-dimensional elasticity couples the two DOFs of a node, so the sparsity
// pattern of K is a 2x2-block pattern: block (I, J) exists exactly when some
// entry of rows {2I, 2I+1} and columns {2J, 2J+1} exists. Storing blocks keeps
// one column index per block instead of one per scalar entry, which is the
// dominant term in the bytes moved per matrix-vector product.
class BsrMatrix {
public:
    static constexpr std::size_t kBlock = 2;

    // Requires an even row count: DOF 2I and 2I+1 belong to the same node.
    [[nodiscard]] static BsrMatrix from_csr(const CsrMatrix& matrix);

    [[nodiscard]] std::size_t block_rows() const noexcept { return block_row_offsets_.empty()
                                                                      ? 0
                                                                      : block_row_offsets_.size() - 1; }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t block_nonzeros() const noexcept { return block_columns_.size(); }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return block_columns_.size() * kBlock * kBlock; }
    [[nodiscard]] const std::vector<std::size_t>& block_row_offsets() const noexcept {
        return block_row_offsets_;
    }
    [[nodiscard]] const std::vector<std::uint32_t>& block_columns() const noexcept {
        return block_columns_;
    }
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }

    void multiply_into(const std::vector<double>& vector, std::vector<double>& result) const;
    [[nodiscard]] std::vector<double> multiply(const std::vector<double>& vector) const;
    // Row-major 2x2 blocks along the diagonal: 4 doubles per block row.
    [[nodiscard]] std::vector<double> diagonal_blocks() const;

    // Storage footprint of one matrix-vector product, used by the data-movement
    // accounting in the scaling study.
    [[nodiscard]] static constexpr std::size_t block_bytes() noexcept {
        return kBlock * kBlock * sizeof(double) + sizeof(std::uint32_t);
    }

private:
    std::size_t rows_{0};
    std::vector<std::size_t> block_row_offsets_{0};
    std::vector<std::uint32_t> block_columns_;
    std::vector<double> values_;
};

struct CsrProfile {
    std::size_t dimension{0};
    std::size_t nonzeros{0};
    std::size_t bandwidth{0};
    std::size_t profile{0};
};

// Each map converts a local stiffness entry (n*n, row-major) directly to an
// index in CsrMatrix::values(), avoiding any row search during numeric assembly.
struct CsrAssemblyPattern {
    CsrMatrix matrix;
    std::vector<std::vector<std::size_t>> element_value_indices;
};

[[nodiscard]] CsrProfile csr_profile(const CsrMatrix& matrix);

}  // namespace lightfem
