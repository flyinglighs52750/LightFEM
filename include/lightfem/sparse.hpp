#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace lightfem {

struct CsrAssemblyPattern;
[[nodiscard]] CsrAssemblyPattern build_csr_assembly_pattern(
    std::size_t dimension, const std::vector<std::array<std::size_t, 8>>& element_dofs);

// Symmetric finite-element matrices are stored as a complete CSR pattern so
// Dirichlet elimination can zero rows and columns without changing the pattern.
class CsrMatrix {
public:
    CsrMatrix() = default;

    [[nodiscard]] static CsrMatrix from_element_dofs(
        std::size_t dimension, const std::vector<std::array<std::size_t, 8>>& element_dofs);

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
    [[nodiscard]] std::vector<double> diagonal() const;

private:
    friend CsrAssemblyPattern build_csr_assembly_pattern(
        std::size_t dimension, const std::vector<std::array<std::size_t, 8>>& element_dofs);

    std::size_t rows_{0};
    std::vector<std::size_t> row_offsets_{0};
    std::vector<std::size_t> column_indices_;
    std::vector<double> values_;
};

struct CsrProfile {
    std::size_t dimension{0};
    std::size_t nonzeros{0};
    std::size_t bandwidth{0};
    std::size_t profile{0};
};

// Each 8x8 map converts a local Q4 stiffness entry directly to an index in
// CsrMatrix::values(), avoiding any row search during numeric assembly.
struct CsrAssemblyPattern {
    CsrMatrix matrix;
    std::vector<std::array<std::size_t, 64>> element_value_indices;
};

[[nodiscard]] CsrProfile csr_profile(const CsrMatrix& matrix);

}  // namespace lightfem
