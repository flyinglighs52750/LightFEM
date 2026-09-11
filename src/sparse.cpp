#include "lightfem/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

#ifdef LIGHTFEM_HAS_OPENMP
#include <omp.h>
#endif

namespace lightfem {

CsrAssemblyPattern build_csr_assembly_pattern(
    const std::size_t dimension, const std::vector<std::array<std::size_t, 8>>& element_dofs) {
    std::vector<std::unordered_set<std::size_t>> adjacency(dimension);
    for (std::size_t row = 0; row < dimension; ++row) adjacency[row].insert(row);
    for (const auto& dofs : element_dofs) {
        for (const std::size_t row : dofs) {
            if (row >= dimension) throw std::out_of_range("element DOF is outside CSR matrix dimension");
            for (const std::size_t column : dofs) {
                if (column >= dimension) throw std::out_of_range("element DOF is outside CSR matrix dimension");
                adjacency[row].insert(column);
            }
        }
    }

    CsrMatrix matrix;
    matrix.rows_ = dimension;
    matrix.row_offsets_.assign(dimension + 1, 0);
    for (std::size_t row = 0; row < dimension; ++row) {
        matrix.row_offsets_[row + 1] = matrix.row_offsets_[row] + adjacency[row].size();
    }
    matrix.column_indices_.resize(matrix.row_offsets_.back());
    for (std::size_t row = 0; row < dimension; ++row) {
        std::vector<std::size_t> columns(adjacency[row].begin(), adjacency[row].end());
        std::sort(columns.begin(), columns.end());
        std::copy(columns.begin(), columns.end(), matrix.column_indices_.begin() +
                                                     static_cast<std::ptrdiff_t>(matrix.row_offsets_[row]));
    }
    matrix.values_.assign(matrix.column_indices_.size(), 0.0);
    CsrAssemblyPattern pattern;
    pattern.matrix = std::move(matrix);
    pattern.element_value_indices.resize(element_dofs.size());
    for (std::size_t element = 0; element < element_dofs.size(); ++element) {
        for (std::size_t row = 0; row < 8; ++row) {
            for (std::size_t column = 0; column < 8; ++column) {
                pattern.element_value_indices[element][8 * row + column] =
                    pattern.matrix.value_index(element_dofs[element][row], element_dofs[element][column]);
            }
        }
    }
    return pattern;
}

CsrMatrix CsrMatrix::from_element_dofs(
    const std::size_t dimension, const std::vector<std::array<std::size_t, 8>>& element_dofs) {
    return build_csr_assembly_pattern(dimension, element_dofs).matrix;
}

std::size_t CsrMatrix::value_index(const std::size_t row, const std::size_t column) const {
    if (row >= rows_ || column >= rows_) throw std::out_of_range("CSR matrix index is out of range");
    const auto begin = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
    const auto end = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1]);
    const auto it = std::lower_bound(begin, end, column);
    if (it == end || *it != column) {
        throw std::logic_error("requested entry is absent from preallocated CSR pattern");
    }
    return static_cast<std::size_t>(std::distance(column_indices_.begin(), it));
}

void CsrMatrix::add(const std::size_t row, const std::size_t column, const double value) {
    values_[value_index(row, column)] += value;
}

void CsrMatrix::set(const std::size_t row, const std::size_t column, const double value) {
    values_[value_index(row, column)] = value;
}

double CsrMatrix::at(const std::size_t row, const std::size_t column) const {
    return values_[value_index(row, column)];
}

bool CsrMatrix::contains(const std::size_t row, const std::size_t column) const {
    if (row >= rows_ || column >= rows_) return false;
    const auto begin = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
    const auto end = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1]);
    const auto it = std::lower_bound(begin, end, column);
    return it != end && *it == column;
}

void CsrMatrix::zero_values() noexcept { std::fill(values_.begin(), values_.end(), 0.0); }

std::vector<double> CsrMatrix::multiply(const std::vector<double>& vector) const {
    if (vector.size() != rows_) throw std::invalid_argument("CSR multiply vector size mismatch");
    std::vector<double> result(rows_, 0.0);
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long row_index = 0; row_index < static_cast<long long>(rows_); ++row_index) {
        const std::size_t row = static_cast<std::size_t>(row_index);
        for (std::size_t index = row_offsets_[row]; index < row_offsets_[row + 1]; ++index) {
            result[row] += values_[index] * vector[column_indices_[index]];
        }
    }
    return result;
}

std::vector<double> CsrMatrix::diagonal() const {
    std::vector<double> result(rows_, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) result[row] = at(row, row);
    return result;
}

CsrProfile csr_profile(const CsrMatrix& matrix) {
    CsrProfile result;
    result.dimension = matrix.rows();
    result.nonzeros = matrix.nonzeros();
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        std::size_t minimum_column = row;
        for (std::size_t index = matrix.row_offsets()[row]; index < matrix.row_offsets()[row + 1]; ++index) {
            const std::size_t column = matrix.column_indices()[index];
            const std::size_t distance = row > column ? row - column : column - row;
            result.bandwidth = std::max(result.bandwidth, distance);
            if (column <= row) minimum_column = std::min(minimum_column, column);
        }
        result.profile += row - minimum_column;
    }
    return result;
}

}  // namespace lightfem
