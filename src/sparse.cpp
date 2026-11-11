#include "lightfem/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

#ifdef LIGHTFEM_HAS_OPENMP
#include <omp.h>
#endif

namespace lightfem {

CsrAssemblyPattern build_csr_assembly_pattern(const std::size_t dimension,
                                              const ElementDofList& element_dofs) {
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
        const std::size_t count = element_dofs[element].size();
        if (count == 0) throw std::invalid_argument("element DOF list must not be empty");
        pattern.element_value_indices[element].resize(count * count);
        for (std::size_t row = 0; row < count; ++row) {
            for (std::size_t column = 0; column < count; ++column) {
                pattern.element_value_indices[element][count * row + column] =
                    pattern.matrix.value_index(element_dofs[element][row], element_dofs[element][column]);
            }
        }
    }
    return pattern;
}

CsrMatrix CsrMatrix::from_element_dofs(const std::size_t dimension,
                                      const ElementDofList& element_dofs) {
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
    std::vector<double> result(rows_, 0.0);
    multiply_into(vector, result);
    return result;
}

void CsrMatrix::multiply_into(const std::vector<double>& vector, std::vector<double>& result) const {
    if (vector.size() != rows_) throw std::invalid_argument("CSR multiply vector size mismatch");
    if (result.size() != rows_) throw std::invalid_argument("CSR multiply result size mismatch");
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long row_index = 0; row_index < static_cast<long long>(rows_); ++row_index) {
        const std::size_t row = static_cast<std::size_t>(row_index);
        // Accumulating in a local keeps the destination out of the inner loop,
        // which otherwise writes to it once per nonzero.
        double sum = 0.0;
        for (std::size_t index = row_offsets_[row]; index < row_offsets_[row + 1]; ++index) {
            sum += values_[index] * vector[column_indices_[index]];
        }
        result[row] = sum;
    }
}

std::vector<double> CsrMatrix::diagonal() const {
    std::vector<double> result(rows_, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) result[row] = at(row, row);
    return result;
}

BsrMatrix BsrMatrix::from_csr(const CsrMatrix& matrix) {
    if (matrix.rows() % kBlock != 0) {
        throw std::invalid_argument("BSR conversion needs an even number of rows (u,v DOF pairs)");
    }
    constexpr std::size_t kMaxColumn = std::numeric_limits<std::size_t>::max();
    BsrMatrix result;
    result.rows_ = matrix.rows();

    const std::size_t block_row_count = matrix.rows() / kBlock;
    result.block_row_offsets_.assign(1, 0);
    result.block_row_offsets_.reserve(block_row_count + 1);
    // Roughly nine blocks per row for a quadrilateral mesh, four values each.
    result.values_.reserve(matrix.nonzeros() * 2);
    result.block_columns_.reserve(matrix.nonzeros() / 2);

    const auto& offsets = matrix.row_offsets();
    const auto& columns = matrix.column_indices();
    const auto& values = matrix.values();

    for (std::size_t block_row = 0; block_row < block_row_count; ++block_row) {
        const std::size_t row0 = kBlock * block_row;
        const std::size_t row1 = row0 + 1;
        std::size_t p0 = offsets[row0];
        std::size_t p1 = offsets[row1];
        const std::size_t e0 = offsets[row0 + 1];
        const std::size_t e1 = offsets[row1 + 1];

        // Both rows are sorted by column, so a merge walk yields the union of
        // their block columns in one pass without sorting or binary search.
        while (p0 < e0 || p1 < e1) {
            const std::size_t c0 = p0 < e0 ? columns[p0] : kMaxColumn;
            const std::size_t c1 = p1 < e1 ? columns[p1] : kMaxColumn;
            const std::size_t block_column = std::min(c0, c1) / kBlock;
            result.block_columns_.push_back(static_cast<std::uint32_t>(block_column));

            double block[kBlock * kBlock] = {0.0, 0.0, 0.0, 0.0};
            while (p0 < e0 && columns[p0] / kBlock == block_column) {
                block[columns[p0] - kBlock * block_column] = values[p0];
                ++p0;
            }
            while (p1 < e1 && columns[p1] / kBlock == block_column) {
                block[kBlock + columns[p1] - kBlock * block_column] = values[p1];
                ++p1;
            }
            result.values_.insert(result.values_.end(), block, block + kBlock * kBlock);
        }
        result.block_row_offsets_.push_back(result.block_columns_.size());
    }
    return result;
}

void BsrMatrix::multiply_into(const std::vector<double>& vector, std::vector<double>& result) const {
    if (vector.size() != rows_) throw std::invalid_argument("BSR multiply vector size mismatch");
    if (result.size() != rows_) result.assign(rows_, 0.0);
    const std::size_t block_row_count = block_rows();
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long block_row_index = 0; block_row_index < static_cast<long long>(block_row_count);
         ++block_row_index) {
        const std::size_t block_row = static_cast<std::size_t>(block_row_index);
        double sum0 = 0.0;
        double sum1 = 0.0;
        const std::size_t begin = block_row_offsets_[block_row];
        const std::size_t end = block_row_offsets_[block_row + 1];
        for (std::size_t index = begin; index < end; ++index) {
            const std::size_t column = kBlock * static_cast<std::size_t>(block_columns_[index]);
            const double* block = &values_[kBlock * kBlock * index];
            const double x0 = vector[column];
            const double x1 = vector[column + 1];
            sum0 += block[0] * x0 + block[1] * x1;
            sum1 += block[2] * x0 + block[3] * x1;
        }
        result[kBlock * block_row] = sum0;
        result[kBlock * block_row + 1] = sum1;
    }
}

std::vector<double> BsrMatrix::multiply(const std::vector<double>& vector) const {
    std::vector<double> result(rows_, 0.0);
    multiply_into(vector, result);
    return result;
}

std::vector<double> BsrMatrix::diagonal_blocks() const {
    std::vector<double> result(kBlock * kBlock * block_rows(), 0.0);
    for (std::size_t block_row = 0; block_row < block_rows(); ++block_row) {
        for (std::size_t index = block_row_offsets_[block_row]; index < block_row_offsets_[block_row + 1];
             ++index) {
            if (static_cast<std::size_t>(block_columns_[index]) != block_row) continue;
            for (std::size_t k = 0; k < kBlock * kBlock; ++k) {
                result[kBlock * kBlock * block_row + k] = values_[kBlock * kBlock * index + k];
            }
            break;
        }
    }
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
