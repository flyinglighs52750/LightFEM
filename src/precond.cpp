#include "lightfem/precond.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lightfem {

namespace {

double elapsed_ms(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

}  // namespace

SubdomainSplit split_by_chunks(const std::size_t node_count, const std::size_t nodes_per_subdomain) {
    SubdomainSplit split;
    const std::size_t size = std::max<std::size_t>(1, nodes_per_subdomain);
    split.node_subdomain.resize(node_count);
    for (std::size_t node = 0; node < node_count; ++node) {
        split.node_subdomain[node] = static_cast<std::uint32_t>(node / size);
    }
    split.count = node_count == 0 ? 0 : (node_count + size - 1) / size;
    return split;
}

SubdomainSplit split_by_coordinates(const std::vector<double>& x, const std::vector<double>& y,
                                    const std::size_t nodes_per_subdomain) {
    SubdomainSplit split;
    const std::size_t node_count = std::min(x.size(), y.size());
    split.node_subdomain.assign(node_count, 0);
    if (node_count == 0) return split;

    const std::size_t target = std::max<std::size_t>(1, nodes_per_subdomain);
    const std::size_t wanted = (node_count + target - 1) / target;
    // Recursive bisection is easier to write with an explicit stack than with
    // recursion: each entry is a contiguous range of the permutation that still
    // has to be cut, and the ranges never overlap.
    std::vector<std::size_t> order(node_count);
    for (std::size_t i = 0; i < node_count; ++i) order[i] = i;

    struct Range {
        std::size_t begin;
        std::size_t end;
    };
    // Halving always yields a power-of-two count, which would overshoot the
    // requested subdomain size by up to a factor of two and make the size
    // parameter hard to sweep. Taking the largest range off a heap each time
    // produces exactly the requested count with balanced sizes.
    const auto by_size = [](const Range& left, const Range& right) {
        return (left.end - left.begin) < (right.end - right.begin);
    };
    std::vector<Range> heap;
    heap.push_back(Range{0, node_count});
    std::size_t subdomains = 1;
    while (subdomains < wanted) {
        std::pop_heap(heap.begin(), heap.end(), by_size);
        const Range range = heap.back();
        heap.pop_back();
        const std::size_t size = range.end - range.begin;
        if (size < 2) {
            heap.push_back(range);
            std::push_heap(heap.begin(), heap.end(), by_size);
            break;  // every remaining range is a single node
        }
        // Cut across the longer side of the bounding box so both halves stay as
        // isotropic as the geometry allows.
        double min_x = x[order[range.begin]];
        double max_x = min_x;
        double min_y = y[order[range.begin]];
        double max_y = min_y;
        for (std::size_t i = range.begin + 1; i < range.end; ++i) {
            const std::size_t node = order[i];
            min_x = std::min(min_x, x[node]);
            max_x = std::max(max_x, x[node]);
            min_y = std::min(min_y, y[node]);
            max_y = std::max(max_y, y[node]);
        }
        const bool cut_on_x = (max_x - min_x) >= (max_y - min_y);
        const std::size_t middle = range.begin + size / 2;
        std::nth_element(order.begin() + static_cast<std::ptrdiff_t>(range.begin),
                         order.begin() + static_cast<std::ptrdiff_t>(middle),
                         order.begin() + static_cast<std::ptrdiff_t>(range.end),
                         [&](const std::size_t a, const std::size_t b) {
                             return cut_on_x ? x[a] < x[b] : y[a] < y[b];
                         });
        heap.push_back(Range{range.begin, middle});
        std::push_heap(heap.begin(), heap.end(), by_size);
        heap.push_back(Range{middle, range.end});
        std::push_heap(heap.begin(), heap.end(), by_size);
        ++subdomains;
    }
    std::size_t next_id = 0;
    for (const Range& range : heap) {
        for (std::size_t i = range.begin; i < range.end; ++i) {
            split.node_subdomain[order[i]] = static_cast<std::uint32_t>(next_id);
        }
        ++next_id;
    }
    split.count = next_id;
    return split;
}

double subdomain_capture_ratio(const CsrMatrix& matrix, const std::vector<std::uint32_t>& dof_subdomain) {
    if (matrix.nonzeros() == 0) return 0.0;
    std::size_t interior = 0;
    const auto& offsets = matrix.row_offsets();
    const auto& columns = matrix.column_indices();
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        const std::uint32_t owner = dof_subdomain[row];
        for (std::size_t index = offsets[row]; index < offsets[row + 1]; ++index) {
            if (dof_subdomain[columns[index]] == owner) ++interior;
        }
    }
    return static_cast<double>(interior) / static_cast<double>(matrix.nonzeros());
}

SubdomainBlockJacobi::SubdomainBlockJacobi(const CsrMatrix& matrix,
                                           std::vector<std::uint32_t> dof_subdomain,
                                           const std::size_t subdomain_count)
    : dof_count_(matrix.rows()), subdomain_count_(subdomain_count) {
    const auto setup_start = std::chrono::steady_clock::now();
    if (dof_subdomain.size() != matrix.rows()) {
        throw std::invalid_argument("subdomain assignment size does not match the matrix");
    }
    stats_.subdomains = subdomain_count_;
    stats_.capture_ratio =
        subdomain_count_ == 0 ? 0.0 : subdomain_capture_ratio(matrix, dof_subdomain);
    if (subdomain_count_ == 0) return;

    std::vector<std::size_t> sizes(subdomain_count_, 0);
    for (const std::uint32_t owner : dof_subdomain) {
        if (owner >= subdomain_count_) throw std::invalid_argument("subdomain id out of range");
        ++sizes[owner];
    }

    dof_begin_.assign(subdomain_count_ + 1, 0);
    for (std::size_t s = 0; s < subdomain_count_; ++s) {
        dof_begin_[s + 1] = dof_begin_[s] + sizes[s];
    }
    dofs_.resize(dof_begin_.back());

    // Visiting the DOFs in matrix order keeps each subdomain's local unknowns in
    // increasing global order, which makes the gather in apply() sequential.
    std::vector<std::size_t> cursor(dof_begin_.begin(), dof_begin_.end() - 1);
    for (std::size_t dof = 0; dof < dof_count_; ++dof) {
        const std::uint32_t owner = dof_subdomain[dof];
        dofs_[cursor[owner]++] = static_cast<std::uint32_t>(dof);
    }
    // Slot of a global DOF within its subdomain, so the dense block can be filled
    // without a search: local_index[dof] - dof_begin_[owner] is the local column.
    std::vector<std::size_t> local_index(dof_count_, 0);
    for (std::size_t s = 0; s < subdomain_count_; ++s) {
        for (std::size_t i = dof_begin_[s]; i < dof_begin_[s + 1]; ++i) local_index[dofs_[i]] = i;
    }

    factor_begin_.assign(subdomain_count_ + 1, 0);
    std::vector<std::size_t> block_offsets(subdomain_count_, 0);
    for (std::size_t s = 0; s < subdomain_count_; ++s) {
        const std::size_t size = dof_begin_[s + 1] - dof_begin_[s];
        factor_begin_[s + 1] = factor_begin_[s] + size * (size + 1) / 2;
        block_offsets[s] = size;
    }
    factor_.assign(factor_begin_.back(), 0.0);
    for (const std::size_t size : block_offsets) stats_.dense_dofs += size;

    const auto& offsets = matrix.row_offsets();
    const auto& columns = matrix.column_indices();
    const auto& values = matrix.values();
    std::vector<char> degenerate(subdomain_count_, 0);

    // Subdomains are disjoint, so every block can be factorized independently.
    // On a 200x200 mesh this loop is a few hundred small dense factorizations,
    // which is enough work to be worth spreading over the available cores.
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
    for (long long subdomain_index = 0; subdomain_index < static_cast<long long>(subdomain_count_);
         ++subdomain_index) {
        const std::size_t s = static_cast<std::size_t>(subdomain_index);
        const std::size_t begin = dof_begin_[s];
        const std::size_t size = dof_begin_[s + 1] - begin;
        if (size == 0) continue;
        std::vector<double> dense(size * size, 0.0);
        for (std::size_t i = 0; i < size; ++i) {
            const std::size_t dof = dofs_[begin + i];
            for (std::size_t index = offsets[dof]; index < offsets[dof + 1]; ++index) {
                const std::size_t column = columns[index];
                if (dof_subdomain[column] != s) continue;
                dense[i * size + (local_index[column] - begin)] = values[index];
            }
        }

        // Unconstrained diagonal entries are strictly positive for an assembled
        // stiffness matrix with Dirichlet data applied, so in-place Cholesky is
        // enough. A singular block would mean the model is a mechanism; rather
        // than abort, that subdomain falls back to an identity so the outer
        // Krylov iteration still makes progress and reports the residual.
        bool factorized = true;
        for (std::size_t i = 0; i < size && factorized; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double sum = dense[i * size + j];
                for (std::size_t k = 0; k < j; ++k) sum -= dense[i * size + k] * dense[j * size + k];
                if (i == j) {
                    if (!(sum > 0.0) || !std::isfinite(sum)) {
                        factorized = false;
                        break;
                    }
                    dense[i * size + i] = std::sqrt(sum);
                } else {
                    dense[i * size + j] = sum / dense[j * size + j];
                }
            }
        }
        if (!factorized) {
            degenerate[s] = 1;
            for (std::size_t i = 0; i < size; ++i) {
                for (std::size_t j = 0; j <= i; ++j) dense[i * size + j] = (i == j) ? 1.0 : 0.0;
            }
        }
        double* target = &factor_[factor_begin_[s]];
        for (std::size_t i = 0; i < size; ++i) {
            for (std::size_t j = 0; j <= i; ++j) target[i * (i + 1) / 2 + j] = dense[i * size + j];
            // Store the reciprocal of the pivot in place. Both triangular sweeps
            // divide by it on every step, and a dependent chain of divisions is
            // what makes a small dense block solve expensive.
            target[i * (i + 1) / 2 + i] = 1.0 / target[i * (i + 1) / 2 + i];
        }
    }
    for (const char flag : degenerate) stats_.degenerate_subdomains += flag == 1 ? 1 : 0;

    stats_.factor_nonzeros = factor_.size();
    std::size_t largest = 0;
    for (const std::size_t size : block_offsets) largest = std::max(largest, size);
    largest_block_ = largest;
    stats_.setup_ms = elapsed_ms(setup_start);
}

void SubdomainBlockJacobi::apply(const std::vector<double>& input, std::vector<double>& output) const {
    if (input.size() != dof_count_) throw std::invalid_argument("preconditioner input size mismatch");
    if (output.size() != dof_count_) output.assign(dof_count_, 0.0);

    // Each thread needs its own scratch for the triangular solves; subdomains are
    // disjoint, so the gather, the solve and the scatter never touch the same
    // output entry twice.
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> work(largest_block_);
#ifdef LIGHTFEM_HAS_OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
        for (long long subdomain_index = 0; subdomain_index < static_cast<long long>(subdomain_count_);
             ++subdomain_index) {
            const std::size_t s = static_cast<std::size_t>(subdomain_index);
            const std::size_t begin = dof_begin_[s];
            const std::size_t size = dof_begin_[s + 1] - begin;
            if (size == 0) continue;
            const double* factor = &factor_[factor_begin_[s]];
            double* scratch = work.data();
            for (std::size_t i = 0; i < size; ++i) scratch[i] = input[dofs_[begin + i]];
            // Forward substitution with L, then back substitution with L^T. Both
            // use the same packed lower factor, so only one copy is stored, and
            // the diagonal holds the reciprocal of the pivot. Row starts advance
            // by i+1 going forward and by k going backward, so no index is
            // recomputed with a multiply.
            std::size_t base = 0;
            for (std::size_t i = 0; i < size; ++i) {
                double sum = scratch[i];
                for (std::size_t k = 0; k < i; ++k) sum -= factor[base + k] * scratch[k];
                scratch[i] = sum * factor[base + i];
                base += i + 1;
            }
            std::size_t row = (size - 1) * size / 2;
            for (std::size_t i = size; i-- > 0;) {
                double sum = scratch[i];
                std::size_t start = row;
                for (std::size_t k = size - 1; k > i; --k) {
                    sum -= factor[start + i] * scratch[k];
                    start -= k;
                }
                scratch[i] = sum * factor[i * (i + 1) / 2 + i];
            }
            for (std::size_t i = 0; i < size; ++i) output[dofs_[begin + i]] = scratch[i];
        }
    }
}

void SubdomainBlockJacobi::account_apply(const double milliseconds) const {
    stats_.apply_ms += milliseconds;
    ++stats_.apply_calls;
}

IncompleteCholesky0::IncompleteCholesky0(const CsrMatrix& matrix)
    : rows_(matrix.rows()), inverse_diagonal_(matrix.rows(), 1.0) {
    const auto setup_start = std::chrono::steady_clock::now();
    const auto& offsets = matrix.row_offsets();
    const auto& matrix_columns = matrix.column_indices();
    const std::vector<double>& matrix_values = matrix.values();

    // Keep only the lower triangle: the upper part is never read, by either the
    // factorization or the two triangular sweeps.
    row_offsets_.assign(rows_ + 1, 0);
    columns_.reserve(matrix.nonzeros() / 2 + rows_);
    values_.reserve(matrix.nonzeros() / 2 + rows_);
    std::vector<std::size_t> diagonal_index(rows_, 0);
    for (std::size_t row = 0; row < rows_; ++row) {
        bool found = false;
        for (std::size_t index = offsets[row]; index < offsets[row + 1]; ++index) {
            const std::size_t column = matrix_columns[index];
            if (column > row) break;
            columns_.push_back(static_cast<std::uint32_t>(column));
            values_.push_back(matrix_values[index]);
            if (column == row) {
                diagonal_index[row] = values_.size() - 1;
                found = true;
            }
        }
        if (!found) {
            // A row with no stored diagonal cannot be factorized; append one so
            // the sweeps still have a defined pivot.
            columns_.push_back(static_cast<std::uint32_t>(row));
            values_.push_back(1.0);
            diagonal_index[row] = values_.size() - 1;
        }
        row_offsets_[row + 1] = values_.size();
    }
    stats_.factor_entries = values_.size();

    double diagonal_scale = 0.0;
    for (std::size_t row = 0; row < rows_; ++row) {
        diagonal_scale = std::max(diagonal_scale, std::abs(values_[diagonal_index[row]]));
    }
    const double floor_value = std::max(1.0e-30, 1.0e-14 * diagonal_scale);

    // Row-wise factorization. Because the matrix is stored symmetrically, the
    // product L[i][k] * L[j][k] can be accumulated by walking the two rows in
    // lock step up to the diagonal; both are sorted, so the intersection costs
    // one pass over the shorter row instead of a search.
    for (std::size_t row = 0; row < rows_; ++row) {
        const std::size_t begin = row_offsets_[row];
        const std::size_t diagonal = diagonal_index[row];
        for (std::size_t index = begin; index < diagonal; ++index) {
            const std::size_t column = columns_[index];
            const std::size_t other_begin = row_offsets_[column];
            const std::size_t other_diagonal = diagonal_index[column];
            double sum = values_[index];
            std::size_t left = begin;
            std::size_t right = other_begin;
            while (left < index && right < other_diagonal) {
                const std::size_t left_column = columns_[left];
                const std::size_t right_column = columns_[right];
                if (left_column < right_column) {
                    ++left;
                } else if (left_column > right_column) {
                    ++right;
                } else {
                    sum -= values_[left] * values_[right];
                    ++left;
                    ++right;
                }
            }
            values_[index] = sum * inverse_diagonal_[column];
        }
        double sum = values_[diagonal];
        for (std::size_t index = begin; index < diagonal; ++index) {
            sum -= values_[index] * values_[index];
        }
        if (!(sum > 0.0) || !std::isfinite(sum)) {
            // Modified IC: clamping the pivot keeps the factorization usable and
            // records that the pattern was too coarse to stay definite.
            ++stats_.breakdowns;
            sum = floor_value;
        }
        const double pivot = std::sqrt(sum);
        values_[diagonal] = pivot;
        inverse_diagonal_[row] = 1.0 / pivot;
    }
    stats_.setup_ms = elapsed_ms(setup_start);
}

void IncompleteCholesky0::apply(const std::vector<double>& input,
                                std::vector<double>& output) const {
    if (input.size() != rows_) throw std::invalid_argument("IC(0) input size mismatch");
    if (output.size() != rows_) output.assign(rows_, 0.0);
    // Forward substitution gathers over the row of L; back substitution
    // scatters with the same row, so no transposed copy of the factor is needed.
    std::vector<double>& work = output;
    work = input;
    for (std::size_t row = 0; row < rows_; ++row) {
        const std::size_t begin = row_offsets_[row];
        const std::size_t end = row_offsets_[row + 1];
        double sum = work[row];
        for (std::size_t index = begin; index + 1 < end; ++index) {
            sum -= values_[index] * work[columns_[index]];
        }
        work[row] = sum * inverse_diagonal_[row];
    }
    for (std::size_t row = rows_; row-- > 0;) {
        const std::size_t begin = row_offsets_[row];
        const std::size_t end = row_offsets_[row + 1];
        work[row] *= inverse_diagonal_[row];
        const double value = work[row];
        for (std::size_t index = begin; index + 1 < end; ++index) {
            work[columns_[index]] -= values_[index] * value;
        }
    }
}

void IncompleteCholesky0::account_apply(const double milliseconds) const {
    stats_.apply_ms += milliseconds;
    ++stats_.apply_calls;
}

}  // namespace lightfem
