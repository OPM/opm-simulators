/*
  Copyright 2026 NORCE Norwegian Research Centre AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_UPWIND_TRANSPORT_SOLVER_HEADER_INCLUDED
#define OPM_UPWIND_TRANSPORT_SOLVER_HEADER_INCLUDED

#include <opm/common/utility/SparseTable.hpp>

#include <opm/simulators/utils/Tarjan.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/// \file
///
/// Solver for the linear system of an implicit upwind transport scheme, cell by cell in the order
/// of the flow (a "reordering" solver).
///
/// Opm::solveReordered() solves `A x_q = b_q` for an assembled sparse matrix `A` with square N x N
/// blocks, e.g. a `Dune::BCRSMatrix`, and all right-hand sides together. Row `r` depends on the
/// cells `c` with a nonzero off-diagonal block `A[r][c]`, i.e. its upwind neighbours; blocks that
/// are stored but zero are ignored.
///
/// The cells are ordered with Opm::tarjan(), so every component only needs values that are already
/// final. A single cell is solved with its diagonal block; the cells of a flow loop depend on each
/// other and are solved directly if the loop is small, by block Gauss-Seidel otherwise.
///
/// For an M-matrix, as from an upwind scheme with storage, the solution is non-negative for
/// non-negative `rhs` and Gauss-Seidel converges. Reordering pays off when the matrix is nearly
/// triangular, as for pure advection; with diffusion or conduction the whole domain becomes one
/// loop.

namespace Opm
{

/// \brief Settings of solveReordered().
struct UpwindTransportSettings
{
    /// Largest flow loop, in cells, that is solved directly instead of by Gauss-Seidel.
    int maxDirectSize = 100;

    /// Gauss-Seidel stops when the largest change in a sweep is at most this times the largest
    /// value, per right-hand side and block component.
    double relTolerance = 1e-4;

    /// Gauss-Seidel also stops when the largest change in a sweep is at most this, so that a
    /// solution that goes to zero converges.
    double absTolerance = 1e-20;

    /// Most Gauss-Seidel sweeps over a flow loop before giving up.
    int maxIterations = 100;
};

/// \brief Statistics of a solveReordered() call.
struct UpwindTransportStatistics
{
    /// Number of components, i.e. cells or flow loops.
    int numComponents = 0;

    /// Number of components with more than one cell, i.e. flow loops.
    int numLoops = 0;

    /// Number of cells in the largest component.
    int largestComponent = 0;

    /// Most sweeps used for any flow loop that was solved by iteration.
    int maxIterations = 0;
};

namespace detail
{

    /// \brief Solve a dense system `a x = b` by Gaussian elimination with partial pivoting.
    ///
    /// \tparam Scalar Floating-point type.
    /// \param a      The m x m matrix, row-major. Destroyed.
    /// \param b      The m x numRhs right-hand sides, row-major. Replaced by the solution.
    /// \param m      Number of rows and columns of `a`.
    /// \param numRhs Number of right-hand sides, i.e. columns of `b`.
    /// \throws std::runtime_error if a pivot is exactly zero.
    template <class Scalar>
    void solveDenseSystem(std::vector<Scalar>& a,
                          std::vector<Scalar>& b,
                          std::size_t m,
                          std::size_t numRhs)
    {
        for (std::size_t col = 0; col < m; ++col) {
            // Gaussian elimination with partial pivoting.
            std::size_t pivot = col;
            for (std::size_t row = col + 1; row < m; ++row) {
                if (std::abs(a[row * m + col]) > std::abs(a[pivot * m + col])) {
                    pivot = row;
                }
            }
            if (a[pivot * m + col] == Scalar {0}) {
                throw std::runtime_error("Opm::solveReordered: singular matrix");
            }
            if (pivot != col) {
                for (std::size_t k = 0; k < m; ++k) {
                    std::swap(a[pivot * m + k], a[col * m + k]);
                }
                for (std::size_t k = 0; k < numRhs; ++k) {
                    std::swap(b[pivot * numRhs + k], b[col * numRhs + k]);
                }
            }
            for (std::size_t row = col + 1; row < m; ++row) {
                const Scalar factor = a[row * m + col] / a[col * m + col];
                if (factor == Scalar {0}) {
                    continue;
                }
                for (std::size_t k = col; k < m; ++k) {
                    a[row * m + k] -= factor * a[col * m + k];
                }
                for (std::size_t k = 0; k < numRhs; ++k) {
                    b[row * numRhs + k] -= factor * b[col * numRhs + k];
                }
            }
        }
        for (std::size_t col = m; col-- > 0;) {
            for (std::size_t k = 0; k < numRhs; ++k) {
                Scalar sum = b[col * numRhs + k];
                for (std::size_t j = col + 1; j < m; ++j) {
                    sum -= a[col * m + j] * b[j * numRhs + k];
                }
                b[col * numRhs + k] = sum / a[col * m + col];
            }
        }
    }

    /// \brief Whether any entry of a block is nonzero. NaN counts as nonzero, so that it is
    ///        rejected by MatrixCoefficients.
    ///
    /// \tparam Block Dense block with `rows`, `cols` and `block[a][b]`, e.g. `Dune::FieldMatrix`.
    /// \param block The block.
    /// \return True if an entry is not zero.
    template <class Block>
    bool isNonzeroBlock(const Block& block)
    {
        for (std::size_t a = 0; a < Block::rows; ++a) {
            for (std::size_t b = 0; b < Block::cols; ++b) {
                if (block[a][b] != 0) {
                    return true;
                }
            }
        }
        return false;
    }

    /// \brief Whether all entries of a block are finite.
    ///
    /// \tparam Block Dense block with `rows`, `cols` and `block[a][b]`, e.g. `Dune::FieldMatrix`.
    /// \param block The block.
    /// \return False if an entry is NaN or infinite.
    template <class Block>
    bool isFiniteBlock(const Block& block)
    {
        for (std::size_t a = 0; a < Block::rows; ++a) {
            for (std::size_t b = 0; b < Block::cols; ++b) {
                if (!std::isfinite(block[a][b])) {
                    return false;
                }
            }
        }
        return true;
    }

    /// \brief Call a function for every nonzero off-diagonal block of a sparse block matrix.
    ///
    /// \tparam Matrix   Sparse block matrix such as `Dune::BCRSMatrix`.
    /// \tparam Function Callable as `f(std::size_t row, std::size_t column, const Block& block)`.
    /// \param matrix The matrix.
    /// \param f      Called row by row, in the order of the entries in each row.
    template <class Matrix, class Function>
    void forEachNonzeroCoupling(const Matrix& matrix, Function&& f)
    {
        for (auto row = matrix.begin(); row != matrix.end(); ++row) {
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                if (entry.index() != row.index() && isNonzeroBlock(*entry)) {
                    f(row.index(), entry.index(), *entry);
                }
            }
        }
    }

    /// \brief The coefficients of solveReordered(), from a sparse block matrix.
    ///
    /// Only the nonzero off-diagonal blocks couple cells. The matrix must outlive this object.
    ///
    /// The interface that ReorderedSystem uses:
    ///  - `Scalar` and `blockSize`, the size N of the blocks.
    ///  - `numCells()`, and `graph()` with the cells each cell depends on.
    ///  - `diagonal(c)`, a reference to the N x N diagonal block.
    ///  - `forEachCoupling(c, f)`, which calls `f(u, block)` for each cell `u` in row `c` of
    ///    `graph()` with the off-diagonal block of the equation.
    ///
    /// \tparam Matrix Square sparse block matrix such as `Dune::BCRSMatrix` with square blocks.
    template <class Matrix>
    class MatrixCoefficients
    {
        using Block = typename Matrix::block_type;
        static_assert(Block::rows == Block::cols, "The matrix blocks must be square");

    public:
        using Scalar = typename Block::field_type;
        static constexpr std::size_t blockSize = Block::rows;

        /// \brief Collect the diagonal blocks and the nonzero couplings of a matrix.
        ///
        /// \param matrix The matrix. Every row must have a diagonal block.
        /// \throws std::invalid_argument if the matrix is not square or a row has no diagonal
        ///         block.
        /// \throws std::runtime_error if a diagonal or coupling block contains NaN or infinite
        ///         values.
        explicit MatrixCoefficients(const Matrix& matrix)
        {
            if (matrix.N() != matrix.M()) {
                throw std::invalid_argument("Opm::solveReordered: the matrix must be square");
            }
            diagonal_.reserve(matrix.N());
            for (std::size_t r = 0; r < matrix.N(); ++r) {
                const auto entry = matrix[r].find(r);
                if (entry == matrix[r].end()) {
                    throw std::invalid_argument("Opm::solveReordered: row " + std::to_string(r)
                                                + " has no diagonal block");
                }
                // Checked up front, since an infinite entry can give a finite solution, e.g.
                // b / inf = 0 or after elimination in a directly solved loop
                if (!isFiniteBlock(*entry)) {
                    throw std::runtime_error("Opm::solveReordered: NaN or infinite entry in the "
                                             "diagonal block of row " + std::to_string(r));
                }
                diagonal_.push_back(&*entry);
            }

            std::vector<int> rowStart(matrix.N() + 1, 0);
            std::vector<int> columns;
            std::vector<const Block*> blocks;
            forEachNonzeroCoupling(
                matrix, [&](std::size_t r, std::size_t c, const Block& block) {
                    if (!isFiniteBlock(block)) {
                        throw std::runtime_error("Opm::solveReordered: NaN or infinite entry in "
                                                 "the block of row " + std::to_string(r)
                                                 + ", column " + std::to_string(c));
                    }
                    ++rowStart[r + 1];
                    columns.push_back(static_cast<int>(c));
                    blocks.push_back(&block);
                });
            for (std::size_t r = 0; r < matrix.N(); ++r) {
                rowStart[r + 1] += rowStart[r];
            }
            graph_ = SparseTable<int>(std::move(columns), std::vector<int>(rowStart));
            couplings_ = SparseTable<const Block*>(std::move(blocks), std::move(rowStart));
        }

        /// \brief The number of cells, i.e. block rows of the matrix.
        ///
        /// \return The number of cells.
        std::size_t numCells() const
        {
            return diagonal_.size();
        }

        /// \brief The dependency graph, as input to Opm::tarjan().
        ///
        /// \return One row per cell `r`, listing the columns `c != r` with a nonzero block
        ///         `A[r][c]`, in the order of the matrix row. Blocks that are stored but zero give
        ///         no edge.
        const SparseTable<int>& graph() const
        {
            return graph_;
        }

        /// \brief The diagonal block of a cell.
        ///
        /// \param c The cell.
        /// \return A reference to the block `A[c][c]` in the matrix.
        const Block& diagonal(std::size_t c) const
        {
            return *diagonal_[c];
        }

        /// \brief Call a function for each cell that a cell depends on.
        ///
        /// \tparam Function Callable as `f(std::size_t u, const Block& block)`.
        /// \param c The cell.
        /// \param f Called for each cell `u` in row `c` of graph(), with the block `A[c][u]`.
        template <class Function>
        void forEachCoupling(std::size_t c, Function&& f) const
        {
            auto block = couplings_[static_cast<int>(c)].begin();
            for (const int u : graph_[static_cast<int>(c)]) {
                f(static_cast<std::size_t>(u), **block++);
            }
        }

    private:
        std::vector<const Block*> diagonal_; //!< Diagonal block of each cell
        SparseTable<int> graph_; //!< Cells that each cell depends on
        SparseTable<const Block*> couplings_; //!< Blocks of graph(), in the same layout
    };


    /// \brief The system of solveReordered() and the work arrays for solving it.
    ///
    /// \tparam Coefficients Source of the coefficients, with the interface of MatrixCoefficients.
    /// \tparam Vectors      Container of block vectors, e.g. `std::vector<Dune::BlockVector<...>>`,
    ///                      with `x[q][c][a]` component `a` of cell `c` for right-hand side `q`.
    template <class Coefficients, class Vectors>
    class ReorderedSystem
    {
        using Scalar = typename Coefficients::Scalar;
        static constexpr std::size_t N = Coefficients::blockSize;

    public:
        /// \brief Set up the system.
        ///
        /// \param coefficients The matrix.
        /// \param rhs          One vector per right-hand side.
        /// \param x            The solution, in the layout of `rhs`. The initial guess for loops
        ///                     solved by Gauss-Seidel.
        /// \param settings     Direct-solve limit, tolerances and iteration limit.
        /// \param name         Name of the calling function, for error messages.
        /// \throws std::invalid_argument if the number or the lengths of the vectors do not match.
        ReorderedSystem(const Coefficients& coefficients,
                        const Vectors& rhs,
                        Vectors& x,
                        const UpwindTransportSettings& settings,
                        const std::string& name)
            : coefficients_(coefficients)
            , rhs_(rhs)
            , x_(x)
            , settings_(settings)
            , name_(name)
            , numCells_(coefficients.numCells())
            , numVectors_(rhs.size())
            , localIndex_(numCells_, -1)
            , cellSolution_(N * numVectors_)
            , change_(N * numVectors_)
            , scale_(N * numVectors_)
        {
            validate_();
        }

        /// \brief Solve the components in the order of Opm::tarjan().
        ///
        /// Single cells are solved with their diagonal block, loops of at most
        /// UpwindTransportSettings::maxDirectSize cells directly, larger loops by Gauss-Seidel.
        ///
        /// \return Statistics about the components.
        /// \throws std::runtime_error if a block or loop is singular, a solution value is not
        ///         finite, or Gauss-Seidel does not converge.
        UpwindTransportStatistics solve()
        {
            const auto components = tarjan(coefficients_.graph());
            UpwindTransportStatistics statistics;
            statistics.numComponents = components.numComponents();
            for (const auto cells : components.components()) {
                const int size = static_cast<int>(cells.size());
                statistics.largestComponent = std::max(statistics.largestComponent, size);
                if (size == 1) {
                    solveCell_(cells);
                } else {
                    ++statistics.numLoops;
                    if (size <= settings_.maxDirectSize) {
                        solveDirectly_(cells);
                    } else {
                        const int sweeps = solveByGaussSeidel_(cells);
                        statistics.maxIterations = std::max(statistics.maxIterations, sweeps);
                    }
                }
            }
            return statistics;
        }

    private:
        /// \brief Check the number and the lengths of the vectors.
        ///
        /// \throws std::invalid_argument if they do not match the matrix and `rhs`.
        void validate_() const
        {
            if (x_.size() != numVectors_) {
                throw std::invalid_argument(name_ + ": sizes do not match");
            }
            for (std::size_t q = 0; q < numVectors_; ++q) {
                if (static_cast<std::size_t>(rhs_[q].size()) != numCells_
                    || static_cast<std::size_t>(x_[q].size()) != numCells_) {
                    throw std::invalid_argument(name_ + ": wrong vector length");
                }
            }
        }

        /// \brief Check a value of the solution before it is stored. This catches NaN and infinite
        ///        values in the right-hand side and the initial values, and overflow.
        ///
        /// \param value     The computed value.
        /// \param q         Its right-hand side.
        /// \param cells     The component it belongs to, for the error message.
        /// \param useInLoop Whether the values of `cells` are inputs (Gauss-Seidel) or unknowns
        ///                  (direct solve).
        /// \throws std::runtime_error with the cause if the value is NaN or infinite.
        void checkFinite_(Scalar value,
                          std::size_t q,
                          std::span<const int> cells,
                          bool useInLoop) const
        {
            if (!std::isfinite(value)) [[unlikely]] {
                throw std::runtime_error(nonFiniteMessage_(q, cells, useInLoop));
            }
        }

        /// \brief The error message for a NaN or infinite solution.
        ///
        /// Only called on failure. Looks for the first input of the component that is not finite.
        ///
        /// \param q         The right-hand side.
        /// \param cells     The component.
        /// \param useInLoop See checkFinite_().
        /// \return The message, naming the right-hand side or initial value that is not finite, or
        ///         an overflow if all of them are finite.
        std::string nonFiniteMessage_(std::size_t q,
                                      std::span<const int> cells,
                                      bool useInLoop) const
        {
            // A sorted copy, so that looking up the cells of a large loop stays cheap
            std::vector<int> sorted(cells.begin(), cells.end());
            std::ranges::sort(sorted);
            const auto inCells = [&sorted](std::size_t u) {
                return std::ranges::binary_search(sorted, static_cast<int>(u));
            };
            std::string cause;
            for (const int cell : cells) {
                const auto c = static_cast<std::size_t>(cell);
                for (std::size_t a = 0; a < N && cause.empty(); ++a) {
                    if (!std::isfinite(rhs_[q][c][a])) {
                        cause = "the right-hand side in row " + std::to_string(c);
                    }
                }
                coefficients_.forEachCoupling(c, [&](std::size_t u, const auto&) {
                    if (!cause.empty()) {
                        return;
                    }
                    if (useInLoop || !inCells(u)) {
                        for (std::size_t b = 0; b < N; ++b) {
                            if (!std::isfinite(x_[q][u][b])) {
                                cause = "the value in row " + std::to_string(u)
                                    + (inCells(u) ? ", an initial value of Gauss-Seidel" : "");
                                return;
                            }
                        }
                    }
                });
                if (!cause.empty()) {
                    break;
                }
            }
            if (cause.empty()) {
                cause = "an overflow in the solve";
            }
            return name_ + ": NaN or infinite solution for right-hand side " + std::to_string(q)
                + " in a component with row " + std::to_string(cells.front()) + ", caused by "
                + cause;
        }

        /// \brief The right-hand side of a cell with the known couplings moved over.
        ///
        /// Couplings to the cells of a loop that is being solved directly (`localIndex_ >= 0`) are
        /// left out, since those cells are unknowns of the dense system.
        ///
        /// \param q The right-hand side.
        /// \param c The cell.
        /// \return `rhs` minus the couplings times the values of the cells that `c` depends on.
        std::array<Scalar, N> rightHandSide_(std::size_t q, std::size_t c) const
        {
            std::array<Scalar, N> sum;
            for (std::size_t a = 0; a < N; ++a) {
                sum[a] = rhs_[q][c][a];
            }
            coefficients_.forEachCoupling(c, [&](std::size_t u, const auto& block) {
                if (localIndex_[u] < 0) {
                    for (std::size_t a = 0; a < N; ++a) {
                        for (std::size_t b = 0; b < N; ++b) {
                            sum[a] -= block[a][b] * x_[q][u][b];
                        }
                    }
                }
            });
            return sum;
        }

        /// \brief Solve with the diagonal block of a cell, for all right-hand sides at once.
        ///
        /// The result is left in `cellSolution_`, N rows by numVectors columns.
        ///
        /// \param c The cell.
        /// \throws std::runtime_error if the diagonal block is singular.
        void solveDiagonal_(std::size_t c)
        {
            if constexpr (N == 1) {
                // One division
                const Scalar diagonal = coefficients_.diagonal(c)[0][0];
                if (diagonal == Scalar {0}) {
                    throw std::runtime_error(name_ + ": singular diagonal in row "
                                             + std::to_string(c));
                }
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    cellSolution_[q] = rightHandSide_(q, c)[0] / diagonal;
                }
            } else {
                matrix_.resize(N * N);
                const auto& diagonal = coefficients_.diagonal(c);
                for (std::size_t a = 0; a < N; ++a) {
                    for (std::size_t b = 0; b < N; ++b) {
                        matrix_[a * N + b] = diagonal[a][b];
                    }
                }
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    const auto sum = rightHandSide_(q, c);
                    for (std::size_t a = 0; a < N; ++a) {
                        cellSolution_[a * numVectors_ + q] = sum[a];
                    }
                }
                try {
                    solveDenseSystem(matrix_, cellSolution_, N, numVectors_);
                } catch (const std::runtime_error&) {
                    throw std::runtime_error(name_ + ": singular diagonal block in row "
                                             + std::to_string(c));
                }
            }
        }

        /// \brief Solve a cell that is not in a flow loop, with its diagonal block.
        ///
        /// \param cells The component, with one cell.
        /// \throws std::runtime_error if the block is singular or the solution is not finite.
        void solveCell_(std::span<const int> cells)
        {
            const auto c = static_cast<std::size_t>(cells[0]);
            solveDiagonal_(c);
            for (std::size_t q = 0; q < numVectors_; ++q) {
                for (std::size_t a = 0; a < N; ++a) {
                    const Scalar value = cellSolution_[a * numVectors_ + q];
                    checkFinite_(value, q, cells, false);
                    x_[q][c][a] = value;
                }
            }
        }

        /// \brief Solve a small flow loop directly, for all of its cells at once.
        ///
        /// \param cells The component.
        /// \throws std::runtime_error if the loop matrix is singular or the solution is not finite.
        void solveDirectly_(std::span<const int> cells)
        {
            const std::size_t m = cells.size();
            const std::size_t size = m * N;
            matrix_.assign(size * size, Scalar {0});
            solution_.assign(size * numVectors_, Scalar {0});
            for (std::size_t i = 0; i < m; ++i) {
                localIndex_[cells[i]] = static_cast<int>(i);
            }
            for (std::size_t i = 0; i < m; ++i) {
                const auto c = static_cast<std::size_t>(cells[i]);
                const auto& diagonal = coefficients_.diagonal(c);
                for (std::size_t a = 0; a < N; ++a) {
                    for (std::size_t b = 0; b < N; ++b) {
                        matrix_[(i * N + a) * size + i * N + b] = diagonal[a][b];
                    }
                }
                coefficients_.forEachCoupling(c, [&](std::size_t u, const auto& block) {
                    if (localIndex_[u] >= 0) {
                        const auto j = static_cast<std::size_t>(localIndex_[u]);
                        for (std::size_t a = 0; a < N; ++a) {
                            for (std::size_t b = 0; b < N; ++b) {
                                matrix_[(i * N + a) * size + j * N + b] += block[a][b];
                            }
                        }
                    }
                });
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    // Cells outside the loop are already solved.
                    const auto sum = rightHandSide_(q, c);
                    for (std::size_t a = 0; a < N; ++a) {
                        solution_[(i * N + a) * numVectors_ + q] = sum[a];
                    }
                }
            }
            for (std::size_t i = 0; i < m; ++i) {
                localIndex_[cells[i]] = -1;
            }
            try {
                solveDenseSystem(matrix_, solution_, size, numVectors_);
            } catch (const std::runtime_error&) {
                throw std::runtime_error(name_ + ": singular matrix for a flow loop of "
                                         + std::to_string(m) + " cells");
            }
            for (std::size_t i = 0; i < m; ++i) {
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    for (std::size_t a = 0; a < N; ++a) {
                        const Scalar value = solution_[(i * N + a) * numVectors_ + q];
                        checkFinite_(value, q, cells, false);
                        x_[q][static_cast<std::size_t>(cells[i])][a] = value;
                    }
                }
            }
        }

        /// \brief Solve a large flow loop by Gauss-Seidel sweeps, in the order of its cells.
        ///
        /// \param cells The component.
        /// \return The number of sweeps.
        /// \throws std::runtime_error if it does not converge in
        ///         UpwindTransportSettings::maxIterations sweeps, or a value is not finite.
        int solveByGaussSeidel_(std::span<const int> cells)
        {
            for (int sweep = 1; sweep <= settings_.maxIterations; ++sweep) {
                std::fill(change_.begin(), change_.end(), Scalar {0});
                std::fill(scale_.begin(), scale_.end(), Scalar {0});
                for (const int cell : cells) {
                    const auto c = static_cast<std::size_t>(cell);
                    solveDiagonal_(c);
                    for (std::size_t q = 0; q < numVectors_; ++q) {
                        for (std::size_t a = 0; a < N; ++a) {
                            auto& value = x_[q][c][a];
                            const Scalar updated = cellSolution_[a * numVectors_ + q];
                            checkFinite_(updated, q, cells, true);
                            const std::size_t k = q * N + a;
                            change_[k] = std::max(change_[k], std::abs(updated - value));
                            scale_[k] = std::max(scale_[k], std::abs(updated));
                            value = updated;
                        }
                    }
                }
                // Each block component converges on its own scale, so small components are not
                // only converged relative to large ones. The absolute tolerance covers zero
                bool converged = true;
                for (std::size_t k = 0; k < change_.size(); ++k) {
                    converged = converged
                        && (change_[k] <= settings_.relTolerance * scale_[k]
                            || change_[k] <= settings_.absTolerance);
                }
                if (converged) {
                    return sweep;
                }
            }
            throw std::runtime_error(name_ + ": no convergence after "
                                     + std::to_string(settings_.maxIterations)
                                     + " sweeps for a flow loop of " + std::to_string(cells.size())
                                     + " cells");
        }

        const Coefficients& coefficients_; //!< The matrix
        const Vectors& rhs_; //!< The right-hand sides
        Vectors& x_; //!< The solution
        const UpwindTransportSettings& settings_; //!< The settings
        std::string name_; //!< Name of the calling function, for error messages
        std::size_t numCells_; //!< Number of cells
        std::size_t numVectors_; //!< Number of right-hand sides

        // Work arrays, reused for every cell and flow loop
        std::vector<int> localIndex_; //!< Index in the loop that is solved directly, else -1
        std::vector<Scalar> matrix_; //!< Dense diagonal block, or matrix of a directly solved loop
        std::vector<Scalar> solution_; //!< Right-hand sides and solution of a directly solved loop
        std::vector<Scalar> cellSolution_; //!< Right-hand sides and solution of a single cell
        std::vector<Scalar> change_; //!< Largest change in a sweep, per rhs and block component
        std::vector<Scalar> scale_; //!< Largest value in a sweep, per rhs and block component
    };

} // namespace detail

/// \brief Solve the block system `A x_q = b_q` described in the file comment.
///
/// \tparam Matrix  Square sparse block matrix with square N x N blocks, such as
///                 `Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>`.
/// \tparam Vectors Container of block vectors, e.g.
///                 `std::vector<Dune::BlockVector<Dune::FieldVector<double, N>>>`.
/// \param matrix   Square sparse matrix with square N x N blocks, such as a `Dune::BCRSMatrix`.
///                 Every row must have a diagonal block.
/// \param rhs      One block vector of length `matrix.N()` per right-hand side.
/// \param x        The solution, in the layout of `rhs`. The values on entry are the initial guess
///                 when a flow loop is solved by iteration; otherwise they are not used.
/// \param settings When flow loops are solved directly, and the tolerances and iteration limit of
///                 the Gauss-Seidel sweeps otherwise.
///
/// \return Statistics about the components that were found and solved.
///
/// \throws std::invalid_argument if the sizes do not match or a row has no diagonal block.
/// \throws std::runtime_error on NaN or infinite values, a singular block or loop, or if
///         Gauss-Seidel does not converge.
template <class Matrix, class Vectors>
UpwindTransportStatistics
solveReordered(const Matrix& matrix,
               const Vectors& rhs,
               Vectors& x,
               const UpwindTransportSettings& settings = {})
{
    const detail::MatrixCoefficients<Matrix> coefficients(matrix);
    return detail::ReorderedSystem<detail::MatrixCoefficients<Matrix>, Vectors>(
               coefficients, rhs, x, settings, "Opm::solveReordered")
        .solve();
}

} // namespace Opm

#endif // OPM_UPWIND_TRANSPORT_SOLVER_HEADER_INCLUDED
