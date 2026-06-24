/*
  Copyright 2024 SINTEF AS
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

#ifndef OPM_DILU2_HEADER_INCLUDED
#define OPM_DILU2_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/unused.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <cstddef>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif

namespace Dune
{

/*! \brief Multicolor variant of the OpenMP-parallel DILU preconditioner.
 *
 *  \details The original MultithreadDILU (DILU.hpp) parallelizes the *exact*
 *  natural-ordering DILU using a level-set (wavefront) schedule. On structured
 *  reservoir grids that schedule has a long critical path (e.g. SPE9: 62 thin
 *  levels), so the threaded triangular solve is dominated by the ~2*nlevels
 *  synchronization barriers per apply and is slower than serial.
 *
 *  MultithreadDILU2 instead reorders the unknowns by a multicolor graph coloring
 *  of the (assumed symmetric) sparsity pattern. Rows of the same color are
 *  mutually non-adjacent, hence independent, so the forward/backward solves and
 *  the factorization need only `#colors` levels (a handful for grid graphs)
 *  instead of one per wavefront. This trades exactness — it is the DILU of the
 *  color-permuted matrix, so convergence per iteration differs slightly — for
 *  far fewer barriers and real thread scaling.
 *
 *  Vectors stay in natural ordering; only the lower/upper split is taken in the
 *  color ordering via `natural_to_color_pos_`. The matrix is not copied.
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 */
template <class M, class X, class Y>
class MultithreadDILU2 : public PreconditionerWithUpdate<X, Y>
{
public:
    using matrix_type = M;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using block_type = typename M::block_type;

    explicit MultithreadDILU2(const M& A)
        : A_(A)
    {
        OPM_TIMEBLOCK(prec_construct);
        computeColoring();
        buildColumnSplit();
        Dinv_.resize(A_.N());
        update();
    }

    void update() override
    {
        OPM_TIMEBLOCK(prec_update);
        // Block diagonal initialization (embarrassingly parallel over all rows).
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (std::size_t i = 0; i < A_.N(); ++i) {
            Dinv_[i] = A_[i][i];
        }

        // Factorization: process colors in order; within a color rows are
        // independent (non-adjacent), so each color is one parallel sweep.
        // Uses the precomputed color-lower neighbor lists (with cached A[j,i]).
        for (std::size_t color = 0; color < color_sets_.size(); ++color) {
            const auto rows = color_sets_[color];
            const int n = rows.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int k = 0; k < n; ++k) {
                const std::size_t i = rows.begin()[k];
                auto Dinv_temp = Dinv_[i];
                for (std::size_t idx = lo_off_[i]; idx < lo_off_[i + 1]; ++idx) {
                    if (lo_blkT_[idx] != nullptr) {
                        // Dinv_temp -= A[i,j] * Dinv[j] * A[j,i]
                        Dinv_temp -= (*lo_blk_[idx]) * Dune::FieldMatrix(Dinv_[lo_col_[idx]]) * (*lo_blkT_[idx]);
                    }
                }
                Dinv_temp.invert();
                Dinv_[i] = Dinv_temp;
            }
        }
    }

    void pre(X& v, Y& d) override
    {
        DUNE_UNUSED_PARAMETER(v);
        DUNE_UNUSED_PARAMETER(d);
    }

    void apply(X& v, const Y& d) override
    {
        OPM_TIMEBLOCK(prec_apply);
        using Xblock = typename X::block_type;
        using Yblock = typename Y::block_type;

        // Lower triangular solve (forward): (D + L) y = d, y stored in v.
        // Uses the precomputed color-lower neighbor lists (no per-apply row scan
        // or color-position lookups).
        {
            OPM_TIMEBLOCK(lower_solve);
            for (std::size_t color = 0; color < color_sets_.size(); ++color) {
                const auto rows = color_sets_[color];
                const int n = rows.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int k = 0; k < n; ++k) {
                    const std::size_t i = rows.begin()[k];
                    Yblock rhs = d[i];
                    for (std::size_t idx = lo_off_[i]; idx < lo_off_[i + 1]; ++idx) {
                        // rhs -= A[i,j] * y[j]
                        lo_blk_[idx]->mmv(v[lo_col_[idx]], rhs);
                    }
                    // y_i = Dinv_i * rhs
                    Dinv_[i].mv(rhs, v[i]);
                }
            }
        }

        // Upper triangular solve (backward): (D + U) v = D y.
        {
            OPM_TIMEBLOCK(upper_solve);
            for (int color = static_cast<int>(color_sets_.size()) - 1; color >= 0; --color) {
                const auto rows = color_sets_[color];
                const int n = rows.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int k = 0; k < n; ++k) {
                    const std::size_t i = rows.begin()[k];
                    Xblock rhs(0.0);
                    for (std::size_t idx = up_off_[i]; idx < up_off_[i + 1]; ++idx) {
                        // rhs += A[i,j] * v[j]
                        up_blk_[idx]->umv(v[up_col_[idx]], rhs);
                    }
                    // v_i = y_i - Dinv_i * rhs   (v_i currently holds y_i)
                    Dinv_[i].mmv(rhs, v[i]);
                }
            }
        }
    }

    void post(X& x) override
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    std::vector<block_type> getDiagonal()
    {
        return Dinv_;
    }

    //! Number of colors (= number of parallel sweeps per triangular solve).
    std::size_t numColors() const
    {
        return color_sets_.size();
    }

    virtual SolverCategory::Category category() const override
    {
        return SolverCategory::sequential;
    }

    virtual bool hasPerfectUpdate() const override
    {
        return true;
    }

private:
    //! \brief The matrix we operate on.
    const M& A_;
    //! \brief The inverse of the (modified) block diagonal, indexed by natural row.
    std::vector<block_type> Dinv_;
    //! \brief Natural row indices grouped by color; rows in a group are independent.
    Opm::SparseTable<std::size_t> color_sets_;
    //! \brief Position of each natural row in the color ordering (defines L/U split).
    std::vector<std::size_t> natural_to_color_pos_;

    // Precomputed per-row neighbor split in color ordering (CSR-like, indexed by
    // natural row). "lower" = neighbors with smaller color position, "upper" =
    // larger. Block pointers point into A_ and stay valid across update() (only
    // values change, not sparsity), so apply()/update() avoid per-call scans,
    // color-position lookups and find() of the transpose entry.
    std::vector<std::size_t> lo_off_;             //!< size N+1, row offsets into lower arrays
    std::vector<std::size_t> lo_col_;             //!< lower neighbor column (natural)
    std::vector<const block_type*> lo_blk_;       //!< &A[i,j] for lower neighbors
    std::vector<const block_type*> lo_blkT_;      //!< &A[j,i] (transpose) or nullptr
    std::vector<std::size_t> up_off_;             //!< size N+1, row offsets into upper arrays
    std::vector<std::size_t> up_col_;             //!< upper neighbor column (natural)
    std::vector<const block_type*> up_blk_;       //!< &A[i,j] for upper neighbors

    //! \brief Precompute the color-lower / color-upper neighbor lists for each row.
    void buildColumnSplit()
    {
        const std::size_t N = A_.N();
        std::vector<std::size_t> loc(N, 0), upc(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t pos_i = natural_to_color_pos_[i];
            for (auto a_ij = A_[i].begin(); a_ij != A_[i].end(); ++a_ij) {
                const std::size_t j = a_ij.index();
                if (j == i) {
                    continue;
                }
                if (natural_to_color_pos_[j] < pos_i) {
                    ++loc[i];
                } else {
                    ++upc[i];
                }
            }
        }
        lo_off_.resize(N + 1);
        up_off_.resize(N + 1);
        lo_off_[0] = 0;
        up_off_[0] = 0;
        for (std::size_t i = 0; i < N; ++i) {
            lo_off_[i + 1] = lo_off_[i] + loc[i];
            up_off_[i + 1] = up_off_[i] + upc[i];
        }
        lo_col_.resize(lo_off_[N]);
        lo_blk_.resize(lo_off_[N]);
        lo_blkT_.resize(lo_off_[N]);
        up_col_.resize(up_off_[N]);
        up_blk_.resize(up_off_[N]);
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t pos_i = natural_to_color_pos_[i];
            std::size_t lp = lo_off_[i];
            std::size_t up = up_off_[i];
            for (auto a_ij = A_[i].begin(); a_ij != A_[i].end(); ++a_ij) {
                const std::size_t j = a_ij.index();
                if (j == i) {
                    continue;
                }
                if (natural_to_color_pos_[j] < pos_i) {
                    lo_col_[lp] = j;
                    lo_blk_[lp] = &(*a_ij);
                    const auto a_ji = A_[j].find(i);
                    lo_blkT_[lp] = (a_ji != A_[j].end()) ? &(*a_ji) : nullptr;
                    ++lp;
                } else {
                    up_col_[up] = j;
                    up_blk_[up] = &(*a_ij);
                    ++up;
                }
            }
        }
    }

    //! \brief Greedy multicolor coloring of the symmetric sparsity pattern.
    void computeColoring()
    {
        const std::size_t N = A_.N();
        std::vector<int> color(N, -1);
        // marker[c] == i means color c is used by a neighbor while coloring row i
        std::vector<std::size_t> marker;
        marker.reserve(16);
        int num_colors = 0;

        for (std::size_t i = 0; i < N; ++i) {
            for (auto a_ij = A_[i].begin(); a_ij != A_[i].end(); ++a_ij) {
                const std::size_t j = a_ij.index();
                if (j != i && color[j] >= 0) {
                    marker[color[j]] = i;
                }
            }
            int c = 0;
            while (c < num_colors && marker[c] == i) {
                ++c;
            }
            if (c == num_colors) {
                marker.push_back(N); // sentinel: not yet used in any iteration
                ++num_colors;
            }
            color[i] = c;
        }

        // Group rows by color into a SparseTable and record color-order positions.
        std::vector<std::size_t> count(num_colors, 0);
        for (std::size_t i = 0; i < N; ++i) {
            ++count[color[i]];
        }
        // Flatten: rows of color 0 first, then color 1, ...
        std::vector<std::size_t> offset(num_colors, 0);
        for (int c = 1; c < num_colors; ++c) {
            offset[c] = offset[c - 1] + count[c - 1];
        }
        std::vector<std::size_t> flat(N);
        std::vector<std::size_t> cursor = offset;
        natural_to_color_pos_.resize(N);
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t pos = cursor[color[i]]++;
            flat[pos] = i;
            natural_to_color_pos_[i] = pos;
        }
        color_sets_ = Opm::SparseTable<std::size_t>(
            flat.data(), flat.data() + flat.size(), count.data(), count.data() + count.size());
    }
};

} // namespace Dune

#endif // OPM_DILU2_HEADER_INCLUDED
