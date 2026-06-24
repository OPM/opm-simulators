/*
  Copyright 2026 SINTEF

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
#ifndef OPM_THREADED_SCALAR_PRODUCT_HEADER_INCLUDED
#define OPM_THREADED_SCALAR_PRODUCT_HEADER_INCLUDED

#include <dune/istl/scalarproducts.hh>

#include <cmath>
#include <cstddef>

namespace Opm {

/*!
   \brief OpenMP-threaded sequential scalar product.

   Drop-in replacement for Dune::SeqScalarProduct that parallelizes the dot
   reduction over block rows. Used by FlexibleSolver in the sequential case so
   the per-iteration dot/norm of the Krylov solver are not serialized.

   Note on floating point: an OpenMP reduction may sum the per-thread partial
   sums in a different order than the sequential loop, so results can differ in
   the last ULPs. This only affects convergence-history bit-reproducibility, not
   correctness; with OMP_NUM_THREADS=1 it is identical to Dune::SeqScalarProduct.
 */
template <class X>
class ThreadedSeqScalarProduct : public Dune::ScalarProduct<X>
{
public:
    using field_type = typename X::field_type;
    using real_type = typename Dune::ScalarProduct<X>::real_type;

    //! Sequential category: no inter-process communication.
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

    //! Below this number of block rows the OpenMP fork/join overhead of a
    //! parallel reduction exceeds the work of the dot product, so we fall back
    //! to the sequential loop. Measured on small black-oil systems (e.g. SPE9,
    //! ~9k blocks) where per-call threading of the cheap dot/norm regressed the
    //! solve; a persistent-parallel-region Krylov solver is the real fix for
    //! threading the cheap vector kernels.
    static constexpr std::size_t threading_threshold = 50000;

    field_type dot(const X& x, const X& y) const override
    {
        field_type sum = 0;
        const std::size_t n = x.size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) if (n >= threading_threshold)
#endif
        for (std::size_t i = 0; i < n; ++i) {
            sum += x[i] * y[i];
        }
        return sum;
    }

    real_type norm(const X& x) const override
    {
        using std::sqrt;
        return sqrt(static_cast<real_type>(this->dot(x, x)));
    }
};

} // namespace Opm

#endif // OPM_THREADED_SCALAR_PRODUCT_HEADER_INCLUDED
