/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.
  Copyright 2020 OPM-OP AS.

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


#ifndef OPM_WELLHELPERS_HEADER_INCLUDED
#define OPM_WELLHELPERS_HEADER_INCLUDED

#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/dynmatrix.hh>

#include <array>

namespace Opm {

class ParallelWellInfo;
struct WellProductionControls;
struct WellInjectionControls;
enum class WellProducerCMode;
enum class WellInjectorCMode;

namespace wellhelpers {

/// \brief A wrapper around the B matrix for distributed wells
///
/// For standard wells the B matrix, is basically a multiplication
/// of the equation of the perforated cells followed by a reduction
/// (summation) of these to the well equations.
///
/// This class does that in the functions mv and mmv (from the DUNE
/// matrix interface.
///
/// \tparam Scalar The scalar used for the computation.
template<typename Scalar>
class ParallelStandardWellB
{
public:
    using Block = Dune::DynamicMatrix<Scalar>;
    using Matrix = Dune::BCRSMatrix<Block>;

    ParallelStandardWellB(const Matrix& B, const ParallelWellInfo& parallel_well_info);

    //! y = A x
    template<class X, class Y>
    void mv (const X& x, Y& y) const;

    //! y = A x
    template<class X, class Y>
    void mmv (const X& x, Y& y) const;

private:
    const Matrix& B_;
    const ParallelWellInfo& parallel_well_info_;
};

double computeHydrostaticCorrection(const double well_ref_depth,
                                    const double vfp_ref_depth,
                                    const double rho, const double gravity);


/// \brief Sums entries of the diagonal Matrix for distributed wells
template<typename Scalar, typename Comm>
void sumDistributedWellEntries(Dune::DynamicMatrix<Scalar>& mat,
                               Dune::DynamicVector<Scalar>& vec,
                               const Comm& comm);


// explicit transpose of a dense matrix due to compilation problems
// used for calculating quasiimpes well weights
template <class DenseMatrix>
DenseMatrix transposeDenseDynMatrix(const DenseMatrix& M);

/// Helper to check whether the well is under zero production rate control
bool rateControlWithZeroProdTarget(const WellProductionControls& controls,
                                   WellProducerCMode mode);

/// Helper to check whether the well is under zero injection rate control
bool rateControlWithZeroInjTarget(const WellInjectionControls& controls,
                                  WellInjectorCMode mode);

} // namespace wellhelpers
} // namespace Opm

#endif
