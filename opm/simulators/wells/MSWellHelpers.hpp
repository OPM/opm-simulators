/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2020 Equinor ASA.

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


#ifndef OPM_MSWELLHELPERS_HEADER_INCLUDED
#define OPM_MSWELLHELPERS_HEADER_INCLUDED

#include <dune/istl/matrix.hh>

namespace Dune {
template<class Matrix> class UMFPack;
}

#include <memory>

namespace Opm {

class DeferredLogger;
class SICD;

namespace mswellhelpers
{

    /// Applies umfpack and checks for singularity
    template <typename MatrixType, typename VectorType>
    VectorType
    applyUMFPack(Dune::UMFPack<MatrixType>& linsolver,
                 VectorType x);



    /// Applies umfpack and checks for singularity
    template <typename VectorType, typename MatrixType>
    Dune::Matrix<typename MatrixType::block_type>
    invertWithUMFPack(const int size,
                      const int bsize,
                      Dune::UMFPack<MatrixType>& linsolver);



    // obtain y = D^-1 * x with a BICSSTAB iterative solver
    template <typename MatrixType, typename VectorType>
    VectorType
    invDX(const MatrixType& D, VectorType x, DeferredLogger& deferred_logger);

    // calculating the friction pressure loss
    // l is the segment length
    // area is the segment cross area
    // diameter is the segment inner diameter
    // w is mass flow rate through the segment
    // density is density
    // roughness is the absolute roughness
    // mu is the average phase viscosity
    template <typename ValueType>
    ValueType frictionPressureLoss(const double l, const double diameter,
                                   const double area, const double roughness,
                                   const ValueType& density,
                                   const ValueType& w, const ValueType& mu);


    template <typename ValueType>
    ValueType valveContrictionPressureLoss(const ValueType& mass_rate,
                                           const ValueType& density,
                                           const double area_con, const double cv);


    template <typename ValueType>
    ValueType velocityHead(const double area, const ValueType& mass_rate,
                           const ValueType& density);


    // calculating the viscosity of oil-water emulsion at local conditons
    template <typename ValueType>
    ValueType emulsionViscosity(const ValueType& water_fraction,
                                const ValueType& water_viscosity,
                                const ValueType& oil_fraction,
                                const ValueType& oil_viscosity,
                                const SICD& sicd);

} // namespace mswellhelpers
} // namespace Opm

#endif
