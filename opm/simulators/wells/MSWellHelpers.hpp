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

#include <cmath>
#include <cassert>
#include <cstddef>

namespace Dune {
template<class Matrix> class UMFPack;
}

namespace Opm {

template<class Scalar> class ParallelWellInfo;
class DeferredLogger;
class SICD;

namespace mswellhelpers
{

    /// \brief A wrapper around the B matrix for distributed MS wells
    ///
    /// For wells the B matrix, is basically a multiplication
    /// of the equation of the perforated cells followed by a reduction
    /// (summation) of these to the well equations.
    ///
    /// This class does that in the functions mv and mmv (from the DUNE
    /// matrix interface.
    ///
    /// \tparam MatrixType The MatrixType of the Matrix B. From this, we
    ///                    deduce the Scalar used for the computation.
    template<class MatrixType>
    class ParallellMSWellB
    {
    public:
        using Scalar = typename MatrixType::field_type;
        ParallellMSWellB(const MatrixType& B,
                         const ParallelWellInfo<Scalar>& parallel_well_info);

        //! y = A x
        template<class X, class Y>
        void mv (const X& x, Y& y) const;

        //! y = A x
        template<class X, class Y>
        void mmv (const X& x, Y& y) const;

    private:
        const MatrixType& B_;
        const ParallelWellInfo<Scalar>& parallel_well_info_;
    };

    /// Applies umfpack and checks for singularity
    template <typename MatrixType, typename VectorType>
    VectorType
    applyUMFPack(Dune::UMFPack<MatrixType>& linsolver,
                 VectorType x);

    /// Coarse CPRW diagonal for a multisegment well: lambda' D(:,p), summed
    /// over every segment row and every column of D.
    ///
    /// This is the Galerkin diagonal for a prolongation that spreads the well's
    /// single coarse unknown over all of its segment pressures, which is the
    /// prolongation the coarse column already assumes. Unlike the row-sum
    /// convention it keeps the segment-to-segment coupling D carries, which is
    /// most of a well once there is one segment per connection.
    ///
    /// \param D                  the well's diagonal block, segment rows by segment columns
    /// \param lambda             well weights, one per conservation equation
    /// \param pressureColumn     the segment-pressure column of a D block
    ///
    /// Returns 1 if the contraction cancels exactly, since a zero would make
    /// the coarse pressure system singular.
    ///
    /// SystemCprwPressureStage::assembleCoarseMatrix runs the same inner
    /// kernel over its merged matrix, but distributes into a whole coarse row
    /// with per-block weights instead of accumulating one well's D into a
    /// single scalar, so the two are not one function.
    template <class DiagMatWell, class WellWeight>
    typename DiagMatWell::field_type
    contractCprWellDiagonal(const DiagMatWell& D,
                            const WellWeight& lambda,
                            const int pressureColumn)
    {
        using Scalar = typename DiagMatWell::field_type;
        using Block = typename DiagMatWell::block_type;

        // lambda indexes the conservation-equation rows of a D block, so the
        // weights must not outrun them.  Catches a future reordering of the
        // well primary variables that puts something else in the first rows.
        assert(lambda.size() <= static_cast<std::size_t>(Block::rows));

        Scalar diag = 0.0;
        for (std::size_t row = 0; row < D.N(); ++row) {
            for (auto col = D[row].begin(), end = D[row].end(); col != end; ++col) {
                for (std::size_t i = 0; i < lambda.size(); ++i) {
                    diag += lambda[i] * (*col)[i][pressureColumn];
                }
            }
        }
        return (std::abs(diag) > 0.0) ? diag : Scalar{1.0};
    }



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
    template <typename ValueType, typename Scalar>
    ValueType frictionPressureLoss(const Scalar l, const Scalar diameter,
                                   const Scalar area, const Scalar roughness,
                                   const ValueType& density,
                                   const ValueType& w, const ValueType& mu);


    template <typename ValueType, typename Scalar>
    ValueType valveContrictionPressureLoss(const ValueType& mass_rate,
                                           const ValueType& density,
                                           const Scalar area_con, const Scalar cv);


    template <typename ValueType, typename Scalar>
    ValueType velocityHead(const Scalar area, const ValueType& mass_rate,
                           const ValueType& density);


    // calculating the viscosity of oil-water emulsion at local conditons
    template <typename ValueType, typename Scalar>
    ValueType emulsionViscosity(const ValueType& water_fraction,
                                const ValueType& water_viscosity,
                                const ValueType& oil_fraction,
                                const ValueType& oil_viscosity,
                                const SICD& sicd);

} // namespace mswellhelpers
} // namespace Opm

#endif
