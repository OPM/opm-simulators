// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::RichardsLocalResidual
 */
#ifndef EWOMS_RICHARDS_LOCAL_RESIDUAL_HH
#define EWOMS_RICHARDS_LOCAL_RESIDUAL_HH

#include "richardsintensivequantities.hh"

#include "richardsextensivequantities.hh"

namespace Opm {

/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the residual for the Richards model.
 */
template <class TypeTag>
class RichardsLocalResidual : public GetPropType<TypeTag, Properties::DiscLocalResidual>
{
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidPhaseIdx = getPropValue<TypeTag, Properties::LiquidPhaseIndex>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Toolbox = Opm::MathToolbox<Evaluation>;

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);

        // partial time derivative of the wetting phase mass
        storage[contiEqIdx] =
            Toolbox::template decay<LhsEval>(intQuants.fluidState().density(liquidPhaseIdx))
            *Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(liquidPhaseIdx))
            *Toolbox::template decay<LhsEval>(intQuants.porosity());
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned focusDofIdx = elemCtx.focusDofIndex();
        unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(liquidPhaseIdx));

        const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        // compute advective mass flux of the liquid phase. This is slightly hacky
        // because it is specific to the element-centered finite volume method.
        const Evaluation& rho = up.fluidState().density(liquidPhaseIdx);
        if (focusDofIdx == upIdx)
            flux[contiEqIdx] = extQuants.volumeFlux(liquidPhaseIdx)*rho;
        else
            flux[contiEqIdx] = extQuants.volumeFlux(liquidPhaseIdx)*Toolbox::value(rho);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeSource
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    { elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx); }
};

} // namespace Opm

#endif
