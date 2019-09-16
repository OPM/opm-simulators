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
 * \copydoc Opm::BlackOilLocalResidual
 */
#ifndef EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH
#define EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilenergymodules.hh"
#include "blackoilfoammodules.hh"

#include <opm/material/fluidstates/BlackOilFluidState.hpp>

namespace Opm {
/*!
 * \ingroup BlackOilModel
 *
 * \brief Calculates the local residual of the black oil model.
 */
template <class TypeTag>
class BlackOilLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static const bool waterEnabled = Indices::waterEnabled;
    static const bool gasEnabled = Indices::gasEnabled;
    static const bool oilEnabled = Indices::oilEnabled;
    static const bool compositionSwitchEnabled = (compositionSwitchIdx >= 0);

    static constexpr bool blackoilConserveSurfaceVolume = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume);
    static constexpr bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef BlackOilSolventModule<TypeTag> SolventModule;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;
    typedef BlackOilEnergyModule<TypeTag> EnergyModule;
    typedef BlackOilFoamModule<TypeTag> FoamModule;

public:
    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified point in time
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        storage = 0.0;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                if (Indices::numPhases == 3) { // add trivial equation for the pseudo phase
                    unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                    if (timeIdx == 0)
                        storage[conti0EqIdx + activeCompIdx] = Opm::variable<LhsEval>(0.0, conti0EqIdx + activeCompIdx);
                    else
                        storage[conti0EqIdx + activeCompIdx] = 0.0;
                }
                continue;
            }

            unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            LhsEval surfaceVolume =
                Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

            storage[conti0EqIdx + activeCompIdx] += surfaceVolume;

            // account for dissolved gas
            if (phaseIdx == oilPhaseIdx && FluidSystem::enableDissolvedGas()) {
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                storage[conti0EqIdx + activeGasCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs())
                    * surfaceVolume;
            }

            // account for vaporized oil
            if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedOil()) {
                unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                storage[conti0EqIdx + activeOilCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rv())
                    * surfaceVolume;
            }
        }

        adaptMassConservationQuantities_(storage, intQuants.pvtRegionIndex());

        // deal with solvents (if present)
        SolventModule::addStorage(storage, intQuants);

        // deal with polymer (if present)
        PolymerModule::addStorage(storage, intQuants);

        // deal with energy (if present)
        EnergyModule::addStorage(storage, intQuants);

        // deal with foam (if present)
        FoamModule::addStorage(storage, intQuants);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        assert(timeIdx == 0);

        flux = 0.0;

        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            unsigned pvtRegionIdx = up.pvtRegionIndex();
            if (upIdx == focusDofIdx)
                evalPhaseFluxes_<Evaluation>(flux, phaseIdx, pvtRegionIdx, extQuants, up.fluidState());
            else
                evalPhaseFluxes_<Scalar>(flux, phaseIdx, pvtRegionIdx, extQuants, up.fluidState());
        }

        // deal with solvents (if present)
        SolventModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with polymer (if present)
        PolymerModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with energy (if present)
        EnergyModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with foam (if present)
        FoamModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);

        // scale the source term of the energy equation
        if (enableEnergy)
            source[Indices::contiEnergyEqIdx] *= GET_PROP_VALUE(TypeTag, BlackOilEnergyScalingFactor);
    }

    /*!
     * \brief Helper function to calculate the flux of mass in terms of conservation
     *        quantities via specific fluid phase over a face.
     */
    template <class UpEval, class FluidState>
    static void evalPhaseFluxes_(RateVector& flux,
                                 unsigned phaseIdx,
                                 unsigned pvtRegionIdx,
                                 const ExtensiveQuantities& extQuants,
                                 const FluidState& upFs)
    {
        const auto& invB = Opm::getInvB_<FluidSystem, FluidState, UpEval>(upFs, phaseIdx, pvtRegionIdx);
        const auto& surfaceVolumeFlux = invB*extQuants.volumeFlux(phaseIdx);
        unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

        if (blackoilConserveSurfaceVolume)
            flux[conti0EqIdx + activeCompIdx] += surfaceVolumeFlux;
        else
            flux[conti0EqIdx + activeCompIdx] += surfaceVolumeFlux*FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);

        if (phaseIdx == oilPhaseIdx) {
            // dissolved gas (in the oil phase).
            if (FluidSystem::enableDissolvedGas()) {
                const auto& Rs = Opm::BlackOil::getRs_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeGasCompIdx] += Rs*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeGasCompIdx] += Rs*surfaceVolumeFlux*FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
        }
        else if (phaseIdx == gasPhaseIdx) {
            // vaporized oil (in the gas phase).
            if (FluidSystem::enableVaporizedOil()) {
                const auto& Rv = Opm::BlackOil::getRv_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeOilCompIdx] += Rv*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeOilCompIdx] += Rv*surfaceVolumeFlux*FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
        }
    }

    /*!
     * \brief Helper function to convert the mass-related parts of a Dune::FieldVector
     *        that stores conservation quantities in terms of "surface-volume" to the
     *        conservation quantities used by the model.
     *
     * Depending on the value of the BlackoilConserveSurfaceVolume property, the model
     * either conserves mass by means of "surface volume" of the components or mass
     * directly. In the former case, this method is a no-op; in the latter, the values
     * passed are multiplied by their respective pure component's density at surface
     * conditions.
     */
    template <class Scalar>
    static void adaptMassConservationQuantities_(Dune::FieldVector<Scalar, numEq>& container, unsigned pvtRegionIdx)
    {
        if (blackoilConserveSurfaceVolume)
            return;

        // convert "surface volume" to mass. this is complicated a bit by the fact that
        // not all phases are necessarily enabled. (we here assume that if a fluid phase
        // is disabled, its respective "main" component is not considered as well.)

        if (waterEnabled) {
            unsigned activeWaterCompIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
            container[conti0EqIdx + activeWaterCompIdx] *=
                FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
        }

        if (gasEnabled) {
            unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
            container[conti0EqIdx + activeGasCompIdx] *=
                FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
        }

        if (oilEnabled) {
            unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
            container[conti0EqIdx + activeOilCompIdx] *=
                FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
        }
    }
};

} // namespace Opm

#endif
