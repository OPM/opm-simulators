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
 * \copydoc Opm::BlackOilBoundaryRateVector
 */
#ifndef EWOMS_BLACK_OIL_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_BLACK_OIL_BOUNDARY_RATE_VECTOR_HH

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include <opm/models/blackoil/blackoilenergymodules.hh>

#include <algorithm>
#include <type_traits>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief Implements a boundary vector for the fully implicit black-oil model.
 */
template <class TypeTag>
class BlackOilBoundaryRateVector : public GetPropType<TypeTag, Properties::RateVector>
{
    using ParentType = GetPropType<TypeTag, Properties::RateVector>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { contiEnergyEqIdx = Indices::contiEnergyEqIdx };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableMICP = Indices::enableMICP };

    static constexpr bool blackoilConserveSurfaceVolume =
        getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>();

    using EnergyModule = BlackOilEnergyModule<TypeTag, enableEnergy>;

public:
    /*!
     * \brief Default constructor
     */
    BlackOilBoundaryRateVector() = default;

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    BlackOilBoundaryRateVector(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector& )
     */
    BlackOilBoundaryRateVector(const BlackOilBoundaryRateVector& value) = default;
    BlackOilBoundaryRateVector& operator=(const BlackOilBoundaryRateVector& value) = default;

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context,
                     unsigned bfIdx,
                     unsigned timeIdx,
                     const FluidState& fluidState)
    {
        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);
        const unsigned focusDofIdx = context.focusDofIndex();
        const unsigned interiorDofIdx = context.interiorScvIndex(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const auto& pBoundary = fluidState.pressure(phaseIdx);
            const Evaluation& pInside = insideIntQuants.fluidState().pressure(phaseIdx);

            RateVector tmp;

            // mass conservation
            if (pBoundary < pInside) {
                // outflux
                LocalResidual::template evalPhaseFluxes_<Evaluation>(tmp,
                                                                     phaseIdx,
                                                                     insideIntQuants.pvtRegionIndex(),
                                                                     extQuants,
                                                                     insideIntQuants.fluidState());
            }
            else if (pBoundary > pInside) {
                using RhsEval = std::conditional_t<std::is_same_v<typename FluidState::Scalar, Evaluation>,
                                                   Evaluation, Scalar>;
                // influx
                LocalResidual::template evalPhaseFluxes_<RhsEval>(tmp,
                                                                  phaseIdx,
                                                                  insideIntQuants.pvtRegionIndex(),
                                                                  extQuants,
                                                                  fluidState);
            }

            for (unsigned i = 0; i < tmp.size(); ++i) {
                (*this)[i] += tmp[i];
            }

            // energy conservation
            if constexpr (enableEnergy) {
                Evaluation density;
                Evaluation specificEnthalpy;
                if (pBoundary > pInside) {
                    if (focusDofIdx == interiorDofIdx) {
                        density = fluidState.density(phaseIdx);
                        specificEnthalpy = fluidState.enthalpy(phaseIdx);
                    }
                    else {
                        density = getValue(fluidState.density(phaseIdx));
                        specificEnthalpy = getValue(fluidState.enthalpy(phaseIdx));
                    }
                }
                else if (focusDofIdx == interiorDofIdx) {
                    density = insideIntQuants.fluidState().density(phaseIdx);
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);
                }
                else {
                    density = getValue(insideIntQuants.fluidState().density(phaseIdx));
                    specificEnthalpy = getValue(insideIntQuants.fluidState().enthalpy(phaseIdx));
                }

                const Evaluation enthalpyRate = density * extQuants.volumeFlux(phaseIdx) * specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate *
                                                       getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>());
            }
        }

        if constexpr (enableSolvent) {
            (*this)[Indices::contiSolventEqIdx] = extQuants.solventVolumeFlux();
            if (blackoilConserveSurfaceVolume) {
                (*this)[Indices::contiSolventEqIdx] *= insideIntQuants.solventInverseFormationVolumeFactor();
            }
            else {
                (*this)[Indices::contiSolventEqIdx] *= insideIntQuants.solventDensity();
            }
        }

        if constexpr (enablePolymer) {
            (*this)[Indices::contiPolymerEqIdx] = extQuants.volumeFlux(FluidSystem::waterPhaseIdx) *
                                                  insideIntQuants.polymerConcentration();
        }

        if constexpr (enableMICP) {
            (*this)[Indices::contiMicrobialEqIdx] = extQuants.volumeFlux(FluidSystem::waterPhaseIdx) *
                                                    insideIntQuants.microbialConcentration();
            (*this)[Indices::contiOxygenEqIdx] = extQuants.volumeFlux(FluidSystem::waterPhaseIdx) *
                                                 insideIntQuants.oxygenConcentration();
            (*this)[Indices::contiUreaEqIdx] = extQuants.volumeFlux(FluidSystem::waterPhaseIdx) *
                                               insideIntQuants.ureaConcentration();
            // since the urea concentration can be much larger than 1, then we apply a scaling factor
            (*this)[Indices::contiUreaEqIdx] *= getPropValue<TypeTag, Properties::BlackOilUreaScalingFactor>();
        }

        // make sure that the right mass conservation quantities are used
        LocalResidual::adaptMassConservationQuantities_(*this, insideIntQuants.pvtRegionIndex());

        // heat conduction
        if constexpr (enableEnergy) {
            EnergyModule::addToEnthalpyRate(*this, extQuants.energyFlux() *
                                                   getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>());
        }

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i) {
            Valgrind::CheckDefined((*this)[i]);
        }
        Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setInFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer
        // unit normal
        std::for_each(this->begin(), this->end(),
                      [](auto& val) { val = std::min(Scalar(0), val); });
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setOutFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context& context,
                    unsigned bfIdx,
                    unsigned timeIdx,
                    const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer
        // unit normal
        std::for_each(this->begin(), this->end(),
                      [](auto& val) { val = std::max(Scalar(0), val); });
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = Scalar(0); }

    /*!
     * \copydoc Specify an energy flux that corresponds to the thermal conduction from
     *          the domain boundary
     *
     * This means that a "thermal flow" boundary is a no-flow condition for mass and thermal
     * conduction for energy.
     */
    template <class Context, class FluidState>
    void setThermalFlow([[maybe_unused]] const Context& context,
                        [[maybe_unused]] unsigned bfIdx,
                        [[maybe_unused]] unsigned timeIdx,
                        [[maybe_unused]] const FluidState& boundaryFluidState)
    {
        // set the mass no-flow condition
        setNoFlow();

        // if we do not conserve energy there is nothing we should do in addition
        if constexpr (enableEnergy) {
            ExtensiveQuantities extQuants;
            extQuants.updateBoundary(context, bfIdx, timeIdx, boundaryFluidState);

            (*this)[contiEnergyEqIdx] += extQuants.energyFlux();

#ifndef NDEBUG
            for (unsigned i = 0; i < numEq; ++i) {
                Valgrind::CheckDefined((*this)[i]);
            }
            Valgrind::CheckDefined(*this);
#endif
        }
    }
};

} // namespace Opm

#endif
