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
 * \copydoc Opm::FlashBoundaryRateVector
 */
#ifndef EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH

#include "flashproperties.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/material/common/Valgrind.hpp>

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief Implements a boundary vector for the fully implicit
 *        compositional multi-phase model which is based on flash
 *        calculations.
 */
template <class TypeTag>
class FlashBoundaryRateVector : public GetPropType<TypeTag, Properties::RateVector>
{
    using ParentType = GetPropType<TypeTag, Properties::RateVector>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;
    using Toolbox = Opm::MathToolbox<Evaluation>;

public:
    FlashBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc
     * ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    FlashBoundaryRateVector(const Evaluation& value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const
     * ImmiscibleBoundaryRateVector& )
     */
    FlashBoundaryRateVector(const FlashBoundaryRateVector& value) = default;
    FlashBoundaryRateVector& operator=(const FlashBoundaryRateVector& value) = default;

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
        unsigned focusDofIdx = context.focusDofIndex();
        unsigned interiorDofIdx = context.interiorScvIndex(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = Evaluation(0.0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Evaluation density;
            if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
                if (focusDofIdx == interiorDofIdx)
                    density = fluidState.density(phaseIdx);
                else
                    density = Opm::getValue(fluidState.density(phaseIdx));
            }
            else if (focusDofIdx == interiorDofIdx)
                density = insideIntQuants.fluidState().density(phaseIdx);
            else
                density = Opm::getValue(insideIntQuants.fluidState().density(phaseIdx));

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Evaluation molarity;
                if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
                    if (focusDofIdx == interiorDofIdx)
                        molarity = fluidState.molarity(phaseIdx, compIdx);
                    else
                        molarity = Opm::getValue(fluidState.molarity(phaseIdx, compIdx));
                }
                else if (focusDofIdx == interiorDofIdx)
                    molarity = insideIntQuants.fluidState().molarity(phaseIdx, compIdx);
                else
                    molarity = Opm::getValue(insideIntQuants.fluidState().molarity(phaseIdx, compIdx));

                // add advective flux of current component in current
                // phase
                (*this)[conti0EqIdx + compIdx] += extQuants.volumeFlux(phaseIdx)*molarity;
            }

            if (enableEnergy) {
                Evaluation specificEnthalpy;
                if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx)) {
                    if (focusDofIdx == interiorDofIdx)
                        specificEnthalpy = fluidState.enthalpy(phaseIdx);
                    else
                        specificEnthalpy = Opm::getValue(fluidState.enthalpy(phaseIdx));
                }
                else if (focusDofIdx == interiorDofIdx)
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);
                else
                    specificEnthalpy = Opm::getValue(insideIntQuants.fluidState().enthalpy(phaseIdx));

                Evaluation enthalpyRate = density*extQuants.volumeFlux(phaseIdx)*specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        // thermal conduction
        EnergyModule::addToEnthalpyRate(*this, EnergyModule::thermalConductionRate(extQuants));

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i) {
            Opm::Valgrind::CheckDefined((*this)[i]);
        }
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

        // we only allow fluxes in the direction opposite to the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::min(0.0, val);
        }
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

        // we only allow fluxes in the same direction as the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::max(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = Evaluation(0.0); }
};

} // namespace Opm

#endif
