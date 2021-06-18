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
 * \copydoc Opm::TwoPhaseImmiscibleFluidSystem
 */
#ifndef OPM_TWO_PHASE_IMMISCIBLE_FLUID_SYSTEM_HPP
#define OPM_TWO_PHASE_IMMISCIBLE_FLUID_SYSTEM_HPP

#include <limits>
#include <cassert>

#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

namespace Opm {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The wetting and the non-wetting phase can be defined individually
 * via <tt>Opm::LiquidPhase<Component></tt> and
 * <tt>Opm::GasPhase<Component></tt>. These phases consist of one pure
 * component. With the help of this adapter class, the phase
 * properties can be accessed. This is suitable for pure two-phase
 * systems without compositional effects.
 */
template <class Scalar, class WettingPhase, class NonwettingPhase>
class TwoPhaseImmiscibleFluidSystem
    : public BaseFluidSystem<Scalar, TwoPhaseImmiscibleFluidSystem<Scalar, WettingPhase, NonwettingPhase> >
{
    // do not try to instanciate this class, it has only static members!
    TwoPhaseImmiscibleFluidSystem()
    {}

    typedef TwoPhaseImmiscibleFluidSystem<Scalar, WettingPhase, NonwettingPhase> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    template <class Evaluation>
    struct ParameterCache : public NullParameterCache<Evaluation>
    {};

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the wetting phase
    static const int wettingPhaseIdx = 0;
    //! Index of the non-wetting phase
    static const int nonWettingPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        static const char* name[] = {
            "wetting",
            "nonwetting"
        };
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return
            (phaseIdx == wettingPhaseIdx)
            ? WettingPhase::isLiquid()
            : NonwettingPhase::isLiquid();
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return
            (phaseIdx == wettingPhaseIdx)
            ? WettingPhase::isCompressible()
            : NonwettingPhase::isCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        return
            (phaseIdx == wettingPhaseIdx)
            ? WettingPhase::isIdealGas()
            : NonwettingPhase::isIdealGas();
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! Index of the wetting phase's component
    static const int wettingCompIdx = 0;
    //! Index of the non-wetting phase's component
    static const int nonWettingCompIdx = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        assert(compIdx < numComponents);

        if (compIdx == wettingCompIdx)
            return WettingPhase::name();
        return NonwettingPhase::name();
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        // let the fluids decide
        return
            (compIdx == wettingCompIdx)
            ? WettingPhase::molarMass()
            : NonwettingPhase::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wettingCompIdx)
            ? WettingPhase::criticalTemperature()
            : NonwettingPhase::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wettingCompIdx)
            ? WettingPhase::criticalPressure()
            : NonwettingPhase::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        // let the fluids decide
        return
            (compIdx == wettingCompIdx)
            ? WettingPhase::acentricFactor()
            : NonwettingPhase::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    //! \copydoc BaseFluidSystem::init
    static void init()
    {
        // two gaseous phases at once do not make sense physically!
        // (But two liquids are fine)
        assert(WettingPhase::isLiquid() || NonwettingPhase::isLiquid());
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        const auto& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == wettingPhaseIdx)
            return WettingPhase::density(temperature, pressure);
        return NonwettingPhase::density(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        const auto& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == wettingPhaseIdx)
            return WettingPhase::viscosity(temperature, pressure);
        return NonwettingPhase::viscosity(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& /*fluidState*/,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(phaseIdx < numPhases);
        assert(compIdx < numComponents);

        if (phaseIdx == compIdx)
            // TODO (?): calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        const auto& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == wettingPhaseIdx)
            return WettingPhase::enthalpy(temperature, pressure);
        return NonwettingPhase::enthalpy(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        const auto& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == wettingPhaseIdx)
            return WettingPhase::thermalConductivity(temperature, pressure);
        return NonwettingPhase::thermalConductivity(temperature, pressure);
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        assert(phaseIdx < numPhases);

        const auto& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == wettingPhaseIdx)
            return WettingPhase::heatCapacity(temperature, pressure);
        return NonwettingPhase::heatCapacity(temperature, pressure);
    }
};

} // namespace Opm

#endif
