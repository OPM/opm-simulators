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
 * \copydoc Opm::FluidSystems::SinglePhase
 */
#ifndef OPM_SINGLE_PHASE_FLUIDSYSTEM_HPP
#define OPM_SINGLE_PHASE_FLUIDSYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/TabulatedComponent.hpp>

#include <opm/material/common/Unused.hpp>

#include <limits>
#include <cassert>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for single phase models.
 *
 * The fluid is defined as a template parameter. For existing
 * components the Opm::LiquidPhase<Component> and
 * Opm::GasPhase<Component> may be used.
 */
template <class Scalar, class Fluid>
class SinglePhase
    : public BaseFluidSystem<Scalar, SinglePhase<Scalar, Fluid> >
{
    typedef SinglePhase<Scalar, Fluid> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    struct ParameterCache : public Opm::NullParameterCache<Evaluation>
    {};

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx OPM_OPTIM_UNUSED)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::name();
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::isLiquid();
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealGas(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isIdealGas();
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx OPM_OPTIM_UNUSED)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::name();
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    //! \copydoc BaseFluidSystem::init
    static void init()
    { }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::density(T, p);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::viscosity(T, p);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& /*fluidState*/,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

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
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::enthalpy(T, p);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::thermalConductivity(T, p);
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::heatCapacity(T, p);
    }
};

}} // namespace Opm, FluidSystems

#endif
