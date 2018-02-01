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
 * \copydoc Opm::FluidSystems::H2ON2LiquidPhase
 */
#ifndef OPM_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HPP
#define OPM_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include <opm/material/binarycoefficients/H2O_N2.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/material/common/Exceptions.hpp>

#include <iostream>
#include <cassert>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A liquid-phase-only fluid system with water and nitrogen as
 *        components.
 */
template <class Scalar>
class H2ON2LiquidPhase
    : public BaseFluidSystem<Scalar, H2ON2LiquidPhase<Scalar> >
{
    typedef H2ON2LiquidPhase<Scalar> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    // convenience typedefs
    typedef Opm::H2O<Scalar> IapwsH2O;
    typedef Opm::TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;
    typedef Opm::N2<Scalar> SimpleN2;

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

    //! Index of the liquid phase
    static const int liquidPhaseIdx = 0;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx OPM_OPTIM_UNUSED)
    {
        assert(phaseIdx == liquidPhaseIdx);

        return "liquid";
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned /*phaseIdx*/)
    {
        //assert(phaseIdx == liquidPhaseIdx);
        return true; //only water phase present
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return false; // not a gas (only liquid phase present)
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! The index of the water component
    static const int H2OIdx = 0;
    //! The index of the component for molecular nitrogen
    static const int N2Idx = 1;

    //! The type of the component for pure water
    typedef TabulatedH2O H2O;
    //typedef SimpleH2O H2O;
    //typedef IapwsH2O H2O;

    //! The type of the component for pure molecular nitrogen
    typedef SimpleN2 N2;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            H2O::name(),
            N2::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == N2Idx)
            ? N2::molarMass()
            : 1e30;
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == N2Idx)
            ? N2::criticalTemperature()
            : 1e30;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == N2Idx)
            ? N2::criticalPressure()
            : 1e30;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::acentricFactor()
            : (compIdx == N2Idx)
            ? N2::acentricFactor()
            : 1e30;
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/50,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/50);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges.
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        LhsEval sumMoleFrac = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoleFrac += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, compIdx));

        assert(phaseIdx == liquidPhaseIdx);

        // assume ideal mixture where each molecule occupies the same volume regardless
        // of whether it is water or nitrogen.
        const LhsEval& clH2O = H2O::liquidDensity(T, p)/H2O::molarMass();

        const auto& xlH2O = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, H2OIdx));
        const auto& xlN2 = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, N2Idx));

        return clH2O*(H2O::molarMass()*xlH2O + N2::molarMass()*xlN2)/sumMoleFrac;
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // assume pure water for the liquid phase
        return H2O::liquidViscosity(T, p);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        assert(0 <= compIdx && compIdx < numComponents);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (compIdx == H2OIdx)
            return H2O::vaporPressure(T)/p;
        return Opm::BinaryCoeff::H2O_N2::henry(T)/p;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned /*compIdx*/)

    {
        assert(phaseIdx == liquidPhaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        assert (phaseIdx == liquidPhaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        // TODO: way to deal with the solutes???
        return H2O::liquidEnthalpy(T, p);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       const unsigned phaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return H2O::liquidThermalConductivity(T, p);
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        assert (phaseIdx == liquidPhaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        return H2O::liquidHeatCapacity(T, p);
    }
};

} // namespace FluidSystems

} // namespace Opm

#endif
