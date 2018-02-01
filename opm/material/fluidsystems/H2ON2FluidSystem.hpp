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
 * \copydoc Opm::FluidSystems::H2ON2
 */
#ifndef OPM_H2O_N2_FLUID_SYSTEM_HPP
#define OPM_H2O_N2_FLUID_SYSTEM_HPP

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
 * \brief A two-phase fluid system with water and nitrogen as components.
 */
template <class Scalar>
class H2ON2
    : public BaseFluidSystem<Scalar, H2ON2<Scalar> >
{
    typedef H2ON2<Scalar> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    // convenience typedefs
    typedef Opm::IdealGas<Scalar> IdealGas;
    typedef Opm::H2O<Scalar> IapwsH2O;
    typedef Opm::TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;
    typedef Opm::N2<Scalar> SimpleN2;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = NullParameterCache<Evaluation>;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int liquidPhaseIdx = 0;
    //! Index of the gas phase
    static const int gasPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {
            "liquid",
            "gas"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        return
            (phaseIdx == gasPhaseIdx)
            ? true
            :H2O::liquidIsCompressible();// the water component decides for the liquid phase...
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return
            (phaseIdx == gasPhaseIdx)
            ? H2O::gasIsIdeal() && N2::gasIsIdeal() // let the components decide
            : false; // not a gas
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

    //! The component index of water
    static const int H2OIdx = 0;
    //! The component index of molecular nitrogen
    static const int N2Idx = 1;

    //! The component for pure water
    typedef TabulatedH2O H2O;
    //typedef SimpleH2O H2O;
    //typedef IapwsH2O H2O;

    //! The component for pure nitrogen
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
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        return (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == N2Idx)
            ? N2::criticalTemperature()
            : 1e30;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        return (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == N2Idx)
            ? N2::criticalPressure()
            : 1e30;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
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
     *        problem specific temperature and pressure ranges
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

    /*!
     * \copydoc BaseFluidSystem::density
     */
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

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            // assume ideal mixture where each molecule occupies the same volume regardless
            // of whether it is water or nitrogen.
            const LhsEval& clH2O = H2O::liquidDensity(T, p)/H2O::molarMass();

            const auto& xlH2O = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, H2OIdx));
            const auto& xlN2 = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, N2Idx));

            return clH2O*(H2O::molarMass()*xlH2O + N2::molarMass()*xlN2)/sumMoleFrac;
        }

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // assume ideal mixture: steam and nitrogen don't "distinguish" each other
        const auto& xgH2O = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
        const auto& xgN2 = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, N2Idx));
        const auto& rho_gH2O = H2O::gasDensity(T, p*xgH2O);
        const auto& rho_gN2 = N2::gasDensity(T, p*xgN2);
        return (rho_gH2O + rho_gN2)/Opm::max(1e-5, sumMoleFrac);
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

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         */
        LhsEval muResult = 0;
        const LhsEval mu[numComponents] = {
            H2O::gasViscosity(T, H2O::vaporPressure(T)),
            N2::gasViscosity(T, p)
        };

        LhsEval sumx = 0.0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, compIdx));
        sumx = Opm::max(1e-10, sumx);

        for (unsigned i = 0; i < numComponents; ++i) {
            LhsEval divisor = 0;
            for (unsigned j = 0; j < numComponents; ++j) {
                LhsEval phiIJ = 1 + Opm::sqrt(mu[i]/mu[j]) * std::pow(molarMass(j)/molarMass(i), 1/4.0);
                phiIJ *= phiIJ;
                phiIJ /= std::sqrt(8*(1 + molarMass(i)/molarMass(j)));
                divisor +=
                    Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, j))
                    /sumx*phiIJ;
            }
            muResult +=
                Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, i))
                /sumx*mu[i]/divisor;
        }
        return muResult;
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return Opm::BinaryCoeff::H2O_N2::henry(T)/p;
        }

        assert(phaseIdx == gasPhaseIdx);

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned /*compIdx*/)

    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);
        return BinaryCoeff::H2O_N2::gasDiffCoeff(T, p);
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            // TODO: correct way to deal with the solutes???
            return H2O::liquidEnthalpy(T, p);
        }

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // assume ideal mixture: Molecules of one component don't discriminate between
        // their own kind and molecules of the other component.
        const auto& XgH2O = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, H2OIdx));
        const auto& XgN2 = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, N2Idx));

        LhsEval hH2O = XgH2O*H2O::gasEnthalpy(T, p);
        LhsEval hN2 = XgN2*N2::gasEnthalpy(T, p);
        return hH2O + hN2;
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
        if (phaseIdx == liquidPhaseIdx) // liquid phase
            return H2O::liquidThermalConductivity(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // return the sum of the partial conductivity of Nitrogen and Steam
        const auto& xH2O = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, H2OIdx));
        const auto& xN2 = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, N2Idx));

        // Assuming Raoult's, Daltons law and ideal gas in order to obtain the
        // partial pressures in the gas phase
        const auto& lambdaN2 = N2::gasThermalConductivity(T, p*xN2);
        const auto& lambdaH2O = H2O::gasThermalConductivity(T, p*xH2O);

        return lambdaN2 + lambdaH2O;
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& xAlphaH2O = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, H2OIdx));
        const auto& xAlphaN2 = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, N2Idx));
        const auto& XAlphaH2O = Opm::decay<LhsEval>(fluidState.massFraction(phaseIdx, H2OIdx));
        const auto& XAlphaN2 = Opm::decay<LhsEval>(fluidState.massFraction(phaseIdx, N2Idx));

        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidHeatCapacity(T, p);

        assert(phaseIdx == gasPhaseIdx);

        // for the gas phase, assume ideal mixture
        LhsEval c_pN2;
        LhsEval c_pH2O;

        c_pN2 = N2::gasHeatCapacity(T, p*xAlphaN2);
        c_pH2O = H2O::gasHeatCapacity(T, p*xAlphaH2O);

        // mingle both components together. this assumes that there is no "cross
        // interaction" between both flavors of molecules.
        return XAlphaH2O*c_pH2O + XAlphaN2*c_pN2;
    }
};

} // namespace FluidSystems

} // namespace Opm

#endif
