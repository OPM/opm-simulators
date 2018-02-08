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
 * \copydoc Opm::FluidSystems::H2OAir
 */
#ifndef OPM_H2O_AIR_SYSTEM_HPP
#define OPM_H2O_AIR_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/binarycoefficients/H2O_Air.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/TabulatedComponent.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <iostream>
#include <cassert>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system with a liquid and a gaseous phase and water and air
 *        as components.
 *
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 */
template <class Scalar,
          //class H2Otype = Opm::SimpleH2O<Scalar>,
          class H2Otype = Opm::TabulatedComponent<Scalar, Opm::H2O<Scalar> >>
class H2OAir
    : public BaseFluidSystem<Scalar, H2OAir<Scalar, H2Otype> >
{
    typedef H2OAir<Scalar,H2Otype> ThisType;
    typedef BaseFluidSystem <Scalar, ThisType> Base;
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    template <class Evaluation>
    struct ParameterCache : public Opm::NullParameterCache<Evaluation>
    {};

    //! The type of the water component used for this fluid system
    typedef H2Otype H2O;
    //! The type of the air component used for this fluid system
    typedef Opm::Air<Scalar> Air;

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! The index of the liquid phase
    static const int liquidPhaseIdx = 0;
    //! The index of the gas phase
    static const int gasPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        switch (phaseIdx) {
        case liquidPhaseIdx: return "liquid";
        case gasPhaseIdx: return "gas";
        };
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
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
        return (phaseIdx == gasPhaseIdx)
            // ideal gases are always compressible
            ? true
            :
            // the water component decides for the liquid phase...
            H2O::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    {
        return
            (phaseIdx == gasPhaseIdx)
            ? H2O::gasIsIdeal() && Air::gasIsIdeal()
            : false;
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
    //! The index of the air component
    static const int AirIdx = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::name();
        case AirIdx: return Air::name();
        };
        throw std::logic_error("Invalid component index "+std::to_string(compIdx));
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == AirIdx)
            ? Air::molarMass()
            : 1e30;
    }


    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == AirIdx)
            ? Air::criticalTemperature()
            : 1e30;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == AirIdx)
            ? Air::criticalPressure()
            : 1e30;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::acentricFactor()
            : (compIdx == AirIdx)
            ? Air::acentricFactor()
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
        if (H2O::isTabulated)
            init(/*tempMin=*/273.15,
                 /*tempMax=*/623.15,
                 /*numTemp=*/50,
                 /*pMin=*/-10,
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
            H2O::init(tempMin, tempMax, nTemp,
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
        LhsEval p;
        if (isCompressible(phaseIdx))
            p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        else {
            // random value which will hopefully cause things to blow
            // up if it is used in a calculation!
            p = - 1e30;
            Valgrind::SetUndefined(p);
        }


        LhsEval sumMoleFrac = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoleFrac += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, compIdx));

        if (phaseIdx == liquidPhaseIdx)
        {
            // assume ideal mixture: Molecules of one component don't discriminate
            // between their own kind and molecules of the other component.
            const LhsEval& clH2O = H2O::liquidDensity(T, p)/H2O::molarMass();

            const auto& xlH2O = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, H2OIdx));
            const auto& xlAir = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, AirIdx));

            return clH2O*(H2O::molarMass()*xlH2O + Air::molarMass()*xlAir)/sumMoleFrac;
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            LhsEval partialPressureH2O =
                Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx))
                *Opm::decay<LhsEval>(fluidState.pressure(gasPhaseIdx));

            LhsEval partialPressureAir =
                Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, AirIdx))
                *Opm::decay<LhsEval>(fluidState.pressure(gasPhaseIdx));

            return H2O::gasDensity(T, partialPressureH2O) + Air::gasDensity(T, partialPressureAir);
        }
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
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

        if (phaseIdx == liquidPhaseIdx)
        {
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            // couldn't find a way to solve the mixture problem
            return H2O::liquidViscosity(T, p);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            /* Wilke method. See:
             *
             * See: R. Reid, et al.: The Properties of Gases and Liquids,
             * 4th edition, McGraw-Hill, 1987, 407-410 or
             * 5th edition, McGraw-Hill, 2000, p. 9.21/22
             */
            LhsEval muResult = 0;
            const LhsEval mu[numComponents] = {
                H2O::gasViscosity(T, H2O::vaporPressure(T)),
                Air::gasViscosity(T, p)
            };

            for (unsigned i = 0; i < numComponents; ++i) {
                LhsEval divisor = 0;
                for (unsigned j = 0; j < numComponents; ++j) {
                    LhsEval phiIJ =
                        1 +
                        Opm::sqrt(mu[i]/mu[j]) * // 1 + (mu[i]/mu[j]^1/2
                        std::pow(molarMass(j)/molarMass(i), 1./4.0);   // (M[i]/M[j])^1/4

                    phiIJ *= phiIJ;
                    phiIJ /= std::sqrt(8*(1 + molarMass(i)/molarMass(j)));
                    divisor += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, j))*phiIJ;
                }
                const auto& xAlphaI = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, i));
                muResult += xAlphaI*mu[i]/divisor;
            }
            return muResult;
        }
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
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

        if (phaseIdx == liquidPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return Opm::BinaryCoeff::H2O_Air::henry(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval binaryDiffusionCoefficient(const FluidState& fluidState,
                                              const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                              unsigned phaseIdx,
                                              unsigned /*compIdx*/)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == liquidPhaseIdx)
            return BinaryCoeff::H2O_Air::liquidDiffCoeff(T, p);

        assert(phaseIdx == gasPhaseIdx);
        return BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);
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

        if (phaseIdx == liquidPhaseIdx)
        {
            // TODO: correct way to deal with the solutes???
            return H2O::liquidEnthalpy(T, p);
        }

        else if (phaseIdx == gasPhaseIdx)
        {
            LhsEval result = 0.0;
            result +=
                H2O::gasEnthalpy(T, p) *
                Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, H2OIdx));

            result +=
                Air::gasEnthalpy(T, p) *
                Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, AirIdx));
            return result;
        }
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const LhsEval& temperature =
            Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure =
            Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidThermalConductivity(temperature, pressure);
        else { // gas phase
            const LhsEval& lambdaDryAir = Air::gasThermalConductivity(temperature, pressure);

            const LhsEval& xAir =
                Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, AirIdx));
            const LhsEval& xH2O =
                Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, H2OIdx));
            LhsEval lambdaAir = xAir*lambdaDryAir;

            // Assuming Raoult's, Daltons law and ideal gas
            // in order to obtain the partial density of water in the air phase
            LhsEval partialPressure  = pressure*xH2O;

            LhsEval lambdaH2O =
                xH2O*H2O::gasThermalConductivity(temperature, partialPressure);
            return lambdaAir + lambdaH2O;
        }
    }
};

} // namespace FluidSystems

} // namespace Opm

#endif
