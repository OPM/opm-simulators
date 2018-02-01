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
 * \copydoc Opm::FluidSystems::BrineCO2
 */
#ifndef OPM_BRINE_CO2_SYSTEM_HPP
#define OPM_BRINE_CO2_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/IdealGas.hpp>

#include <opm/material/components/Brine.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/TabulatedComponent.hpp>

#include <opm/material/binarycoefficients/H2O_CO2.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/binarycoefficients/H2O_N2.hpp>

#include <opm/material/common/Unused.hpp>

#include <iostream>

namespace Opm {
namespace FluidSystems {

/*!
 * \brief A two-phase fluid system with water and CO2.
 *
 * This fluid system uses the a tabulated CO2 component to achieve
 * high thermodynamic accuracy and thus requires the tables of the
 * sampling to be supplied as template argument.
 */
template <class Scalar, class CO2Tables>
class BrineCO2
    : public BaseFluidSystem<Scalar, BrineCO2<Scalar, CO2Tables> >
{
    typedef Opm::H2O<Scalar> H2O_IAPWS;
    typedef Opm::Brine<Scalar, H2O_IAPWS> Brine_IAPWS;
    typedef Opm::TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;
    typedef Opm::TabulatedComponent<Scalar, Brine_IAPWS> Brine_Tabulated;

    typedef H2O_Tabulated H2O;

public:
    template <class Evaluation>
    struct ParameterCache : public Opm::NullParameterCache<Evaluation>
    {};

    //! The binary coefficients for brine and CO2 used by this fluid system
    typedef Opm::BinaryCoeff::Brine_CO2<Scalar, CO2Tables> BinaryCoeffBrineCO2;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! The number of phases considered by the fluid system
    static const int numPhases = 2;

    //! The index of the liquid phase
    static const int liquidPhaseIdx = 0;
    //! The index of the gas phase
    static const int gasPhaseIdx = 1;

    /*!
     * \copydoc BaseFluidSystem::phaseName
     */
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {
            "liquid",
            "gas"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \copydoc BaseFluidSystem::isLiquid
     */
    static bool isLiquid(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != gasPhaseIdx;
    }

    /*!
     * \copydoc BaseFluidSystem::isIdealGas
     */
    static bool isIdealGas(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == gasPhaseIdx)
            return CO2::gasIsIdeal();
        return false;
    }

    /*!
     * \copydoc BaseFluidSystem::isIdealMixture
     */
    static bool isIdealMixture(unsigned phaseIdx OPM_OPTIM_UNUSED)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
     * \copydoc BaseFluidSystem::isCompressible
     */
    static bool isCompressible(unsigned phaseIdx OPM_OPTIM_UNUSED)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! The index of the brine component
    static const int BrineIdx = 0;
    //! The index of the CO2 component
    static const int CO2Idx = 1;

    //! The type of the component for brine used by the fluid system
    typedef Brine_Tabulated Brine;
    //! The type of the component for pure CO2 used by the fluid system
    typedef Opm::CO2<Scalar, CO2Tables> CO2;

    /*!
     * \copydoc BaseFluidSystem::componentName
     */
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            Brine::name(),
            CO2::name(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \copydoc BaseFluidSystem::molarMass
     */
    static Scalar molarMass(unsigned compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx==BrineIdx)
            ? Brine::molarMass()
            : CO2::molarMass();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::init
     */
    static void init()
    {
        init(/*startTemp=*/273.15, /*endTemp=*/623.15, /*tempSteps=*/50,
             /*startPressure=*/1e4, /*endPressure=*/40e6, /*pressureSteps=*/50);
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
            H2O_Tabulated::init(tempMin, tempMax, nTemp,
                                pressMin, pressMax, nPress);
        }

        // set the salinity of brine to the one used by the CO2 tables
        Brine_IAPWS::salinity = CO2Tables::brineSalinity;

        if (Brine::isTabulated) {
            Brine_Tabulated::init(tempMin, tempMax, nTemp,
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

        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == liquidPhaseIdx) {
            // use normalized composition for to calculate the density
            // (the relations don't seem to take non-normalized
            // compositions too well...)
            LhsEval xlBrine = Opm::min(1.0, Opm::max(0.0, Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, BrineIdx))));
            LhsEval xlCO2 = Opm::min(1.0, Opm::max(0.0,  Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, CO2Idx))));
            LhsEval sumx = xlBrine + xlCO2;
            xlBrine /= sumx;
            xlCO2 /= sumx;

            LhsEval result = liquidDensity_(temperature,
                                            pressure,
                                            xlBrine,
                                            xlCO2);

            Valgrind::CheckDefined(result);
            return result;
        }

        assert(phaseIdx == gasPhaseIdx);

        // use normalized composition for to calculate the density
        // (the relations don't seem to take non-normalized
        // compositions too well...)
        LhsEval xgBrine = Opm::min(1.0, Opm::max(0.0, Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, BrineIdx))));
        LhsEval xgCO2 = Opm::min(1.0, Opm::max(0.0,  Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, CO2Idx))));
        LhsEval sumx = xgBrine + xgCO2;
        xgBrine /= sumx;
        xgCO2 /= sumx;

        LhsEval result = gasDensity_(temperature,
                                     pressure,
                                     xgBrine,
                                     xgCO2);
        Valgrind::CheckDefined(result);
        return result;
    }

    /*!
     * \copydoc BaseFluidSystem::viscosity
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == liquidPhaseIdx) {
            // assume pure brine for the liquid phase. TODO: viscosity
            // of mixture
            LhsEval result = Brine::liquidViscosity(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }

        assert(phaseIdx == gasPhaseIdx);
        LhsEval result = CO2::gasViscosity(temperature, pressure);
        Valgrind::CheckDefined(result);
        return result;
    }

    /*!
     * \copydoc BaseFluidSystem::fugacityCoefficient
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == gasPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        assert(temperature > 0);
        assert(pressure > 0);

        // calulate the equilibrium composition for the given
        // temperature and pressure. TODO: calculateMoleFractions()
        // could use some cleanup.
        LhsEval xlH2O, xgH2O;
        LhsEval xlCO2, xgCO2;
        BinaryCoeffBrineCO2::calculateMoleFractions(temperature,
                                                    pressure,
                                                    Brine_IAPWS::salinity,
                                                    /*knownPhaseIdx=*/-1,
                                                    xlCO2,
                                                    xgH2O);

        // normalize the phase compositions
        xlCO2 = Opm::max(0.0, Opm::min(1.0, xlCO2));
        xgH2O = Opm::max(0.0, Opm::min(1.0, xgH2O));

        xlH2O = 1.0 - xlCO2;
        xgCO2 = 1.0 - xgH2O;

        if (compIdx == BrineIdx) {
            Scalar phigH2O = 1.0;
            return phigH2O * xgH2O / xlH2O;
        }
        else {
            assert(compIdx == CO2Idx);

            Scalar phigCO2 = 1.0;
            return phigCO2 * xgCO2 / xlCO2;
        };
    }

    /*!
     * \copydoc BaseFluidSystem::diffusionCoefficient
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned /*compIdx*/)
    {
        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == liquidPhaseIdx)
            return BinaryCoeffBrineCO2::liquidDiffCoeff(temperature, pressure);

        assert(phaseIdx == gasPhaseIdx);
        return BinaryCoeffBrineCO2::gasDiffCoeff(temperature, pressure);
    }

    /*!
     * \copydoc BaseFluidSystem::enthalpy
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == liquidPhaseIdx) {
            const LhsEval& XlCO2 = Opm::decay<LhsEval>(fluidState.massFraction(phaseIdx, CO2Idx));
            const LhsEval& result = liquidEnthalpyBrineCO2_(temperature,
                                                            pressure,
                                                            Brine_IAPWS::salinity,
                                                            XlCO2);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            const LhsEval& XCO2 = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, CO2Idx));
            const LhsEval& XBrine = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, BrineIdx));

            LhsEval result = 0;
            result += XBrine * Brine::gasEnthalpy(temperature, pressure);
            result += XCO2 * CO2::gasEnthalpy(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    /*!
     * \copydoc BaseFluidSystem::thermalConductivity
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& /*fluidState*/,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        // TODO way too simple!
        if (phaseIdx == liquidPhaseIdx)
            return  0.6; // conductivity of water[W / (m K ) ]

        // gas phase
        return 0.025; // conductivity of air [W / (m K ) ]
    }

    /*!
     * \copydoc BaseFluidSystem::heatCapacity
     *
     * We employ the heat capacity of the pure phases.
     *
     * Todo: Include compositional effects.
     *
     * \param fluidState An arbitrary fluid state
     * \param paramCache The object which caches parameters which are expensive to compute
     * \param phaseIdx The index of the fluid phase to consider
     * \tparam FluidState the fluid state class
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const LhsEval& temperature = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if(phaseIdx == liquidPhaseIdx)
            return H2O::liquidHeatCapacity(temperature, pressure);
        else
            return CO2::gasHeatCapacity(temperature, pressure);
    }

private:
    template <class LhsEval>
    static LhsEval gasDensity_(const LhsEval& T,
                               const LhsEval& pg,
                               const LhsEval& xgH2O,
                               const LhsEval& xgCO2)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pg);
        Valgrind::CheckDefined(xgH2O);
        Valgrind::CheckDefined(xgCO2);

        return CO2::gasDensity(T, pg);
    }

    /***********************************************************************/
    /*                                                                     */
    /* Total brine density with dissolved CO2                              */
    /* rho_{b,CO2} = rho_w + contribution(salt) + contribution(CO2)        */
    /*                                                                     */
    /***********************************************************************/
    template <class LhsEval>
    static LhsEval liquidDensity_(const LhsEval& T,
                                  const LhsEval& pl,
                                  const LhsEval& xlH2O,
                                  const LhsEval& xlCO2)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(xlH2O);
        Valgrind::CheckDefined(xlCO2);

        if(T < 273.15) {
            std::ostringstream oss;
            oss << "Liquid density for Brine and CO2 is only "
                "defined above 273.15K (is "<<T<<"K)";
            throw NumericalIssue(oss.str());
        }
        if(pl >= 2.5e8) {
            std::ostringstream oss;
            oss << "Liquid density for Brine and CO2 is only "
                "defined below 250MPa (is "<<pl<<"Pa)";
            throw NumericalIssue(oss.str());
        }

        const LhsEval& rho_brine = Brine::liquidDensity(T, pl);
        const LhsEval& rho_pure = H2O::liquidDensity(T, pl);
        const LhsEval& rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
        const LhsEval& contribCO2 = rho_lCO2 - rho_pure;

        return rho_brine + contribCO2;
    }

    template <class LhsEval>
    static LhsEval liquidDensityWaterCO2_(const LhsEval& temperature,
                                          const LhsEval& pl,
                                          const LhsEval& /*xlH2O*/,
                                          const LhsEval& xlCO2)
    {
        Scalar M_CO2 = CO2::molarMass();
        Scalar M_H2O = H2O::molarMass();

        const LhsEval& tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const LhsEval& rho_pure = H2O::liquidDensity(temperature, pl);
        // calculate the mole fraction of CO2 in the liquid. note that xlH2O is available
        // as a function parameter, but in the case of a pure gas phase the value of M_T
        // for the virtual liquid phase can become very large
        const LhsEval xlH2O = 1.0 - xlCO2;
        const LhsEval& M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const LhsEval& V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }

    template <class LhsEval>
    static LhsEval liquidEnthalpyBrineCO2_(const LhsEval& T,
                                           const LhsEval& p,
                                           Scalar S, // salinity
                                           const LhsEval& X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        /* same function as enthalpy_brine, only extended by CO2 content */

        /*Numerical coefficents from PALLISER*/
        static Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static Scalar a[4][3] = {
            { 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        LhsEval theta, h_NaCl;
        LhsEval h_ls1, d_h;
        LhsEval delta_h;
        LhsEval delta_hCO2, hg, hw;

        theta = T - 273.15;

        // Regularization
        Scalar scalarTheta = Opm::scalarValue(theta);
        Scalar S_lSAT = f[0] + scalarTheta*(f[1] + scalarTheta*(f[2] + scalarTheta*f[3]));
        if (S > S_lSAT)
            S = S_lSAT;

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        Scalar m = 1E3/58.44 * S/(1-S);
        int i = 0;
        int j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * Opm::pow(theta, static_cast<Scalar>(i)) * std::pow(m, j);
            }
        }
        /* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without CO2 */
        h_ls1 =(1-S)*hw + S*h_NaCl + S*delta_h; /* kJ/kg */

        /* heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
           In the relevant temperature ranges CO2 dissolution is
           exothermal */
        delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

        /* enthalpy contribution of CO2 (kJ/kg) */
        hg = CO2::gasEnthalpy(T, p)/1E3 + delta_hCO2;

        /* Enthalpy of brine with dissolved CO2 */
        return (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1E3; /*J/kg*/
    }
};

} // namespace FluidSystems
} // namespace Opm

#endif
