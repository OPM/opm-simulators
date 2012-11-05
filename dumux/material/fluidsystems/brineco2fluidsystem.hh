// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \copydoc Dumux::FluidSystems::BrineCO2
 */
#ifndef DUMUX_BRINE_CO2_SYSTEM_HH
#define DUMUX_BRINE_CO2_SYSTEM_HH

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/brine.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/binarycoefficients/h2o_co2.hh>
#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux {
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
    typedef Dumux::H2O<Scalar> H2O_IAPWS;
    typedef Dumux::Brine<Scalar, H2O_IAPWS> Brine_IAPWS;
    typedef Dumux::TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;
    typedef Dumux::TabulatedComponent<Scalar, Brine_IAPWS> Brine_Tabulated;

    typedef H2O_Tabulated H2O;

public:
    //! The binary coefficients for brine and CO2 used by this fluid system
    typedef Dumux::BinaryCoeff::Brine_CO2<Scalar, CO2Tables> BinaryCoeffBrineCO2;

    //! \copydoc BaseFluidSystem::ParameterCache
    typedef Dumux::NullParameterCache ParameterCache;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! The number of phases considered by the fluid system
    static const int numPhases = 2;

    //! The index of the liquid phase
    static const int lPhaseIdx = 0;
    //! The index of the gas phase
    static const int gPhaseIdx = 1;

    /*!
     * \copydoc BaseFluidSystem::phaseName
     */
    static const char *phaseName(int phaseIdx)
    {
        static const char *name[] = {
            "l",
            "g"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \copydoc BaseFluidSystem::isLiquid
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != gPhaseIdx;
    }

    /*!
     * \copydoc BaseFluidSystem::isIdealGas
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == gPhaseIdx)
            return CO2::gasIsIdeal();
        return false;
    }

    /*!
     * \copydoc BaseFluidSystem::isIdealMixture
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
     * \copydoc BaseFluidSystem::isCompressible
     */
    static bool isCompressible(int phaseIdx)
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
    typedef Dumux::CO2<Scalar, CO2Tables> CO2;

    /*!
     * \copydoc BaseFluidSystem::componentName
     */
    static const char *componentName(int compIdx)
    {
        static const char *name[] = {
            Brine::name(),
            CO2::name(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \copydoc BaseFluidSystem::molarMass
     */
    static Scalar molarMass(int compIdx)
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
        init(/*startTemp=*/273.15, /*endTemp=*/623.15, /*tempSteps=*/100,
             /*startPressure=*/1e4, /*endPressure=*/40e6, /*pressureSteps=*/200);
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
    static void init(Scalar tempMin, Scalar tempMax, int nTemp,
                     Scalar pressMin, Scalar pressMax, int nPress)
    {
        int myRank = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif
        if (H2O::isTabulated) {
            if (myRank == 0)
                std::cout << "Initializing tables for the pure-water properties.\n";
            H2O_Tabulated::init(tempMin, tempMax, nTemp,
                                pressMin, pressMax, nPress);
        }

        // set the salinity of brine to the one used by the CO2 tables
        Brine_IAPWS::salinity = CO2Tables::brineSalinity;

        if (Brine::isTabulated) {
            if (myRank == 0)
                std::cout << "Initializing tables for the brine fluid properties.\n";
            Brine_Tabulated::init(tempMin, tempMax, nTemp,
                                  pressMin, pressMax, nPress);
        }
    }

    /*!
     * \copydoc BaseFluidSystem::density
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx) {
            // use normalized composition for to calculate the density
            // (the relations don't seem to take non-normalized
            // compositions too well...)
            Scalar xlBrine = std::min(1.0, std::max(0.0, fluidState.moleFraction(lPhaseIdx, BrineIdx)));
            Scalar xlCO2 = std::min(1.0, std::max(0.0, fluidState.moleFraction(lPhaseIdx, CO2Idx)));
            Scalar sumx = xlBrine + xlCO2;
            xlBrine /= sumx;
            xlCO2 /= sumx;

            Scalar result = liquidDensity_(temperature,
                                           pressure,
                                           xlBrine,
                                           xlCO2);

            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == gPhaseIdx);

            // use normalized composition for to calculate the density
            // (the relations don't seem to take non-normalized
            // compositions too well...)
            Scalar xgBrine = std::min(1.0, std::max(0.0, fluidState.moleFraction(gPhaseIdx, BrineIdx)));
            Scalar xgCO2 = std::min(1.0, std::max(0.0, fluidState.moleFraction(gPhaseIdx, CO2Idx)));
            Scalar sumx = xgBrine + xgCO2;
            xgBrine /= sumx;
            xgCO2 /= sumx;

            Scalar result = gasDensity_(temperature,
                                        pressure,
                                        xgBrine,
                                        xgCO2);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    /*!
     * \copydoc BaseFluidSystem::viscosity
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx) {
            // assume pure brine for the liquid phase. TODO: viscosity
            // of mixture
            Scalar result = Brine::liquidViscosity(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }

        assert(phaseIdx == gPhaseIdx);
        Scalar result = CO2::gasViscosity(temperature, pressure);
        Valgrind::CheckDefined(result);
        return result;
    }

    /*!
     * \copydoc BaseFluidSystem::fugacityCoefficient
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == gPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        assert(temperature > 0);
        assert(pressure > 0);

        // calulate the equilibrium composition for the given
        // temperature and pressure. TODO: calculateMoleFractions()
        // could use some cleanup.
        Scalar xlH2O, xgH2O;
        Scalar xlCO2, xgCO2;
        BinaryCoeffBrineCO2::calculateMoleFractions(temperature,
                                                    pressure,
                                                    Brine_IAPWS::salinity,
                                                    /*knownPhaseIdx=*/-1,
                                                    xlCO2,
                                                    xgH2O);

        // normalize the phase compositions
        xlCO2 = std::max(0.0, std::min(1.0, xlCO2));
        xgH2O = std::max(0.0, std::min(1.0, xgH2O));

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
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == lPhaseIdx)
            return BinaryCoeffBrineCO2::liquidDiffCoeff(temperature, pressure);

        assert(phaseIdx == gPhaseIdx);
        return BinaryCoeffBrineCO2::gasDiffCoeff(temperature, pressure);
    }

    /*!
     * \copydoc BaseFluidSystem::enthalpy
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx) {
            Scalar XlCO2 = fluidState.massFraction(phaseIdx, CO2Idx);
            Scalar result = liquidEnthalpyBrineCO2_(temperature,
                                                    pressure,
                                                    Brine_IAPWS::salinity,
                                                    XlCO2);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            Scalar XCO2 = fluidState.massFraction(gPhaseIdx, CO2Idx);
            Scalar XBrine = fluidState.massFraction(gPhaseIdx, BrineIdx);

            Scalar result = 0;
            result += XBrine * Brine::gasEnthalpy(temperature, pressure);
            result += XCO2 * CO2::gasEnthalpy(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    /*!
     * \copydoc BaseFluidSystem::thermalConductivity
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        // TODO way too simple!
        if (phaseIdx == lPhaseIdx)
            return  0.6; // conductivity of water[W / (m K ) ]

        // gas phase
        return 0.025; // conductivity of air [W / (m K ) ]
    }

private:
    static Scalar gasDensity_(Scalar T,
                              Scalar pg,
                              Scalar xgH2O,
                              Scalar xgCO2)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pg);
        Valgrind::CheckDefined(xgH2O);
        Valgrind::CheckDefined(xgCO2);

        Scalar gasDensity = CO2::gasDensity(T, pg);
        return gasDensity;
    }

    /***********************************************************************/
    /*                                                                     */
    /* Total brine density with dissolved CO2                              */
    /* rho_{b,CO2} = rho_w + contribution(salt) + contribution(CO2)        */
    /*                                                                     */
    /***********************************************************************/
    static Scalar liquidDensity_(Scalar T,
                                 Scalar pl,
                                 Scalar xlH2O,
                                 Scalar xlCO2)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(xlH2O);
        Valgrind::CheckDefined(xlCO2);

        if(T < 273.15) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined above 273.15K (is " << T << "K)");
        }
        if(pl >= 2.5e8) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined below 250MPa (is " << pl << "Pa)");
        }

        Scalar rho_brine = Brine::liquidDensity(T, pl);
        Scalar rho_pure = H2O::liquidDensity(T, pl);
        Scalar rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
        Scalar contribCO2 = rho_lCO2 - rho_pure;

        return rho_brine + contribCO2;
    }

    static Scalar liquidDensityWaterCO2_(Scalar temperature,
                                         Scalar pl,
                                         Scalar xlH2O,
                                         Scalar xlCO2)
    {
        const Scalar M_CO2 = CO2::molarMass();
        const Scalar M_H2O = H2O::molarMass();

        const Scalar tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const Scalar rho_pure = H2O::liquidDensity(temperature, pl);
        xlH2O = 1.0 - xlCO2; // xlH2O is available, but in case of a pure gas phase
                             // the value of M_T for the virtual liquid phase can become very large
        const Scalar M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const Scalar V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }

    static Scalar liquidEnthalpyBrineCO2_(Scalar T,
                                          Scalar p,
                                          Scalar S,
                                          Scalar X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        /* same function as enthalpy_brine, only extended by CO2 content */

        /*Numerical coefficents from PALLISER*/
        static const Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = {
            { 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        Scalar theta, h_NaCl;
        Scalar m, h_ls1, d_h;
        Scalar S_lSAT, delta_h;
        int i, j;
        Scalar delta_hCO2, hg, hw;

        theta = T - 273.15;

        S_lSAT = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;
        /*Regularization*/
        if (S>S_lSAT) {
            S = S_lSAT;
        }

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(S/(1-S));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
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

} // end namepace
} // end namepace

#endif
