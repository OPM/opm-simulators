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
 * \copydoc Dumux::BinaryCoeff::Brine_CO2
 */
#ifndef DUMUX_BINARY_COEFF_BRINE_CO2_HH
#define DUMUX_BINARY_COEFF_BRINE_CO2_HH

#include <dumux/material/components/brine.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for brine and CO2.
 */
template<class Scalar, class CO2Tables, bool verbose = true>
class Brine_CO2 {
    typedef Dumux::H2O<Scalar> H2O;
    typedef Dumux::CO2<Scalar, CO2Tables> CO2;
    typedef Dumux::Brine<Scalar,H2O> Brine;
    typedef Dumux::IdealGas<Scalar> IdealGas;
    static const int lPhaseIdx = 0; // index of the liquid phase
    static const int gPhaseIdx = 1; // index of the gas phase

public:
    /*!
     * \brief Binary diffusion coefficent [m^2/s] of water in the CO2 phase.
     *
     * According to "Diffusion of Water in Liquid and Supercritical Carbon Dioxide: An NMR Study",Bin Xu et al., 2002
     * \param temperature the temperature [K]
     * \param pressure the phase pressure [Pa]
     */
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //Diffusion coefficient of water in the CO2 phase
        Scalar const PI=3.141593;
        Scalar const k = 1.3806504e-23; // Boltzmann constant
        Scalar const c = 4; // slip parameter, can vary between 4 (slip condition) and 6 (stick condition)
        Scalar const R_h = 1.72e-10; // hydrodynamic radius of the solute
        Scalar mu = CO2::gasViscosity(temperature, pressure); // CO2 viscosity
        Scalar D = k / (c * PI * R_h) * (temperature / mu);
        return D;
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] of CO2 in the brine phase.
     *
     * \param temperature the temperature [K]
     * \param pressure the phase pressure [Pa]
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //Diffusion coefficient of CO2 in the brine phase
        return 2e-9;
    }

    /*!
     * \brief Returns the _mol_ (!) fraction of CO2 in the liquid
     *        phase and the mol_ (!) fraction of H2O in the gas phase
     *        for a given temperature, pressure, CO2 density and brine
     *        salinity.
     *
     *        Implemented according to "Spycher and Pruess 2005"
     *        applying the activity coefficient expression of "Duan and Sun 2003"
     *        and the correlations for pure water given in "Spycher, Pruess and Ennis-King 2003"
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     * \param salinity the salinity [kg NaCl / kg solution]
     * \param knownPhaseIdx indicates which phases are present
     * \param xlCO2 mole fraction of CO2 in brine [mol/mol]
     * \param ygH2O mole fraction of water in the gas phase [mol/mol]
     */
    static void calculateMoleFractions(const Scalar temperature,
                                       const Scalar pg,
                                       const Scalar salinity,
                                       const int knownPhaseIdx,
                                       Scalar &xlCO2,
                                       Scalar &ygH2O)
    {
        Scalar A = computeA_(temperature, pg);

        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = salinityToMolFrac_(salinity);

        // if both phases are present the mole fractions in each phase can be calculate
        // with the mutual solubility function
        if (knownPhaseIdx < 0) {
            Scalar molalityNaCl = molFracToMolality_(x_NaCl); // molality of NaCl //CHANGED
            Scalar m0_CO2 = molalityCO2inPureWater_(temperature, pg); // molality of CO2 in pure water
            Scalar gammaStar = activityCoefficient_(temperature, pg, molalityNaCl);// activity coefficient of CO2 in brine
            Scalar m_CO2 = m0_CO2 / gammaStar; // molality of CO2 in brine
            xlCO2 = m_CO2 / (molalityNaCl + 55.508 + m_CO2); // mole fraction of CO2 in brine
            ygH2O = A * (1 - xlCO2 - x_NaCl); // mole fraction of water in the gas phase
        }

        // if only liquid phase is present the mole fraction of CO2 in brine is given and
        // and the virtual equilibrium mole fraction of water in the non-existing gas phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == lPhaseIdx) {
            ygH2O = A * (1 - xlCO2 - x_NaCl);

        }

        // if only gas phase is present the mole fraction of water in the gas phase is given and
        // and the virtual equilibrium mole fraction of CO2 in the non-existing liquid phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == gPhaseIdx) {
            //y_H2o = fluidstate.
            xlCO2 = 1 - x_NaCl - ygH2O / A;
        }

    }

    /*!
     * \brief Henry coefficent \f$\mathrm{[N/m^2]}\f$ for CO2 in brine.
     */
    static Scalar henry(Scalar temperature)
    { return fugacityCoefficientCO2(temperature, /*pressure=*/1e5)*1e5; }

    /*!
     * \brief Returns the fugacity coefficient of the CO2 component in a water-CO2 mixture
     *
     * (given in Spycher, Pruess and Ennis-King (2003))
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar fugacityCoefficientCO2(Scalar T, Scalar pg)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pg);

        Scalar V = 1 / (CO2::gasDensity(T, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        Scalar a_CO2 = (7.54e7 - 4.13e4 * T); // mixture parameter of  Redlich-Kwong equation
        static const Scalar b_CO2 = 27.8; // mixture parameter of Redlich-Kwong equation
        static const Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Scalar lnPhiCO2, phiCO2;

        lnPhiCO2 = std::log(V / (V - b_CO2));
        lnPhiCO2 += b_CO2 / (V - b_CO2);
        lnPhiCO2 -= 2 * a_CO2 / (R * std::pow(T, 1.5) * b_CO2) * log((V + b_CO2) / V);
        lnPhiCO2 +=
            a_CO2 * b_CO2
            / (R
               * std::pow(T, 1.5)
               * b_CO2
               * b_CO2)
            * (std::log((V + b_CO2) / V)
               - b_CO2 / (V + b_CO2));
        lnPhiCO2 -= std::log(pg_bar * V / (R * T));

        phiCO2 = exp(lnPhiCO2); // fugacity coefficient of CO2
        return phiCO2;

    }

    /*!
     * \brief Returns the fugacity coefficient of the H2O component in a water-CO2 mixture
     *
     * (given in Spycher, Pruess and Ennis-King (2003))
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar fugacityCoefficientH2O(Scalar T, Scalar pg)
    {
        Scalar V = 1 / (CO2::gasDensity(T, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        Scalar a_CO2 = (7.54e7 - 4.13e4 * T);// mixture parameter of  Redlich-Kwong equation
        static const Scalar a_CO2_H2O = 7.89e7;// mixture parameter of Redlich-Kwong equation
        static const Scalar b_CO2 = 27.8;// mixture parameter of Redlich-Kwong equation
        static const Scalar b_H2O = 18.18;// mixture parameter of Redlich-Kwong equation
        static const Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Scalar lnPhiH2O, phiH2O;

        lnPhiH2O = log(V / (V - b_CO2)) + b_H2O / (V - b_CO2) - 2 * a_CO2_H2O
                / (R * pow(T, 1.5) * b_CO2) * log((V + b_CO2) / V) + a_CO2
                * b_H2O / (R * pow(T, 1.5) * b_CO2 * b_CO2) * (log((V + b_CO2)
                / V) - b_CO2 / (V + b_CO2)) - log(pg_bar * V / (R * T));
        phiH2O = exp(lnPhiH2O); // fugacity coefficient of H2O
        return phiH2O;
    }

private:
    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
     *
     * \param salinity the salinity [kg NaCl / kg solution]
     */
    static Scalar salinityToMolFrac_(Scalar salinity) {

        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = salinity;
        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return x_NaCl;
    }

    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction (mol NaCl / mol solution)
     *
     * \param x_NaCl mole fraction of NaCL in brine [mol/mol]
     */
    static Scalar molFracToMolality_(Scalar x_NaCl) {

        // conversion from mol fraction to molality (dissolved CO2 neglected)
        const Scalar mol_NaCl = 55.508 * x_NaCl / (1 - x_NaCl);

        return mol_NaCl;
    }

    /*!
     * \brief Returns the equilibrium molality of CO2 (mol CO2 / kg water) for a
     * CO2-water mixture at a given pressure and temperature
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar molalityCO2inPureWater_(Scalar temperature, Scalar pg) {
        Scalar A = computeA_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        Scalar B = computeB_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        Scalar yH2OinGas = (1 - B) / (1. / A - B); // equilibrium mol fraction of H2O in the gas phase
        Scalar xCO2inWater = B * (1 - yH2OinGas); // equilibrium mol fraction of CO2 in the water phase
        Scalar molalityCO2 = (xCO2inWater * 55.508) / (1 - xCO2inWater); // CO2 molality
        return molalityCO2;
    }

    /*!
     * \brief Returns the activity coefficient of CO2 in brine for a
     *           molal description. According to "Duan and Sun 2003"
     *           given in "Spycher and Pruess 2005"
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     * \param molalityNaCl molality of NaCl (mol NaCl / kg water)
     */
    static Scalar activityCoefficient_(Scalar temperature,
                                       Scalar pg,
                                       Scalar molalityNaCl)
    {
        Scalar lambda = computeLambda_(temperature, pg); // lambda_{CO2-Na+}
        Scalar xi = computeXi_(temperature, pg); // Xi_{CO2-Na+-Cl-}
        Scalar lnGammaStar = 2 * lambda * molalityNaCl + xi * molalityNaCl
                * molalityNaCl;
        Scalar gammaStar = exp(lnGammaStar);
        return gammaStar; // molal activity coefficient of CO2 in brine
    }

    /*!
     * \brief Returns the paramater A for the calculation of
     * them mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003)
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar computeA_(Scalar T, Scalar pg)
    {
        Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_H2O = 18.1; // average partial molar volume of H2O [cm^3/mol]
        const Scalar R = IdealGas::R * 10;
        Scalar k0_H2O = equilibriumConstantH2O_(T); // equilibrium constant for H2O at 1 bar
        Scalar phi_H2O = fugacityCoefficientH2O(T, pg); // fugacity coefficient of H2O for the water-CO2 system
        Scalar pg_bar = pg / 1.e5;
        Scalar A = k0_H2O / (phi_H2O * pg_bar) * exp(deltaP * v_av_H2O / (R * T));
        return A;

    }

    /*!
     * \brief Returns the paramater B for the calculation of
     * the mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003)
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar computeB_(Scalar T, Scalar pg) {
        Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_CO2 = 32.6; // average partial molar volume of CO2 [cm^3/mol]
        const Scalar R = IdealGas::R * 10;
        Scalar k0_CO2 = equilibriumConstantCO2_(T); // equilibrium constant for CO2 at 1 bar
        Scalar phi_CO2 = fugacityCoefficientCO2(T, pg); // fugacity coefficient of CO2 for the water-CO2 system
        Scalar pg_bar = pg / 1.e5;
        Scalar B = phi_CO2 * pg_bar / (55.508 * k0_CO2) * exp(-(deltaP
                * v_av_CO2) / (R * T));
        return B;
    }

    /*!
     * \brief Returns the parameter lambda, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005)
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar computeLambda_(Scalar T, Scalar pg)
    {
        Scalar lambda;
        static const Scalar c[6] = { -0.411370585, 6.07632013E-4, 97.5347708,
                -0.0237622469, 0.0170656236, 1.41335834E-5 };

        Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        lambda = c[0] + c[1] * T + c[2] / T + c[3] * pg_bar / T + c[4] * pg_bar
                / (630.0 - T) + c[5] * T * std::log(pg_bar);

        return lambda;
    }

    /*!
     * \brief Returns the parameter xi, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005)
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar computeXi_(Scalar T, Scalar pg)
    {
        Scalar xi;
        static const Scalar c[4] = { 3.36389723E-4, -1.98298980E-5,
                2.12220830E-3, -5.24873303E-3 };

        Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        xi = c[0] + c[1] * T + c[2] * pg_bar / T + c[3] * pg_bar / (630.0 - T);

        return xi;
    }

    /*!
     * \brief Returns the equilibrium constant for CO2, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003)
     * \param T the temperature [K]
     */
    static Scalar equilibriumConstantCO2_(Scalar T)
    {
        Scalar TinC = T - 273.15; //temperature in °C
        static const Scalar c[3] = { 1.189, 1.304e-2, -5.446e-5 };
        Scalar logk0_CO2 = c[0] + c[1] * TinC + c[2] * TinC * TinC;
        Scalar k0_CO2 = pow(10, logk0_CO2);
        return k0_CO2;
    }

    /*!
     * \brief Returns the equilibrium constant for H2O, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003)
     * \param T the temperature [K]
     */
    static Scalar equilibriumConstantH2O_(Scalar T)
    {
        Scalar TinC = T - 273.15; //temperature in °C
        static const Scalar c[4] = { -2.209, 3.097e-2, -1.098e-4, 2.048e-7 };
        Scalar logk0_H2O = c[0] + c[1] * TinC + c[2] * TinC * TinC + c[3]
                * TinC * TinC * TinC;
        Scalar k0_H2O = pow(10, logk0_H2O);
        return k0_H2O;
    }

};

} // end namepace BinaryCoeff
} // end namepace Dumux

#endif
