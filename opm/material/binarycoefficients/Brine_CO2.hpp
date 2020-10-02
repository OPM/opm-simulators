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
 * \copydoc Opm::BinaryCoeff::Brine_CO2
 */
#ifndef OPM_BINARY_COEFF_BRINE_CO2_HPP
#define OPM_BINARY_COEFF_BRINE_CO2_HPP

#include <opm/material/IdealGas.hpp>

namespace Opm {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for brine and CO2.
 */
template<class Scalar, class H2O, class CO2, bool verbose = true>
class Brine_CO2 {
    typedef Opm::IdealGas<Scalar> IdealGas;
    static const int liquidPhaseIdx = 0; // index of the liquid phase
    static const int gasPhaseIdx = 1; // index of the gas phase

public:
    /*!
     * \brief Binary diffusion coefficent [m^2/s] of water in the CO2 phase.
     *
     * According to "Diffusion of Water in Liquid and Supercritical Carbon Dioxide: An NMR Study",Bin Xu et al., 2002
     * \param temperature the temperature [K]
     * \param pressure the phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        //Diffusion coefficient of water in the CO2 phase
        Scalar k = 1.3806504e-23; // Boltzmann constant
        Scalar c = 4; // slip parameter, can vary between 4 (slip condition) and 6 (stick condition)
        Scalar R_h = 1.72e-10; // hydrodynamic radius of the solute
        const Evaluation& mu = CO2::gasViscosity(temperature, pressure); // CO2 viscosity
        return k / (c * M_PI * R_h) * (temperature / mu);
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] of CO2 in the brine phase.
     *
     * \param temperature the temperature [K]
     * \param pressure the phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
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
    template <class Evaluation>
    static void calculateMoleFractions(const Evaluation& temperature,
                                       const Evaluation& pg,
                                       Scalar salinity,
                                       const int knownPhaseIdx,
                                       Evaluation& xlCO2,
                                       Evaluation& ygH2O)
    {
        Evaluation A = computeA_(temperature, pg);

        /* salinity: conversion from mass fraction to mol fraction */
        Scalar x_NaCl = salinityToMolFrac_(salinity);

        // if both phases are present the mole fractions in each phase can be calculate
        // with the mutual solubility function
        if (knownPhaseIdx < 0) {
            Scalar molalityNaCl = moleFracToMolality_(x_NaCl); // molality of NaCl //CHANGED
            Evaluation m0_CO2 = molalityCO2inPureWater_(temperature, pg); // molality of CO2 in pure water
            Evaluation gammaStar = activityCoefficient_(temperature, pg, molalityNaCl);// activity coefficient of CO2 in brine
            Evaluation m_CO2 = m0_CO2 / gammaStar; // molality of CO2 in brine
            xlCO2 = m_CO2 / (molalityNaCl + 55.508 + m_CO2); // mole fraction of CO2 in brine
            ygH2O = A * (1 - xlCO2 - x_NaCl); // mole fraction of water in the gas phase
        }

        // if only liquid phase is present the mole fraction of CO2 in brine is given and
        // and the virtual equilibrium mole fraction of water in the non-existing gas phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == liquidPhaseIdx)
            ygH2O = A * (1 - xlCO2 - x_NaCl);

        // if only gas phase is present the mole fraction of water in the gas phase is given and
        // and the virtual equilibrium mole fraction of CO2 in the non-existing liquid phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == gasPhaseIdx)
            //y_H2o = fluidstate.
            xlCO2 = 1 - x_NaCl - ygH2O / A;
    }

    /*!
     * \brief Henry coefficent \f$\mathrm{[N/m^2]}\f$ for CO2 in brine.
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& temperature)
    { return fugacityCoefficientCO2(temperature, /*pressure=*/1e5)*1e5; }

    /*!
     * \brief Returns the fugacity coefficient of the CO2 component in a water-CO2 mixture
     *
     * (given in Spycher, Pruess and Ennis-King (2003))
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation fugacityCoefficientCO2(const Evaluation& temperature, const Evaluation& pg)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(pg);

        Evaluation V = 1 / (CO2::gasDensity(temperature, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        Evaluation pg_bar = pg / 1.e5; // gas phase pressure in bar
        Evaluation a_CO2 = (7.54e7 - 4.13e4 * temperature); // mixture parameter of  Redlich-Kwong equation
        Scalar b_CO2 = 27.8; // mixture parameter of Redlich-Kwong equation
        Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Evaluation lnPhiCO2;

        lnPhiCO2 = Opm::log(V / (V - b_CO2));
        lnPhiCO2 += b_CO2 / (V - b_CO2);
        lnPhiCO2 -= 2 * a_CO2 / (R * Opm::pow(temperature, 1.5) * b_CO2) * log((V + b_CO2) / V);
        lnPhiCO2 +=
            a_CO2 * b_CO2
            / (R
               * Opm::pow(temperature, 1.5)
               * b_CO2
               * b_CO2)
            * (Opm::log((V + b_CO2) / V)
               - b_CO2 / (V + b_CO2));
        lnPhiCO2 -= Opm::log(pg_bar * V / (R * temperature));

        return Opm::exp(lnPhiCO2); // fugacity coefficient of CO2
    }

    /*!
     * \brief Returns the fugacity coefficient of the H2O component in a water-CO2 mixture
     *
     * (given in Spycher, Pruess and Ennis-King (2003))
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation fugacityCoefficientH2O(const Evaluation& temperature, const Evaluation& pg)
    {
        const Evaluation& V = 1 / (CO2::gasDensity(temperature, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        const Evaluation& pg_bar = pg / 1.e5; // gas phase pressure in bar
        const Evaluation& a_CO2 = (7.54e7 - 4.13e4 * temperature);// mixture parameter of  Redlich-Kwong equation
        Scalar a_CO2_H2O = 7.89e7;// mixture parameter of Redlich-Kwong equation
        Scalar b_CO2 = 27.8;// mixture parameter of Redlich-Kwong equation
        Scalar b_H2O = 18.18;// mixture parameter of Redlich-Kwong equation
        Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Evaluation lnPhiH2O;

        lnPhiH2O =
            Opm::log(V/(V - b_CO2))
            + b_H2O/(V - b_CO2) - 2*a_CO2_H2O
            / (R*Opm::pow(temperature, 1.5)*b_CO2)*Opm::log((V + b_CO2)/V)
            + a_CO2*b_H2O/(R*Opm::pow(temperature, 1.5)*b_CO2*b_CO2)
            *(Opm::log((V + b_CO2)/V) - b_CO2/(V + b_CO2))
            - Opm::log(pg_bar*V/(R*temperature));
        return Opm::exp(lnPhiH2O); // fugacity coefficient of H2O
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
    static Scalar moleFracToMolality_(Scalar x_NaCl)
    {
        // conversion from mol fraction to molality (dissolved CO2 neglected)
        return 55.508 * x_NaCl / (1 - x_NaCl);
    }

    /*!
     * \brief Returns the equilibrium molality of CO2 (mol CO2 / kg water) for a
     * CO2-water mixture at a given pressure and temperature
     *
     * \param temperature The temperature [K]
     * \param pg The gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation molalityCO2inPureWater_(const Evaluation& temperature, const Evaluation& pg)
    {
        const Evaluation& A = computeA_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        const Evaluation& B = computeB_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        const Evaluation& yH2OinGas = (1 - B) / (1. / A - B); // equilibrium mol fraction of H2O in the gas phase
        const Evaluation& xCO2inWater = B * (1 - yH2OinGas); // equilibrium mol fraction of CO2 in the water phase
        return (xCO2inWater * 55.508) / (1 - xCO2inWater); // CO2 molality
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
    template <class Evaluation>
    static Evaluation activityCoefficient_(const Evaluation& temperature,
                                           const Evaluation& pg,
                                           Scalar molalityNaCl)
    {
        const Evaluation& lambda = computeLambda_(temperature, pg); // lambda_{CO2-Na+}
        const Evaluation& xi = computeXi_(temperature, pg); // Xi_{CO2-Na+-Cl-}
        const Evaluation& lnGammaStar =
            2*molalityNaCl*lambda + xi*molalityNaCl*molalityNaCl;
        return Opm::exp(lnGammaStar);
    }

    /*!
     * \brief Returns the paramater A for the calculation of
     * them mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003)
     *
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation computeA_(const Evaluation& temperature, const Evaluation& pg)
    {
        const Evaluation& deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        Scalar v_av_H2O = 18.1; // average partial molar volume of H2O [cm^3/mol]
        Scalar R = IdealGas::R * 10;
        const Evaluation& k0_H2O = equilibriumConstantH2O_(temperature); // equilibrium constant for H2O at 1 bar
        const Evaluation& phi_H2O = fugacityCoefficientH2O(temperature, pg); // fugacity coefficient of H2O for the water-CO2 system
        const Evaluation& pg_bar = pg / 1.e5;
        return k0_H2O/(phi_H2O*pg_bar)*Opm::exp(deltaP*v_av_H2O/(R*temperature));
    }

    /*!
     * \brief Returns the paramater B for the calculation of
     * the mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003)
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation computeB_(const Evaluation& temperature, const Evaluation& pg)
    {
        const Evaluation& deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_CO2 = 32.6; // average partial molar volume of CO2 [cm^3/mol]
        const Scalar R = IdealGas::R * 10;
        const Evaluation& k0_CO2 = equilibriumConstantCO2_(temperature); // equilibrium constant for CO2 at 1 bar
        const Evaluation& phi_CO2 = fugacityCoefficientCO2(temperature, pg); // fugacity coefficient of CO2 for the water-CO2 system
        const Evaluation& pg_bar = pg / 1.e5;
        return phi_CO2*pg_bar/(55.508*k0_CO2)*Opm::exp(-(deltaP*v_av_CO2)/(R*temperature));
    }

    /*!
     * \brief Returns the parameter lambda, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005)
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation computeLambda_(const Evaluation& temperature, const Evaluation& pg)
    {
        static const Scalar c[6] =
            { -0.411370585, 6.07632013E-4, 97.5347708, -0.0237622469, 0.0170656236, 1.41335834E-5 };

        Evaluation pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return
            c[0]
            + c[1]*temperature
            + c[2]/temperature
            + c[3]*pg_bar/temperature
            + c[4]*pg_bar/(630.0 - temperature)
            + c[5]*temperature*Opm::log(pg_bar);
    }

    /*!
     * \brief Returns the parameter xi, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005)
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    template <class Evaluation>
    static Evaluation computeXi_(const Evaluation& temperature, const Evaluation& pg)
    {
        static const Scalar c[4] =
            { 3.36389723E-4, -1.98298980E-5, 2.12220830E-3, -5.24873303E-3 };

        Evaluation pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return c[0] + c[1]*temperature + c[2]*pg_bar/temperature + c[3]*pg_bar/(630.0 - temperature);
    }

    /*!
     * \brief Returns the equilibrium constant for CO2, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003)
     * \param temperature the temperature [K]
     */
    template <class Evaluation>
    static Evaluation equilibriumConstantCO2_(const Evaluation& temperature)
    {
        Evaluation temperatureCelcius = temperature - 273.15;
        static const Scalar c[3] = { 1.189, 1.304e-2, -5.446e-5 };
        Evaluation logk0_CO2 = c[0] + temperatureCelcius*(c[1] + temperatureCelcius*c[2]);
        Evaluation k0_CO2 = Opm::pow(10.0, logk0_CO2);
        return k0_CO2;
    }

    /*!
     * \brief Returns the equilibrium constant for H2O, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003)
     * \param temperature the temperature [K]
     */
    template <class Evaluation>
    static Evaluation equilibriumConstantH2O_(const Evaluation& temperature)
    {
        Evaluation temperatureCelcius = temperature - 273.15;
        static const Scalar c[4] = { -2.209, 3.097e-2, -1.098e-4, 2.048e-7 };
        Evaluation logk0_H2O =
            c[0] + temperatureCelcius*(c[1] + temperatureCelcius*(c[2] + temperatureCelcius*c[3]));
        return Opm::pow(10.0, logk0_H2O);
    }

};

} // namespace BinaryCoeff
} // namespace Opm

#endif
