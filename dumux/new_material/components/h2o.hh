/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_H2O_HH
#define DUMUX_H2O_HH

#include <dumux/new_material/idealgas.hh>
#include <dune/common/exceptions.hh>
#include <dumux/exceptions.hh>

#include "component.hh"

#include <cmath>
#include <iostream>

namespace Dune
{
/*!
 * \brief Properties of pure water \f$H_2O\f$.
 *
 * See: 
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class H2O : public Component<Scalar, H2O<Scalar> >
{
    typedef Component<Scalar, H2O<Scalar> > ParentType;
    typedef Dune::IdealGas<Scalar> IdealGas;

    static const Scalar R = 461.526;  // specific gas constant of water
public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "H2O"; } 

    /*!
     * \brief The mass in [kg] of one mole of water.
     */
    static Scalar molarMass()
    { return 18e-3; } 

    /*!
     * \brief Returns the critical temperature [K] of water
     */
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of water
     */
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K]at water's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at water's triple point.
     */
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [N/m^2] of pure water
     *        at a given temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar vaporPressure(Scalar T)
    { 
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account
        
        static const Scalar n[10] = {
            0.11670521452767e4,  -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5,  -0.32325550322333e7,  0.14915108613530e2,
            -0.48232657361591e4,  0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };
         
        Scalar sigma = T + n[8]/(T - n[9]);
        
        Scalar A = (sigma + n[0])*sigma + n[1];
        Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        Scalar tmp = 2*C/(std::sqrt(B*B - 4*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Specific enthalpy of water steam [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar gasEnthalpy(Scalar temperature, 
                                    Scalar pressure)
    {
        // gas is only present if the vapor pressure is larger than the
        // partial pressure
        assert(vaporPressure(temperature) >= pressure);
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return 540.0 * gamma_tau_region2(temperature, pressure) * R;    /* J/kg */
    }

    /*!
     * \brief Specific enthalpy of liquid water [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        // liquid is only present if the vapor pressure is smaller than the
        // partial pressure
        assert(vaporPressure(temperature) <= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return 1386.0 * gamma_tau_region1(temperature, pressure) * R;    /* J/kg */
    }

    /*!
     * \brief Specific internal energy of liquid water [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        // liquid is only present if the vapor pressure is smaller than the
        // partial pressure
        assert(vaporPressure(temperature) <= pressure);
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar tau = 1386.0 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 16.53e6;    /* reduced pressure */
        return
            R * temperature *
            ( tau*gamma_tau_region1(temperature, pressure) - 
              pi*gamma_pi_region1(temperature, pressure));
    }

    /*!
     * \brief Specific internal energy of steam and water vapor [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
    */
    static Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        // gas is only present if the vapor pressure is larger than the
        // partial pressure
        assert(vaporPressure(temperature) >= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar tau = 540.0/temperature;
        Scalar pi = pressure/1.0e6;
        return
            R * temperature *
            ( tau*gamma_tau_region2(temperature, pressure) - 
              pi*gamma_pi_region2(temperature, pressure));
    }

    /*!
     * \brief The density of steam at a given pressure and temperature [kg/m^3].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // gas is only present if the vapor pressure is larger than the
        // partial pressure
        assert(vaporPressure(temperature) >= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Density of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar specificVolume = gamma_pi_region2(temperature, pressure) * R*temperature / 1e6;
        return 1/specificVolume;
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature [kg/m^3].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        // liquid is only present if the vapor pressure is smaller than the
        // partial pressure
        assert(vaporPressure(temperature) <= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Density of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar specificVolume = gamma_pi_region1(temperature, pressure) * R*temperature / 16.53e6;
        return 1/specificVolume;
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of steam.
     *
     * This method is only valid if pressure is below or at the vapour
     * pressure of water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        // gas is only present if the vapor pressure is larger than the
        // partial pressure
        assert(vaporPressure(temperature) >= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return viscosityIAPWS_(temperature, H2O::gasDensity(temperature, pressure));
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        // liquid is only present if the vapor pressure is smaller than the
        // partial pressure
        assert(vaporPressure(temperature) <= pressure);      
        if (temperature > 623.15 || pressure > 100e6)
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return viscosityIAPWS_(temperature, H2O::liquidDensity(temperature, pressure));
    };

private:
    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar viscosityIAPWS_(Scalar temperature, Scalar rho)
    {
        Scalar muBar;
        Scalar rhoBar = rho/322.0;
        Scalar TBar = temperature/H2O::criticalTemperature();

        // muBar = muBar_1
        const Scalar Hij[6][7] = {
            { 5.20094e-1, 2.22531e-1,-2.81378e-1, 1.61913e-1,-3.25372e-2,          0,          0 },
            { 8.50895e-2, 9.99115e-1,-9.06851e-1, 2.57399e-1,          0,          0,          0 },
            {-1.08374   , 1.88797   ,-7.72479e-1,          0,          0,          0,          0 },
            {-2.89555e-1, 1.26613   ,-4.89837e-1,          0, 6.98452e-2,          0,-4.35673e-3 },
            {          0,          0,-2.57040e-1,          0,          0, 8.72102e-3,          0 },
            {          0, 1.20573e-1,          0,          0,          0,          0,-5.93264e-4 }
        };

        Scalar tmp, tmp2, tmp3 = 1;
        muBar = 0;
        for (int i = 0; i <= 5; ++i) {
            tmp = 0;
            tmp2 = 1;
            for (int j = 0; j <= 6; ++j) {
                tmp += Hij[i][j]*tmp2;
                tmp2 *= (rhoBar - 1);
            };
            muBar += tmp3 * tmp;
            tmp3 *= 1.0/TBar - 1;
        };
        muBar *= rhoBar;
        muBar = std::exp(muBar);

        // muBar *= muBar_0
        muBar  *= 100*std::sqrt(TBar);
        const Scalar H[4] = {
            1.67752, 2.20462, 0.6366564, -0.241605
        };
        
        tmp = 0, tmp2 = 1;
        for (int i = 0; i < 4; ++i) {
            tmp += H[i]/tmp2;
            tmp2 *= TBar;
        };
        muBar /= tmp;

        return 1e-6*muBar;
    }

    static Scalar n_region1(int i)
    {
        static const Scalar n[34] =  {
            0.14632971213167,    -0.84548187169114,    -0.37563603672040e1, 
            0.33855169168385e1,  -0.95791963387872,     0.15772038513228,
           -0.16616417199501e-1,  0.81214629983568e-3,  0.28319080123804e-3,
           -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,
           -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,
           -0.30001780793026e-3,  0.47661393906987e-4, -0.44141845330846e-5,
           -0.72694996297594e-15,-0.31679644845054e-4, -0.28270797985312e-5,
           -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6, 
           -0.14341729937924e-12,-0.40516996860117e-6, -0.12734301741641e-8,
           -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,
            0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,
           -0.93537087292458e-25
        };
        return n[i];
    }

    static short int I_region1(int i)
    {
        static const short int I[34] = {
            0,  0,  0,
            0,  0,  0,
            0,  0,  1,
            1,  1,  1,
            1,  1,  2,
            2,  2,  2,
            2,  3,  3,
            3,  4,  4,
            4,  5,  8,
            8,  21, 23,
            29, 30, 31,
            32
        };
        return I[i];
    }

    static short int J_region1(int i)
    {
        static const short int J[34] = {
             -2,  -1,    0,
              1,   2,    3,
              4,   5,   -9,
             -7,  -1,    0,
              1,   3,   -3,
              0,   1,    3,
             17,  -4,    0,
              6,  -5,   -2,
             10,  -8,  -11,
             -6,  -29, -31,
            -38,  -39, -40,
            -41
        };
        return J[i];
    }

    /*!
     * The gibbs free energy for IAPWS region 1 (i.e. liquid)
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma_region1(Scalar temperature, Scalar pressure)
    {
        Scalar tau = 1386.0 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 16.53e6;    /* reduced pressure */
        
        Scalar result = 0;
        for (int i = 0; i < 34; ++i) {
            result += n_region1(i)*pow(7.1 - pi, I_region1(i))*pow(tau - 1.222, J_region1(i));
        }

        return result;
    }


    /*!
     * The partial derivative of the gibbs free energy to the
     * normalized temperature for IAPWS region 1 (i.e. liquid)
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma_tau_region1(Scalar temperature, Scalar pressure)
    {
        Scalar tau = 1386.0 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 16.53e6;    /* reduced pressure */
        
        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result += 
                n_region1(i) *
                std::pow(7.1 - pi, I_region1(i)) *
                std::pow(tau - 1.222,  J_region1(i)-1) *
                J_region1(i);
        }

        return result;
    }

    /*!
     * The partial derivative of the gibbs free energy to the
     * normalized pressure for IAPWS region 1 (i.e. liquid)
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma_pi_region1(Scalar temperature, Scalar pressure)
    {
        Scalar tau = 1386.0 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 16.53e6;    /* reduced pressure */
        
        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result += 
                -n_region1(i) *
                I_region1(i) *
                std::pow(7.1 - pi, I_region1(i) - 1) *
                std::pow(tau - 1.222,  J_region1(i));
        }

        return result;
    }


    static Scalar n_g_region2(int i)
    {
        static const Scalar n[9] =  {
            -0.96927686500217e1, 0.10086655968018e2, -0.56087911283020e-2,
            0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1,
            -0.43839511319450e1, -0.28408632460772, 0.21268463753307e-1
        };
        return n[i];
    }

    static Scalar n_r_region2(int i)
    {
        static const Scalar n[43] =  {
            -0.17731742473213e-2, -0.17834862292358e-1,  -0.45996013696365e-1,
            -0.57581259083432e-1, -0.50325278727930e-1,  -0.33032641670203e-4,
            -0.18948987516315e-3, -0.39392777243355e-2,  -0.43797295650573e-1,
            -0.26674547914087e-4,  0.20481737692309e-7,   0.43870667284435e-6,
            -0.32277677238570e-4, -0.15033924542148e-2,  -0.40668253562649e-1,
            -0.78847309559367e-9,  0.12790717852285e-7,   0.48225372718507e-6,
             0.22922076337661e-5, -0.16714766451061e-10, -0.21171472321355e-2,
            -0.23895741934104e2,  -0.59059564324270e-17, -0.12621808899101e-5,
            -0.38946842435739e-1,  0.11256211360459e-10, -0.82311340897998e1,
             0.19809712802088e-7,  0.10406965210174e-18, -0.10234747095929e-12,
            -0.10018179379511e-8, -0.80882908646985e-10,  0.10693031879409,
            -0.33662250574171,     0.89185845355421e-24,  0.30629316876232e-12,
            -0.42002467698208e-5, -0.59056029685639e-25,  0.37826947613457e-5,
            -0.12768608934681e-14, 0.73087610595061e-28,  0.55414715350778e-16,
            -0.94369707241210e-6
        };
        return n[i];
    }

    static Scalar I_r_region2(int i)
    {
        static const short int I[43] = {
            1, 1, 1,
            1, 1, 2,
            2, 2, 2,
            2, 3, 3,
            3, 3, 3,
            4, 4, 4,
            5, 6, 6,
            6, 7, 7,
            7, 8, 8,
            9, 10, 10,
            10, 16, 16,
            18, 20, 20,
            20, 21, 22,
            23, 24, 24,
            24
        };
        return I[i];
    }

    static Scalar J_g_region2(int i)
    {
        static const short int J[9] = {
            0,  1, -5,
            -4, -3, -2,
            -1,  2,  3
        };
        return J[i];
    }

    static Scalar J_r_region2(int i)
    {
        static const short int J[43] = {
            0, 1, 2,
            3, 6, 1,
            2, 4, 7,
            36, 0, 1,
            3, 6, 35,
            1, 2, 3,
            7, 3, 16,
            35, 0, 11,
            25, 8, 36,
            13, 4, 10,
            14, 29, 50,
            57, 20, 35,
            48, 21, 53,
            39, 26, 40,
            58
        };
        return J[i];
    }

    /*!
     * The partial derivative of the gibbs free energy to the
     * normalized temperature for IAPWS region 2 (i.e. sub-critical
     * steam)
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma_tau_region2(Scalar temperature, Scalar pressure)
    {
        Scalar tau = 540.0 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 1e6;    /* reduced pressure */
        
        // ideal gas part
        Scalar result = 0;
        for (int i = 0; i < 9; i++) {
            Scalar n_i = n_g_region2(i);
            Scalar J_i = J_g_region2(i);
            result += 
                n_i *
                J_i *
                std::pow(tau,  J_i - 1);
        }

        // residual part
        for (int i = 0; i < 43; i++) {
            Scalar n_i = n_r_region2(i);
            Scalar I_i = I_r_region2(i);
            Scalar J_i = J_r_region2(i);
            result += 
                n_i *
                std::pow(pi,  I_i) *
                J_i *
                std::pow(tau - 0.5,  J_i - 1);
        }

        return result;
    }

    /*!
     * The partial derivative of the gibbs free energy to the
     * normalized pressure for IAPWS region 2 (i.e. sub-critical
     * steam)
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma_pi_region2(Scalar temperature, Scalar pressure)
    {
        Scalar tau = 540 / temperature;   /* reduced temperature */
        Scalar pi = pressure / 1e6;    /* reduced pressure */
        
        // ideal gas part
        Scalar result = 1/pi;
        
        // residual part
        for (int i = 0; i < 43; i++) {
            Scalar n_i = n_r_region2(i);
            Scalar I_i = I_r_region2(i);
            Scalar J_i = J_r_region2(i);
            result += 
                n_i *
                I_i *
                std::pow(pi, I_i - 1) *
                std::pow(tau - 0.5,  J_i);
        }
        
        return result;
    }
};

} // end namepace

#endif
