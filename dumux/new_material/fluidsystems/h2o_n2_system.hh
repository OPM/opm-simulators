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
 * \brief A fluid system with water and gas as phases and \f$H_2O\f$ and \f$N_2\f$
 *        as components.
 */
#ifndef DUMUX_H2O_N2_SYSTEM_HH
#define DUMUX_H2O_N2_SYSTEM_HH

#include <dumux/new_material/idealgas.hh>
#include <dumux/new_material/components/n2.hh>
#include <dumux/new_material/components/h2o.hh>
#include <dumux/new_material/components/simpleh2o.hh>
#include <dumux/new_material/components/tabulatedcomponent.hh>

#include <dumux/common/properties.hh>

#include <dumux/new_material/binarycoefficients/h2o_n2.hh>

#define USE_SIMPLE_WATER 0

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
};

/*!
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class H2O_N2_System
{
    typedef H2O_N2_System<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef Dumux::IdealGas<Scalar> IdealGas;

    typedef Dumux::SimpleH2O<Scalar>                     SimpleH2O;
    typedef Dumux::H2O<Scalar>                           H2O_IAPWS;
    typedef Dumux::TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;

public:
#if ! USE_SIMPLE_WATER
    typedef H2O_Tabulated                             H2O;
    //typedef H2O_IAPWS                                 H2O;
#else
    typedef SimpleH2O                                 H2O;
#endif
    typedef Dumux::N2<Scalar>                          N2;

    static const int numComponents = 2;
    static const int numPhases = 2;

    static const int lPhaseIdx = 0; // index of the liquid phase 
    static const int gPhaseIdx = 1; // index of the gas phase 

    static const int wPhaseIdx = lPhaseIdx; // index of the wetting phase 
    static const int nPhaseIdx = gPhaseIdx; // index of the non-wetting phase 

    static const int H2OIdx = 0;
    static const int N2Idx = 1;
    
    static void init()
    {
#if ! USE_SIMPLE_WATER
        std::cout << "Initializing tables for the H2O fluid properties.\n";
        H2O_Tabulated::init(273.15, 623.15, 100,
                            -10,      20e6, 200);
#endif

#if 0
        // this used is to print the quanties of water into a CSV file
        // which can be inspected by matlab to check the
        // regularizations, etc.
        int nT = 300;
        Scalar minT = 273.15;
        Scalar maxT = 620.0;
        int nP = 600;
        Scalar minP = 0.0;
        Scalar maxP = 10e6;
        std::cerr << "'T' 'p' 'u_l' 'u_g' 'h_l' 'h_g' 'u_l,tab' 'u_g,tab' 'h_l,tab' 'h_g,tab'\n";
        for (int iT = 0; iT < nT; ++ iT) {
            Scalar T = Scalar(iT)/(nT - 1)*(maxT - minT) + minT;
            for (int iP = 0; iP < nP; ++ iP) {
                Scalar p = Scalar(iP)/(nP - 1)*(maxP - minP) + minP;
                std::cerr << T << " "
                          << p << " "
                          << H2O_IAPWS::liquidInternalEnergy(T, p) << " "
                          << H2O_IAPWS::gasInternalEnergy(T, p) << " "
                          << H2O_IAPWS::liquidEnthalpy(T, p) << " "
                          << H2O_IAPWS::gasEnthalpy(T, p) << " "
                          << H2O_Tabulated::liquidInternalEnergy(T, p) << " "
                          << H2O_Tabulated::gasInternalEnergy(T, p) << " "
                          << H2O_Tabulated::liquidEnthalpy(T, p) << " "
                          << H2O_Tabulated::gasEnthalpy(T, p) << "\n";
            }
            std::cerr << "\n";
        };
#endif
    }

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case N2Idx: return N2::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case N2Idx: return N2::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given the gas phase's composition, temperature and
     *        pressure, compute the partial presures of all components
     *        in [Pa] and assign it to the FluidState.
     *
     * This is required for models which cannot calculate the the
     * partial pressures of the components in the gas phase from the
     * degasPressure(). To use this method, the FluidState must have a
     * setPartialPressure(componentIndex, pressure) method.
     */
    template <class FluidState>
    static void computePartialPressures(Scalar temperature,
                                        Scalar pg,
                                        FluidState &fluidState)
    {
        Scalar X1 = fluidState.massFrac(gPhaseIdx, H2OIdx);

        // We use the newton method for this. For the initial value we
        // assume all components to be an ideal gas
        Scalar pH2O = 
            fluidState.moleFrac(gPhaseIdx, H2OIdx) * pg;
        Scalar eps = pg*1e-9;

        Scalar deltaP = pH2O;
        Valgrind::CheckDefined(pH2O);
        Valgrind::CheckDefined(deltaP);
        for (int i = 0; i < 5 && std::abs(deltaP/pg) > 1e-9; ++i) {
            Scalar f = 
                H2O::gasDensity(temperature, pH2O)*(1 - 1/X1) +
                N2::gasDensity(temperature, pg - pH2O);
            
            Scalar df_dp;
            df_dp  =
                H2O::gasDensity(temperature, pH2O + eps)*(1 - 1/X1) +
                N2::gasDensity(temperature, pg - (pH2O + eps));
            df_dp -= 
                H2O::gasDensity(temperature, pH2O - eps)*(1 - 1/X1) +
                N2::gasDensity(temperature, pg - (pH2O - eps));
            df_dp /= 
                2*eps;
            
            deltaP = - f/df_dp;
            
            pH2O += deltaP;
            Valgrind::CheckDefined(pH2O);
            Valgrind::CheckDefined(deltaP);
        }
        
        fluidState.setPartialPressure(H2OIdx, pH2O);
        fluidState.setPartialPressure(N2Idx, pg - pH2O);
    };
   
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    { 
        switch (phaseIdx) {
        case lPhaseIdx: 
        {
            // See: Ochs 2008
            // \todo: proper citation
            Scalar rhoWater = H2O::liquidDensity(temperature, pressure);
            Scalar cWater = rhoWater/H2O::molarMass();
            return 
                fluidState.massFrac(lPhaseIdx, H2OIdx)*rhoWater
                + 
                fluidState.moleFrac(lPhaseIdx, N2Idx)*cWater*N2::molarMass();
        }
        case gPhaseIdx:
        {
            Scalar pH2O = fluidState.partialPressure(H2OIdx);
            Scalar pN2 = fluidState.partialPressure(N2Idx);

            // using the two partial pressures we can calculate the
            // density of the phase. This assumes that Dalton's law is
            // valid
            return 
                H2O::gasDensity(temperature, pH2O) +
                N2::gasDensity(temperature, pN2);
        };
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its viscosity.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature, 
                                 Scalar pressure,
                                 const FluidState &fluidState)
    { 
        if (phaseIdx == lPhaseIdx) {
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            return H2O::liquidViscosity(temperature, 
                                        pressure);
        }
        else {
            return N2::gasViscosity(temperature, 
                                    pressure);
            /* Wilke method. See:
             *
             * S.O.Ochs: "Development of a multiphase multicomponent
             * model for PEMFC - Technical report: IRTG-NUPUS",
             * University of Stuttgart, 2008
             *
             * and:
             *
             * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
             * edition, McGraw-Hill, 1987, 407-410
             */
            Scalar muResult = 0; 
            const Scalar mu[numComponents] = {
                H2O::gasViscosity(temperature,
                                  H2O::vaporPressure(temperature)),
                N2::gasViscosity(temperature, 
                                 pressure)
            };
            // molar masses
            const Scalar M[numComponents] = {
                H2O::molarMass(),
                N2::molarMass()
            };

            for (int i = 0; i < numComponents; ++i) {
                Scalar divisor = 0;
                for (int j = 0; j < numComponents; ++j) {
                    Scalar phiIJ = 1 + sqrt(mu[i]/mu[j] *
                                            pow(M[i]/M[j], 1/4.0));
                    phiIJ *= phiIJ;
                    phiIJ /= sqrt(8*(1 + M[i]/M[j]));
                    divisor += fluidState.moleFrac(phaseIdx, j)*phiIJ;
                }
                muResult += fluidState.moleFrac(phaseIdx, i)*mu[i] / divisor;
            }

            return muResult;
        }
    } 

    /*!
     * \brief Return the pressure which a component degases from the
     *        liquid scaled to \f$100\%\f$ of the component.
     *
     * For solutions with only traces in a solvent this boils down to
     * the inverse Henry constant for the solutes and the partial
     * pressure for the solvent.
     */
    static Scalar degasPressure(int compIdx, 
                                Scalar temperature,
                                Scalar pg)
    {        
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return BinaryCoeff::H2O_N2::henry(temperature);
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }
    
    /*!
     * \brief Given a component's pressure and temperature, return its
     *        density in a phase [kg/m^3].
     */
    static Scalar componentDensity(int phaseIdx,
                                   int compIdx,
                                   Scalar temperature, 
                                   Scalar pressure)
    {
        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::liquidDensity(temperature, pressure);
            else if (compIdx == N2Idx)
                return N2::liquidDensity(temperature, pressure);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        else if (phaseIdx == gPhaseIdx) {
            if (compIdx == H2OIdx) {
                return H2O::gasDensity(temperature, pressure);
            }
            else if (compIdx == N2Idx)
                return N2::gasDensity(temperature, pressure);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    };

    /*!
     * \brief Given a component's density and temperature, return the
     *        corresponding pressure in a phase [Pa].
     */
    static Scalar componentPressure(int phaseIdx,
                                    int compIdx,
                                    Scalar temperature, 
                                    Scalar density)
    {
        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::liquidPressure(temperature, density);
            else if (compIdx == N2Idx)
                return N2::liquidPressure(temperature, density);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        else if (phaseIdx == gPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::gasPressure(temperature, density);
            else if (compIdx == N2Idx)
                return N2::gasPressure(temperature, density);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficent for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    { 
        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);
        
#ifndef NDEBUG
        if (compIIdx == compJIdx || 
            phaseIdx > numPhases - 1 ||
            compJIdx > numComponents - 1)
        {
            DUNE_THROW(Dune::InvalidStateException, 
                       "Binary diffusion coefficient of components " 
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
#endif

        
        switch (phaseIdx) {
        case lPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(temperature, 
                                                                        pressure);
                }
            default:
                DUNE_THROW(Dune::InvalidStateException, 
                           "Binary diffusion coefficients of trace "
                           "substances in liquid phase is undefined!\n");
            }
        case gPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(temperature,
                                                                     pressure);
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException, 
                   "Binary diffusion coefficient of components " 
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &fluidState)
    { 
        if (phaseIdx == lPhaseIdx)  {
            Scalar cN2 = fluidState.concentration(lPhaseIdx, N2Idx);
            Scalar pN2 = N2::gasPressure(temperature, cN2*N2::molarMass());

            // TODO: correct way to deal with the solutes??? 
            return 
                fluidState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidEnthalpy(temperature, pressure)
                +
                fluidState.massFrac(lPhaseIdx, N2Idx)*
                N2::gasEnthalpy(temperature, pN2);
        }
        else {
            Scalar pWater = fluidState.partialPressure(H2OIdx);
            Scalar pN2 = fluidState.partialPressure(N2Idx);

            Scalar result = 0;
            result += 
                H2O::gasEnthalpy(temperature, pWater) *
                fluidState.massFrac(gPhaseIdx, H2OIdx);
            result += 
                N2::gasEnthalpy(temperature, pN2) *
                fluidState.massFrac(gPhaseIdx, N2Idx);
            
            return result;
        }
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific internal energy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                      Scalar temperature,
                                      Scalar pressure,
                                      const FluidState &fluidState)
    { 
        if (phaseIdx == lPhaseIdx)  {
            Scalar cN2 = fluidState.concentration(lPhaseIdx, N2Idx);
            Scalar pN2 = N2::gasPressure(temperature, cN2*N2::molarMass());

            // TODO: correct way to deal with the solutes??? 
            return 
                fluidState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidInternalEnergy(temperature, pressure)
                +
                fluidState.massFrac(lPhaseIdx, N2Idx)*
                N2::gasInternalEnergy(temperature, pN2);
        }
        else {
            Scalar pWater = fluidState.partialPressure(H2OIdx);
            Scalar pN2 = fluidState.partialPressure(N2Idx);
            
            Scalar result = 0;
            result += 
                H2O::gasInternalEnergy(temperature, pWater)*
                fluidState.massFrac(gPhaseIdx, H2OIdx);
            result += 
                N2::gasInternalEnergy(temperature, pN2)*
                fluidState.massFrac(gPhaseIdx, N2Idx);
            
            return result;
        }
    }
};

} // end namepace

#endif
