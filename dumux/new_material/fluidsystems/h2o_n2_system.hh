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

#include <dumux/new_material/binarycoefficients/h2o_n2.hh>

namespace Dune
{

/*!
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class H2O_N2_System
{
    typedef H2O_N2_System<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef Dune::IdealGas<Scalar> IdealGas;

    typedef Dune::H2O<Scalar>                           H2O_IAPWS;
    typedef Dune::TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;

public:
    typedef H2O_Tabulated                             H2O;
    //typedef H2O_IAPWS                                 H2O;
    //typedef Dune::SimpleH2O<Scalar>                   H2O;
    typedef Dune::N2<Scalar>                          N2;

    static const int numComponents = 2;
    static const int numPhases = 2;

    static const int lPhaseIdx = 0; // index of the liquid phase 
    static const int gPhaseIdx = 1; // index of the gas phase 

    static const int H2OIdx = 0;
    static const int N2Idx = 1;
    
    H2O_N2_System()
    {
    }

    static void init()
    {
        std::cout << "Initializing tables for the H2O fluid properties.\n";
        H2O_Tabulated::init(273.15, 623.15, 100,
                            -10,      20e6, 200);
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
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case N2Idx: return N2::molarMass();
        };
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the vapor pressure of a component in [Pa].
     */
    static Scalar vaporPressure(int compIdx, 
                                Scalar temperature)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return N2::vaporPressure(temperature);
        };
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }
    
    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class PhaseState>
    static Scalar phaseDensity(int phaseIdx,
                               const PhaseState &phaseState)
    { 
        switch (phaseIdx) {
        case lPhaseIdx: 
        {
            Scalar pVap = 0.9*H2O::vaporPressure(phaseState.temperature());
            Scalar pressure = phaseState.phasePressure(lPhaseIdx);
            if (pressure < pVap)
                pressure = pVap;

            // See: Ochs 2008
            // \todo: proper citation
            Scalar rhoWater = H2O::liquidDensity(phaseState.temperature(), 
                                                 pressure);
            Scalar cWater = rhoWater/H2O::molarMass();
            return 
                phaseState.moleFrac(lPhaseIdx, H2OIdx)*rhoWater
                + 
                phaseState.moleFrac(lPhaseIdx, N2Idx)*cWater*N2::molarMass();
        }
        case gPhaseIdx:
        {
            // assume ideal gas
            Scalar avgMolarMass = 
                phaseState.moleFrac(gPhaseIdx, N2Idx)*N2::molarMass();
            Scalar pWater = phaseState.partialPressure(H2OIdx);

            Scalar pVap = 1.1*H2O::vaporPressure(phaseState.temperature());
            if (pWater > pVap)
                pWater = pVap;
            
            Scalar rhoWater = H2O::gasDensity(phaseState.temperature(), pWater);
            return 
                rhoWater +
                IdealGas::density(avgMolarMass, 
                                  phaseState.temperature(), 
                                  phaseState.phasePressure(gPhaseIdx) - pWater);
        };
        }
        DUNE_THROW(InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class PhaseState>
    static Scalar phaseViscosity(int phaseIdx,
                                 const PhaseState &phaseState)
    { 
        if (phaseIdx == lPhaseIdx) {
            Scalar pVap = 0.9*H2O::vaporPressure(phaseState.temperature());
            Scalar pressure = phaseState.phasePressure(lPhaseIdx);
            if (pressure < pVap)
                pressure = pVap;
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            return H2O::liquidViscosity(phaseState.temperature(), 
                                        pressure);
        }
        else {
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
                H2O::gasViscosity(phaseState.temperature(), 
                                  H2O::vaporPressure(phaseState.temperature())),
                N2::gasViscosity(phaseState.temperature(), 
                                 phaseState.phasePressure(gPhaseIdx))
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
                    divisor += phaseState.moleFrac(phaseIdx, j)*phiIJ;
                }
                muResult += phaseState.moleFrac(phaseIdx, i)*mu[i] / divisor;
            }

            return muResult;
        }
    } 

    /*!
     * \brief Returns the derivative of the equilibrium partial
     *        pressure \f$\partial p^\kappa_g / \partial x^\kappa_l\$
     *        to the mole fraction of a component in the liquid phase.
     *
     * For solutions with only traces in a solvent this boils down to
     * the inverse Henry constant for the solutes and the partial
     * pressure for the solvent.
     */
    template <class PhaseState>
    static Scalar dPg_dxl(int compIdx, 
                          const PhaseState &phaseState)
    {        
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(phaseState.temperature());
        case N2Idx: return BinaryCoeff::H2O_N2::henry(phaseState.temperature());
        };
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent of a component in a phase.
     */
    template <class PhaseState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            const PhaseState &phaseState)
    { 
        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);
        
#ifndef NDEBUG
        if (compIIdx == compJIdx || 
            phaseIdx > numPhases - 1 ||
            compJIdx > numComponents - 1)
        {
            DUNE_THROW(InvalidStateException, 
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
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(phaseState.temperature(), 
                                                                        phaseState.phasePressure(phaseIdx));
                }
            default:
                DUNE_THROW(InvalidStateException, 
                           "Binary diffusion coefficients of trace "
                           "substances in liquid phase is undefined!\n");
            }
        case gPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(phaseState.temperature(),
                                                                     phaseState.phasePressure(phaseIdx));
                }
            }
        }

        DUNE_THROW(InvalidStateException, 
                   "Binary diffusion coefficient of components " 
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    template <class PhaseState>
    static Scalar enthalpy(int phaseIdx,
                           const PhaseState &phaseState)
    { 
        Scalar temperature = phaseState.temperature();
        if (phaseIdx == lPhaseIdx)  {
            Scalar pVap = 1.1*H2O::vaporPressure(temperature);
            Scalar pWater = phaseState.phasePressure(lPhaseIdx);
            if (pWater < pVap)
                pWater = pVap;

            Scalar cN2 = phaseState.concentration(lPhaseIdx, N2Idx);
            Scalar pN2 = IdealGas::pressure(temperature, cN2);

            // TODO: correct way to deal with the solutes??? 
            return 
                phaseState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidEnthalpy(temperature, pWater)
                +
                phaseState.massFrac(lPhaseIdx, N2Idx)*
                N2::gasEnthalpy(temperature, pN2);
        }
        else {
            Scalar pVap = 0.9*H2O::vaporPressure(temperature);
            Scalar pWater = phaseState.partialPressure(H2OIdx);
            if (pWater > pVap)
                pWater = pVap;

            Scalar pN2 = phaseState.partialPressure(N2Idx);

            Scalar result = 0;
            result += 
                H2O::gasEnthalpy(temperature, pWater) *
                phaseState.massFrac(gPhaseIdx, H2OIdx);
            result += 
                N2::gasEnthalpy(temperature, pN2) *
                phaseState.massFrac(gPhaseIdx, N2Idx);
            
            return result;
        }
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase's
     *        internal energy [J/kg].
     */
    template <class PhaseState>
    static Scalar internalEnergy(int phaseIdx,
                                 const PhaseState &phaseState)
    { 
        Scalar temperature = phaseState.temperature();
        if (phaseIdx == lPhaseIdx)  {
            Scalar pVap = 1.1*H2O::vaporPressure(temperature);
            Scalar pWater = phaseState.phasePressure(lPhaseIdx);
            if (pWater < pVap)
                pWater = pVap;

            Scalar cN2 = phaseState.concentration(lPhaseIdx, N2Idx);
            Scalar pN2 = IdealGas::pressure(temperature, cN2);

            // TODO: correct way to deal with the solutes??? 
            return 
                phaseState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidInternalEnergy(temperature, pWater)
                +
                phaseState.massFrac(lPhaseIdx, N2Idx)*
                N2::gasInternalEnergy(temperature, pN2);
        }
        else {
            Scalar pVap = 0.9*H2O::vaporPressure(temperature);
            Scalar pWater = phaseState.partialPressure(H2OIdx);
            if (pWater > pVap)
                pWater = pVap;

            Scalar pN2 = phaseState.partialPressure(N2Idx);
            
            Scalar result = 0;
            result += 
                H2O::gasInternalEnergy(temperature, pWater)*
                phaseState.massFrac(gPhaseIdx, H2OIdx);
            result += 
                N2::gasInternalEnergy(temperature, pN2)*
                phaseState.massFrac(gPhaseIdx, N2Idx);
            
            return result;
        }
    }
};

} // end namepace

#endif
