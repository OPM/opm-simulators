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
 * \brief A fluid system with water and air as phases and \f$H_2O\f$, \f$N_2\f$ and
 *        \f$O_2\f$ as components.
 */
#ifndef DUMUX_H2O_N2_O2_SYSTEM_HH
#define DUMUX_H2O_N2_O2_SYSTEM_HH

#include <dumux/new_material/components/h2o.hh>
#include <dumux/new_material/components/n2.hh>
#include <dumux/new_material/components/o2.hh>
#include <dumux/new_material/components/tabulatedcomponent.hh>

#include <dumux/new_material/binarycoefficients/h2o_n2.hh>
#include <dumux/new_material/binarycoefficients/h2o_o2.hh>
#include <dumux/new_material/binarycoefficients/n2_o2.hh>

namespace Dune
{

/*!
 * \brief A fluid system with water, molecular nitrogen and molecular
 *        oxygen as components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class H2O_N2_O2_System
{
    typedef H2O_N2_O2_System<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef Dune::IdealGas<Scalar> IdealGas;

    typedef Dune::H2O<Scalar>                           H2O_IAPWS;
    typedef Dune::TabulatedComponent<Scalar, H2O_IAPWS> H2O_Tabulated;

public:
    typedef H2O_Tabulated     H2O;
    //typedef H2O_IAPWS         H2O;
    typedef Dune::N2<Scalar>  N2;
    typedef Dune::O2<Scalar>  O2;

    static const int numComponents = 3;
    static const int numPhases = 2;
    
    // phase indices
    static const int lPhaseIdx = 0; // index of the liquid phase 
    static const int gPhaseIdx = 1; // index of the gas phase 
    
    // component indices
    static const int H2OIdx = 0; // component index of water
    static const int N2Idx = 1; // component index of nitrogen
    static const int O2Idx = 2; // component index of oxygen

    H2O_N2_O2_System()
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
        case O2Idx: return O2::name();
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
        case O2Idx: return O2::molarMass();
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
        case O2Idx: return O2::vaporPressure(temperature);
        };
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }
    
    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class Flash>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature, 
                               Scalar phasePressure,
                               const Flash &flash)
    { 
        switch (phaseIdx) {
        case lPhaseIdx: 
        {
            // See: Ochs 2008
            // \todo: proper citation
            Scalar rhoWater = H2O::liquidDensity(temperature, 
                                                 phasePressure);
            Scalar cWater = rhoWater/H2O::molarMass();
            return 
                (1 - 
                 flash.moleFrac(lPhaseIdx, N2Idx) -
                 flash.moleFrac(lPhaseIdx, O2Idx)) * rhoWater
                + 
                flash.moleFrac(lPhaseIdx, 1)*cWater*N2::molarMass() +
                flash.moleFrac(lPhaseIdx, 2)*cWater*O2::molarMass();
        }
        case gPhaseIdx:
        {
            // assume ideal gas
            Scalar avgMolarMass = 
                flash.moleFrac(gPhaseIdx, H2OIdx)*H2O::molarMass() + 
                flash.moleFrac(gPhaseIdx, N2Idx)*N2::molarMass() +
                flash.moleFrac(gPhaseIdx, O2Idx)*O2::molarMass();
            return IdealGas::density(avgMolarMass, temperature, phasePressure);
        };
        }
        DUNE_THROW(InvalidStateException, "Invalid phase index " << phaseIdx);

    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class Flash>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature, 
                                 Scalar phasePressure,
                                 const Flash &flash)
    { 
        if (phaseIdx == lPhaseIdx)
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            return H2O::liquidViscosity(temperature, phasePressure);
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
                H2O::gasViscosity(temperature, H2O::vaporPressure(temperature)),
                N2::gasViscosity(temperature, phasePressure),
                O2::gasViscosity(temperature, phasePressure)
            };
            // molar masses
            const Scalar M[numComponents] = {
                H2O::molarMass(),
                N2::molarMass(),
                O2::molarMass()
            };

            for (int i = 0; i < numComponents; ++i) {
                Scalar divisor = 0;
                for (int j = 0; j < numComponents; ++j) {
                    Scalar phiIJ = 1 + sqrt(mu[i]/mu[j] *
                                            pow(M[i]/M[j], 1/4.0));
                    phiIJ *= phiIJ;
                    phiIJ /= sqrt(8*(1 + M[i]/M[j]));
                    divisor += flash.moleFrac(phaseIdx, j)*phiIJ;
                }
                muResult += flash.moleFrac(phaseIdx, i)*mu[i] / divisor;
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
    template <class Flash>
    static Scalar dPg_dxl(int compIdx, 
                          Scalar temperature,
                          const Flash &flash)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return BinaryCoeff::H2O_N2::henry(temperature);
        case O2Idx: return BinaryCoeff::H2O_O2::henry(temperature);
        };
        DUNE_THROW(InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent of a component in a phase.
     */
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature, 
                            Scalar phasePressure)
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
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(temperature, phasePressure);
                case O2Idx: return BinaryCoeff::H2O_O2::liquidDiffCoeff(temperature, phasePressure);
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
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(temperature, phasePressure);
                case O2Idx: return BinaryCoeff::H2O_O2::gasDiffCoeff(temperature, phasePressure);
                }
               
            case N2Idx:
                // compJIdx == O2Idx
                return BinaryCoeff::N2_O2::gasDiffCoeff(temperature, phasePressure);
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
    template <class Flash>
    static Scalar enthalpy(int phaseIdx,
                           Scalar temperature, 
                           Scalar phasePressure,
                           const Flash &flash)
    { 
        if (phaseIdx == gPhaseIdx) {
            Scalar result = 0;
            result += H2O::gasEnthalpy(temperature, flash.partialPressure(0))*flash.massFrac(gPhaseIdx, 0);
            result += N2::gasEnthalpy(temperature, flash.partialPressure(1) )*flash.massFrac(gPhaseIdx, 1);
            result += O2::gasEnthalpy(temperature, flash.partialPressure(2) )*flash.massFrac(gPhaseIdx, 2);
            
            return result;
        }
        else {
            // TODO (?): solutes are not yet considered!
            return H2O::liquidEnthalpy(temperature, phasePressure);            
        }
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase's
     *        specific internal energy [J/kg].
     */
    template <class Flash>
    static Scalar internalEnergy(int phaseIdx,
                                 Scalar temperature, 
                                 Scalar phasePressure,
                                 const Flash &flash)
    { 
        if (phaseIdx == lPhaseIdx) 
            return enthalpy(phaseIdx, temperature, phasePressure, flash);
        else {
            return
                enthalpy(phaseIdx, temperature, phasePressure, flash)
                - phasePressure/density(phaseIdx, temperature, phasePressure, flash);
        }
        DUNE_THROW(InvalidStateException, 
                   "Invalid phase index: " << phaseIdx);
    }

};

} // end namepace

#endif
