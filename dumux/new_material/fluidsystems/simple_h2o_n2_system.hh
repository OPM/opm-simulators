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
#ifndef DUMUX_SIMPLE_H2O_N2_SYSTEM_HH
#define DUMUX_SIMPLE_H2O_N2_SYSTEM_HH

#include <dumux/new_material/idealgas.hh>
#include <dumux/new_material/components/n2.hh>
#include <dumux/new_material/components/h2o.hh>
#include <dumux/new_material/components/simpleh2o.hh>
#include <dumux/new_material/components/tabulatedcomponent.hh>

#include <dumux/new_material/binarycoefficients/h2o_n2.hh>

namespace Dumux
{

/*!
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class Simple_H2O_N2_System
{
    typedef Simple_H2O_N2_System<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    typedef Dumux::SimpleH2O<Scalar>                   H2O;
    typedef Dumux::N2<Scalar>                          N2;

    static const int numComponents = 2;
    static const int numPhases = 2;

    static const int lPhaseIdx = 0; // index of the liquid phase 
    static const int gPhaseIdx = 1; // index of the gas phase 

    static const int H2OIdx = 0;
    static const int N2Idx = 1;
    
    static void init()
    {
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
     * \brief Return the molar mass of a component in [kg/mol].
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
     * \brief Return the vapor pressure of a component in [Pa].
     */
    static Scalar vaporPressure(int compIdx, 
                                Scalar temperature)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return N2::vaporPressure(temperature);
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }
    
    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               const FluidState &phaseState)
    { 
        switch (phaseIdx) {
        case lPhaseIdx: 
        {
            Scalar pressure = phaseState.phasePressure(lPhaseIdx);
            return H2O::liquidDensity(phaseState.temperature(),
                                      pressure);
        }
        case gPhaseIdx:
        {
            // assume ideal gas
            Scalar avgMolarMass = 
                phaseState.moleFrac(gPhaseIdx, H2OIdx)*H2O::molarMass() +
                phaseState.moleFrac(gPhaseIdx, N2Idx)*N2::molarMass();

            return 
                IdealGas::density(avgMolarMass, 
                                  phaseState.temperature(), 
                                  phaseState.phasePressure(gPhaseIdx));
        };
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 const FluidState &phaseState)
    { 
        if (phaseIdx == lPhaseIdx) {
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(phaseState.temperature(), 
                                        phaseState.phasePressure(lPhaseIdx));
        }
        else {
            // assume pure nitrogen for the gas phase
            return N2::gasViscosity(phaseState.temperature(), 
                                    phaseState.phasePressure(gPhaseIdx));
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
    template <class FluidState>
    static Scalar dPg_dxl(int compIdx, 
                          const FluidState &phaseState)
    {        
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(phaseState.temperature());
        case N2Idx: return BinaryCoeff::H2O_N2::henry(phaseState.temperature());
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent of a component in a phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            const FluidState &phaseState)
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
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(phaseState.temperature(), 
                                                                        phaseState.phasePressure(phaseIdx));
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
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(phaseState.temperature(),
                                                                     phaseState.phasePressure(phaseIdx));
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException, 
                   "Binary diffusion coefficient of components " 
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar enthalpy(int phaseIdx,
                           const FluidState &phaseState)
    { 
        if (phaseIdx == lPhaseIdx) {
            Scalar temperature = phaseState.temperature();
            Scalar pressure = phaseState.phasePressure(lPhaseIdx);
            
            return 
                //phaseState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidEnthalpy(temperature, pressure);
        }
        else {
            Scalar pWater = phaseState.partialPressure(H2OIdx);

            Scalar result = 0;
            result += 
                H2O::gasEnthalpy(phaseState.temperature(), pWater) *
                phaseState.massFrac(gPhaseIdx, H2OIdx);
            result +=  
                N2::gasEnthalpy(phaseState.temperature(), 
                                phaseState.partialPressure(N2Idx)) *
                phaseState.massFrac(gPhaseIdx, N2Idx);
            
            return result;
        }
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase's
     *        internal energy [J/kg].
     */
    template <class FluidState>
    static Scalar internalEnergy(int phaseIdx,
                                 const FluidState &phaseState)
    { 
        if (phaseIdx == lPhaseIdx) 
            return enthalpy(phaseIdx, phaseState);
        else {
            return
                enthalpy(phaseIdx, phaseState)
                - IdealGas::R*phaseState.temperature(); // = pressure * spec. volume for an ideal gas
        }
    }
};

} // end namepace

#endif
