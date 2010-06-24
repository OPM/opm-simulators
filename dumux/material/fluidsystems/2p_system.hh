// $Id$
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
#ifndef DUMUX_2P_SYSTEM_HH
#define DUMUX_2P_SYSTEM_HH

#include "liquidphase.hh"
#include "gasphase.hh"

namespace Dumux
{

/*!
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class FluidSystem2P
{
    typedef FluidSystem2P<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;

public:

    enum
    {
        wPhaseIdx = 0,
        nPhaseIdx = 1
    };

    static const int numPhases = 2;

    FluidSystem2P()
    {
    }

    static void init()
    {
    }

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int phaseIdx)
    {
        switch (phaseIdx) {
        case wPhaseIdx: return WettingPhase::name();
        case nPhaseIdx: return NonwettingPhase::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << phaseIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int phaseIdx)
    {
        switch (phaseIdx) {
        case wPhaseIdx: return WettingPhase::molarMass();
        case nPhaseIdx: return NonwettingPhase::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << phaseIdx);
    }

    /*!
     * \brief Return the density of a component for a phase [kg/m^3].
     */
    static Scalar componentDensity(int phaseIdx,
                                   int compIdx,
                                   Scalar temperature,
                                   Scalar pressure)
    {
        if (phaseIdx != compIdx)
            return 0;
        switch (phaseIdx) {
        case wPhaseIdx:
            return WettingPhase::density(temperature, pressure);
        case nPhaseIdx:
            return NonwettingPhase::density(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }


    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &phaseState)
    {
        switch (phaseIdx) {
        case wPhaseIdx:
            return WettingPhase::density(temperature, pressure);
        case nPhaseIdx:
            return NonwettingPhase::density(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &phaseState)
    {
        switch (phaseIdx) {
        case wPhaseIdx:
            return WettingPhase::viscosity(temperature, pressure);
        case nPhaseIdx:
            return NonwettingPhase::viscosity(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                           Scalar temperature,
                           Scalar pressure,
                           const FluidState &phaseState)
    {
        switch (phaseIdx) {
        case wPhaseIdx:
            return WettingPhase::enthalpy(temperature, pressure);
        case nPhaseIdx:
            return NonwettingPhase::enthalpy(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        internal energy of the phase [J/kg].
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &phaseState)
    {
        switch (phaseIdx) {
        case wPhaseIdx:
            return WettingPhase::internalEnergy(temperature, pressure);
        case nPhaseIdx:
            return NonwettingPhase::internalEnergy(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

};
} // end namepace

#endif
