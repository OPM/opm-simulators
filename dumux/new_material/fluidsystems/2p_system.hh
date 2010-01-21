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

#include <dumux/new_decoupled/2p/2pproperties.hh>

namespace Dune
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
        DUNE_THROW(InvalidStateException, "Invalid component index " << phaseIdx);
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
        DUNE_THROW(InvalidStateException, "Invalid component index " << phaseIdx);
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
        case wPhaseIdx:
        {
            return WettingPhase::density(phaseState.temperature(), phaseState.phasePressure(phaseIdx));
        }
        case nPhaseIdx:
        {
            return NonwettingPhase::density( phaseState.temperature(), phaseState.phasePressure(phaseIdx));
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
        switch (phaseIdx) {
        case wPhaseIdx:
        {
            return WettingPhase::viscosity(phaseState.temperature(), phaseState.phasePressure(phaseIdx));
        }
        case nPhaseIdx:
        {
            return NonwettingPhase::viscosity(phaseState.temperature(), phaseState.phasePressure(phaseIdx));
        };
        }
        DUNE_THROW(InvalidStateException, "Invalid phase index " << phaseIdx);
    } 

    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    template <class PhaseState>
    static Scalar enthalpy(int phaseIdx,
                           const PhaseState &phaseState)
    { 
        switch (phaseIdx) {
        case wPhaseIdx:
        {
            return WettingPhase::enthalpy(phaseState.temperature(), phaseState.phasePressure(phaseIdx));
        }
        case nPhaseIdx:
        {
            return NonwettingPhase::enthalpy(phaseState.temperature(), phaseState.phasePressure(phaseIdx));
        };
        }
        DUNE_THROW(InvalidStateException, "Invalid phase index " << phaseIdx);
    }

};

template <class Scalar, class Component> class GasPhase
{
public:
    /*!
     * \brief A human readable name for the compoent.
     */
    static const char *name()
    { return Component::name(); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static Scalar molarMass()
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static Scalar criticalTemperature()
    {  DUNE_THROW(NotImplemented, "Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static Scalar criticalPressure()
    {  DUNE_THROW(NotImplemented, "Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static Scalar tripleTemperature()
    {  DUNE_THROW(NotImplemented, "Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static Scalar triplePressure()
    {  DUNE_THROW(NotImplemented, "Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    { return Component::vaporPressure(T); }

    /*!
     * \brief The density [kg/m^3] of the component at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::gasDensity(temperature, pressure); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure the pure component in gas.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    { return Component::gasEnthalpy(temperature, pressure); }

    /*!
     * \brief The dynamic viscosity [Pa s] of the pure component at a given pressure and temperature.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::gasViscosity(temperature, pressure); }
};

template <class Scalar, class Component> class LiquidPhase
{
public:
    /*!
     * \brief A human readable name for the compoent.
     */
    static const char *name()
    { return Component::name(); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static Scalar molarMass()
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static Scalar criticalTemperature()
    {  DUNE_THROW(NotImplemented, "Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static Scalar criticalPressure()
    {  DUNE_THROW(NotImplemented, "Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static Scalar tripleTemperature()
    {  DUNE_THROW(NotImplemented, "Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static Scalar triplePressure()
    {  DUNE_THROW(NotImplemented, "Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    {  return Component::vaporPressure(T); }

    /*!
     * \brief The density [kg/m^3] of the component at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::liquidDensity(temperature, pressure); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure the pure component in liquid.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::liquidEnthalpy(temperature, pressure); }

    /*!
     * \brief The dynamic liquid viscosity [N/m^3*s] of the pure component.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::liquidViscosity(temperature, pressure); }
};

} // end namepace

#endif
