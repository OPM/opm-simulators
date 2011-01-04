// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief An adapter call for a pure two-phase system,
 * where the wetting and the non-wetting phase can be defined.
 */
#ifndef DUMUX_2P_SYSTEM_HH
#define DUMUX_2P_SYSTEM_HH

#include "liquidphase.hh"
#include "gasphase.hh"

#include <dune/common/exceptions.hh>
#include "defaultcomponents.hh"

#include <dumux/common/propertysystem.hh>

namespace Dumux {
// forward defintions of the property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(WettingPhase);
NEW_PROP_TAG(NonwettingPhase);
};


/*!
 * \ingroup Fluidsystems
 *
 * \brief This is an adapter class for a pure two-phase system.
 *
 * In this adapter class, two fluid phases of the type liquidphase or gasphase
 * can be defined. These phases consist of one pure component. With the help of
 * this adapter class, the phase properties can be accessed. This is suitable
 * for pure two-phase systems without compositional effects.
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
     *
     * \param phaseIdx index of the phase
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
     *
     * \param phaseIdx index of the phase
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
     *
     * \param phaseIdx index of the phase
     * \param compIdx index of the component
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \return density of the phase, if phaseIdx = compIdx, otherwise 0
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
     * \brief Return the density of a phase [kg/m^3].
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState The fluid state of the two-phase model
     * \tparam FluidState the fluid state class of the two-phase model
     * \return returns the density of the phase
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
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
     * \brief Return the viscosity of a phase [Pa*s].
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState The fluid state of the two-phase model
     * \tparam FluidState the fluid state class of the two-phase model
     * \return returns the viscosity of the phase [Pa*s]
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
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
     * \brief Return the specific enthalpy of a phase [J/kg].
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState The fluid state of the two-phase model
     * \tparam FluidState the fluid state class of the two-phase model
     * \return returns the specific enthalpy of the phase [J/kg]
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &fluidState)
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
     * \brief Return the specific internal energy of a phase [J/kg].
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState The fluid state of the two-phase model
     * \tparam FluidState the fluid state class of the two-phase model
     * \return returns the specific internal energy of the phase [J/kg]
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                      Scalar temperature,
                                      Scalar pressure,
                                      const FluidState &fluidState)
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
