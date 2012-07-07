/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
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
 * \brief A fluid system with one phase and an arbitrary number of components.
 */
#ifndef DUMUX_ISFLUID_TRAIL_SYSTEM_HH
#define DUMUX_ISFLUID_TRAIL_SYSTEM_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/boxmodels/1p2c/1p2cproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
};

/*!
 * \brief A fluid system with one phase and an arbitrary number of components.
 */
template <class TypeTag, bool verbose=true>
class ISFluid_Trail_System
{
    typedef ISFluid_Trail_System<TypeTag, verbose> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

    enum {
        isfluid = Indices::konti,
        trail = Indices::transport
    };

public:
    static void init()
    {}

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        switch(compIdx)
        {
        case isfluid:
            return "ISFluid";
        case trail:
            return "Trail";
        default:
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 1.03e-9; // in [kg /microm^3]
                            // Umwandlung in Pascal notwendig -> density*10^6
                            // will man wieder mit kg/m^3 rechnen muss g auch wieder geändert werden

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar molarDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 3.035e-16;   // in [mol/microm^3] = 303.5 mol/m^3

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 0.00069152; // in [Pas]

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent of a component in a phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    {
        return 0.088786695; // in [microm^2/s] = 3.7378e-12 m^2/s
    }
};

} // end namepace

#endif
