/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
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
 * \brief A fluid system with one phase and an arbitrary number of components.
 */
#ifndef DUMUX_ISFLUID_TRAIL_SYSTEM_HH
#define DUMUX_ISFLUID_TRAIL_SYSTEM_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/boxmodels/1p2c/1p2cproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(OnePTwoCIndices);
};

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system with one phase and two components representing an
 *        intersticial fluid that contains a medical trail.
 *
 * A fluid system with one phase and two components representing an
 * intersticial fluid that contains a medical trail. This is used for
 * in the 1p2c model.
 */
template <class TypeTag, bool verbose=true>
class ISFluid_Trail_System
{
    typedef ISFluid_Trail_System<TypeTag, verbose> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

public:
    enum {
        // component indices
        isFluidIdx = 0,
        trailIdx = 1
    };

    static void init()
    {}

    /*!
     * \brief Return the human readable name of a component
     * \param compIdx the component index
     */
    static const char *componentName(int compIdx)
    {
        switch(compIdx)
        {
        case isFluidIdx:
            return "ISFluid";
        case trailIdx:
            return "Trail";
        default:
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Return the molar mass of a component [kg/mol].
     *
     * \param compIdx the component index
     *
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case isFluidIdx:
            // TODO: this is just a rough guess
            return 22e-3; // [kg/mol]
        case trailIdx:
            return 567e-3; // [kg/mol]
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }


    /*!
     * \brief Return the (constant) phase density [kg/m^3].
     *
     * \param phaseIdx index of the phase, only 0 is accepted (single-phase system)
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState the fluid state
     *
     * \tparam FluidState the fluid state class
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 1.03e3; // in [kg /m^3]

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the (constant) dynamic viscosity of a phase [Pa/s].
     *
     * \param phaseIdx index of the phase, only 0 is accepted (single-phase system)
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState the fluid state
     *
     * \tparam FluidState the fluid state class
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 0.00069152; // in [Pa*s]

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the (constant) binary diffusion coefficient [m^2/s].
     *
     * \param phaseIdx index of the phase, only 0 is accepted (single-phase system)
     * \param compIIdx index of the first component in the binary mixture
     * \param compJIdx index of second component in the binary mixture
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState the fluid state
     *
     * \tparam FluidState the fluid state class
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    {
        // 3.7378e-12
        return 8.8786695-14; // in [m^2/s]
    }
};

} // end namepace

#endif
