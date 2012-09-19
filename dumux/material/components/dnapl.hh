// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Felix Bode                                        *
 *   Copyright (C) 2010 by Katherina Baber                                   *
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
 * \brief A simple implementation of a dense non-aqueous phase liquid (DNAPL).
 *
 * The parameters are chosen to roughly correspond to those of
 * trichloroethylene (TCE) at standard conditions.
 */
#ifndef DUMUX_DNAPL_HH
#define DUMUX_DNAPL_HH

#include <dumux/material/idealgas.hh>
#include "component.hh"

namespace Dumux {
/*!
 * \ingroup Components
 *
 * \brief A simple implementation of a dense non-aqueous phase liquid (DNAPL).
 *
 * The parameters are chosen to roughly correspond to those of
 * trichloroethylene (TCE) at standard conditions.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class DNAPL : public Component<Scalar, DNAPL<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the TCE.
     */
    static const char *name()
    { return "DNAPL_TCE"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of TCE.
     */
    static constexpr Scalar molarMass()
    { return 131.39e-3; /* [kg/mol] */ };

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of TCE.
     */
    static Scalar criticalTemperature()
    { DUNE_THROW(Dune::NotImplemented, "criticalTemperature for TCE"); };

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of TCE.
     */
    static Scalar criticalPressure()
    { DUNE_THROW(Dune::NotImplemented, "criticalPressure for TCE"); };

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at TCE's triple point.
     */
    static Scalar tripleTemperature()
    { DUNE_THROW(Dune::NotImplemented, "tripleTemperature for TCE"); };

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at TCE's triple point.
     */
    static Scalar triplePressure()
    { DUNE_THROW(Dune::NotImplemented, "triplePressure for TCE"); };

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure TCE
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { return 3900; /* [Pa] (at 20C) */ };

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief The density of steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas<Scalar>::density(molarMass(),
                                         temperature,
                                         pressure);
    };

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The density of pure TCE at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static constexpr Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { return 1460.0; /* [kg/m^3] */ }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure TCE.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static constexpr Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { return 5.7e-4; /* [Pa s] */ };
};

} // end namepace

#endif
