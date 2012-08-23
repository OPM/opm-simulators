// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
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
 * \ingroup Components
 * \brief  A simple implementation of a light non-aqueous phase liquid (LNAPL).
 */
#ifndef DUMUX_LNAPL_HH
#define DUMUX_LNAPL_HH

#include "component.hh"

namespace Dumux {
/*!
 * \ingroup Components
 *
 * \brief A simple implementation of a LNAPL, e.g. a kind of oil
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class LNAPL : public Component<Scalar, LNAPL<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the LNAPL.
     */
    static const char *name()
    { return "LNAPL"; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Rough estimate of the density of oil \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static constexpr Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { return 890; }

    /*!
     * \brief Rough estimate of the viscosity of oil in \f$\mathrm{[Pa*s]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static constexpr Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { return 8e-3; };
};

} // end namepace

#endif
