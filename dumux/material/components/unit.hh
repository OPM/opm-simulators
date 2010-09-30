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
 * \ingroup Components
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_UNIT_HH
#define DUMUX_UNIT_HH


#include "component.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief Rough estimate for testing purposes of water.
 *
 * \tparam Scalar  The type used for scalar values
 */
template <class Scalar>
class Unit : public Component<Scalar, Unit<Scalar> >
{
    typedef Component<Scalar, Unit<Scalar> > ParentType;

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "Unit"; }

    /*!
     * \brief Rough estimate of the density of oil [kg/m^3].
     *
     * \param temperature temperature of component
     * \param pressure pressure of component
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1.0;
    }

    /*!
     * \brief Rough estimate of the viscosity of oil kg/(ms).
     *
     * \param temperature temperature of component
     * \param pressure pressure of component
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1.0;
    };

};

} // end namepace

#endif
