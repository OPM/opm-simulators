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
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_OIL_HH
#define DUMUX_OIL_HH

#include <dune/common/exceptions.hh>

#include "component.hh"

namespace Dune
{
/*!
 * \brief Rough estimate for testing purposes of some oil.
 */
template <class Scalar>
class Oil : public Component<Scalar, Oil<Scalar> >
{
    typedef Component<Scalar, Oil<Scalar> > ParentType;

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "Oil"; }

    /*!
     * \brief Rough estimate of the density of oil [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 890;
    }

    /*!
     * \brief Rough estimate of the viscosity of oil kg/(ms).
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 8e-3;
    };

};

} // end namepace

#endif
