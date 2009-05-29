/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
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
 * \file PlayTypeSate.hh Specification of the state API for a twophase
 *                       playtype hystesis model.
 */
#ifndef DUMUX_PLAY_TYPE_STATE_HH
#define DUMUX_PLAY_TYPE_STATE_HH

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/statehelpermacros.hh>

namespace Dune
{
namespace Api
{
BEGIN_API_DEF(PlayTypeParams)
{
    typedef typename Implementation::Scalar Scalar;
    typedef typename Implementation::CapPressureParams PCParams;
    require<TwophaseSatParams>(impl);

    const PCParams *tmp;
    tmp = &const_impl.micParams();
    tmp = &const_impl.mdcParams();

    Scalar tmp2 = 0.5;
    tmp2   = const_impl.SweRef();
    tmp2   = const_impl.deltaSwe();
    bool b;
    b = const_impl.isImbib();
}
END_API_DEF;

BEGIN_API_DEF(PlayTypeState)
{
    typedef typename Implementation::Scalar Scalar;
    typedef typename Implementation::CapPressureParams PCParams;
    require<PlayTypeParams>(impl);

    Scalar tmp = 0.5;
    impl.setSweRef(tmp);
    impl.setImbib(true);
}
END_API_DEF;
}; // namespace Api
}; // namespace Dune

#endif
