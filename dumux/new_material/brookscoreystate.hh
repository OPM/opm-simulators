/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \brief Specification of the state API for the Brooks-Corey
 *        capillary pressure model.
 */
#ifndef DUMUX_BROOKS_COREY_STATE_HH
#define DUMUX_BROOKS_COREY_STATE_HH

#include <dumux/auxiliary/apis.hh>

#include <dumux/auxiliary/statehelpermacros.hh>

namespace Dune
{
namespace Api
{
BEGIN_API_DEF(BrooksCoreyParams)
{
    typedef typename Implementation::Scalar Scalar;
    Scalar tmp;
    tmp = const_impl.alpha();
    tmp = const_impl.pe();
}
END_API_DEF;

BEGIN_API_DEF(BrooksCoreyState)
{
    typedef typename Implementation::Scalar Scalar;
    Scalar tmp = 0.5;
    impl.setAlpha(tmp);
    impl.setPe(tmp);
}
END_API_DEF;
}; // namespace Api

/*!
 * \brief A reference implementation of the state API class for the
 *        Brooks-Corey Sw-pC relation.
 */
template <class ScalarT>
class BrooksCoreyState
{
public:
    typedef ScalarT Scalar;

    BrooksCoreyState(Scalar pe, Scalar alpha)
        : pe_(pe), alpha_(alpha)
    {
    }

    //! The entry pressure
    PROPERTY(Scalar, pe, setPe);

    //! The alpha shape parameter
    PROPERTY(Scalar, alpha, setAlpha);
};
}; // namespace Dune

#endif
