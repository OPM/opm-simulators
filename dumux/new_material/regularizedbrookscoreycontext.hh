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
 * \brief Specification of the context API for the Brooks-Corey
 *        capillary pressure model.
 */
#ifndef DUMUX_REGULARIZED_BROOKS_COREY_CONTEXT_HH
#define DUMUX_REGULARIZED_BROOKS_COREY_CONTEXT_HH

#include <dumux/new_material/statehelpermacros.hh>

namespace Dune
{
/*!
 * \brief A reference implementation of the context API class for the
 *        Brooks-Corey Sw-pC relation.
 */
template <class ScalarT>
class RegularizedBrooksCoreyContext
{
public:
    typedef ScalarT Scalar;
    
    RegularizedBrooksCoreyContext(Scalar pe = 0, Scalar alpha = 0)
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
