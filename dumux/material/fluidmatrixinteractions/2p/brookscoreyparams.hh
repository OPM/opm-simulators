// $Id$
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
 * \ingroup fluidmatrixinteractions
 *  \defgroup fluidmatrixinteractionsparams FluidMatrixInteractions Parameters
 */

/*!
 * \file
 *
 * \brief Specification of the material parameters
 *       for the Brooks Corey constitutive relations.
 */
#ifndef DUMUX_BROOKS_COREY_PARAMS_HH
#define DUMUX_BROOKS_COREY_PARAMS_HH

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \brief Specification of the material parameters
 *       for the Brooks Corey constitutive relations.
 *
 *        \ingroup fluidmatrixinteractionsparams
 *
 *\see BrooksCorey
 */
template <class ScalarT>
class BrooksCoreyParams
{
public:
    typedef ScalarT Scalar;

    BrooksCoreyParams()
    {
        Valgrind::SetUndefined(*this);
    }

    BrooksCoreyParams(Scalar pe, Scalar alpha)
        : pe_(pe), alpha_(alpha)
    {
    }

    /*!
     * \brief Returns the entry pressure [Pa]
     */
    Scalar pe() const
    { return pe_; }

    /*!
     * \brief Set the entry pressure [Pa]
     */
    void setPe(Scalar v)
    { pe_ = v; }


    /*!
     * \brief Returns the alpha shape parameter
     */
    Scalar alpha() const
    { return alpha_; }

    /*!
     * \brief Set the alpha shape parameter
     */
    void setAlpha(Scalar v)
    { alpha_ = v; }

private:
    Scalar pe_;
    Scalar alpha_;
};
}; // namespace Dumux

#endif
