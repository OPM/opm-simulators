// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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

    BrooksCoreyParams(Scalar pe, Scalar lambda)
        : pe_(pe), lambda_(lambda)
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
     * \brief Returns the lambda shape parameter
     */
    Scalar lambda() const
    { return lambda_; }

    /*!
     * \brief Set the lambda shape parameter
     */
    void setLambda(Scalar v)
    { lambda_ = v; }

private:
    Scalar pe_;
    Scalar lambda_;
};
} // namespace Dumux

#endif
