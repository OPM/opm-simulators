// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief Base class for Peng-Robinson parameters of a
 *        single-component fluid or a mixture
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_HH
#define DUMUX_PENG_ROBINSON_PARAMS_HH

#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \brief Stores and provides access to the Peng-Robinson parameters
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
template <class Scalar>
class PengRobinsonParams
{
public:
    /*!
     * \brief Returns the attractive parameter 'a' of the
     *        Peng-Robinson fluid.
     */
    Scalar a() const
    { return a_; }

    /*!
     * \brief Returns the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     */
    Scalar b() const
    { return b_; }

    /*!
     * \brief If run under valgrind, this method produces an warning
     *        if the parameters where not determined correctly.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        Valgrind::CheckDefined(a_);
        Valgrind::CheckDefined(b_);
#endif
    };

    /*!
     * \brief Set the attractive parameter 'a' of the Peng-Robinson
     *        fluid.
     */
    void setA(Scalar value)
    { a_ = value; }

    /*!
     * \brief Set the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     */
    void setB(Scalar value)
    { b_ = value; }

protected:
    Scalar a_;
    Scalar b_;
};

} // end namepace

#endif
