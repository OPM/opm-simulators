// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \brief Some exceptions thrown by opm-material.
 */
#ifndef OPM_EXCEPTIONS_HH
#define OPM_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup Exception
 * \brief Exception thrown if a fixable numerical problem occurs.
 *
 * (e.g. time step too big, etc.)
 */
class NumericalProblem : public Dune::Exception
{
public:
    // copy constructor
    NumericalProblem(const NumericalProblem &v)
        : Dune::Exception(v)
    {}

    // default constructor
    NumericalProblem()
    {}

    // constructor with error message
    NumericalProblem(const std::string &s)
    { this->message(s); }
};

/*!
 * \ingroup Exception
 * \brief Exception thrown if a run-time parameter is not specified correctly.
 */
class ParameterException : public Dune::Exception
{
public:
    // copy constructor
    ParameterException(const ParameterException &v)
        : Dune::Exception(v)
    {}

    // default constructor
    ParameterException()
    {}

    // constructor with error message
    ParameterException(const std::string &s)
    { this->message(s); }
};

} // namespace Opm

#endif
