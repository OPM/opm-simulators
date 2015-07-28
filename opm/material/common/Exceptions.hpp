// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \brief Provides some exception classes.
 */
#ifndef OPM_MATERIAL_EXCEPTIONS_HPP
#define OPM_MATERIAL_EXCEPTIONS_HPP

#include <stdexcept>

// some additional exception classes
namespace Opm {

/*!
 * \brief Exception that indicates that a functionality has not been implemented yet.
 *
 * A better name for this class would be 'NotImplemented', but this would conflict with
 * opm-core's exception class of the same name.
 */
class NotAvailable : public std::logic_error
{
public:
    explicit NotAvailable(const std::string &message)
        : std::logic_error(message)
    {}
};

/*!
 * \brief Exception that indicates that something went wrong during a numerical procedure
 *
 * For example this exception should be thrown if a value is out of range, an unexpected
 * NaN was encountered, etc. Usually such errors are non-fatal, and the numerical
 * proceedure can be restarted with more conservertive parameters like a smaller time
 * step.
 */
class NumericalIssue : public std::runtime_error
{
public:
    explicit NumericalIssue(const std::string &message)
        : std::runtime_error(message)
    {}
};
}

#endif // OPM_MATERIAL_EXCEPTIONS_HPP
