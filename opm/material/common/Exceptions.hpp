/*****************************************************************************
 *   Copyright (C) 2013 by Andreas Lauser                                    *
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
 * A better name for this class would be 'NotAvailable', but this would conflict with
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
