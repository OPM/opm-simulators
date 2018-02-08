// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Provides the opm-material specific exception classes.
 */
#ifndef OPM_MATERIAL_EXCEPTIONS_HPP
#define OPM_MATERIAL_EXCEPTIONS_HPP

#if HAVE_OPM_COMMON
#include <opm/common/Exceptions.hpp>
#endif

#include <dune/common/exceptions.hh>

#include <stdexcept>

// the opm-material specific exception classes
namespace Opm {
class NumericalIssue
#if HAVE_OPM_COMMON
    : public Opm::NumericalProblem
#else
    : public std::runtime_error
#endif
{
public:
    explicit NumericalIssue(const std::string &message)
#if HAVE_OPM_COMMON
      : Opm::NumericalProblem(message)
#else
      : std::runtime_error(message)
#endif
    {}
};
}

#endif // OPM_MATERIAL_EXCEPTIONS_HPP
