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

#ifndef OPM_SIMULATOR_UTILS_HPP
#define OPM_SIMULATOR_UTILS_HPP

#include <string>
#include <vector>

namespace Opm {

/*!
* \brief Given a time step size in seconds, return it in a format which is more
*        easily parsable by humans.
*
* e.g. 874000.0 will become "10.12 days"
*/
std::string humanReadableTime(double timeInSeconds, bool isAmendment = true);

/*!
 * \brief Read explicitly defined time steps from file.
 * \param file File to read
 */
template<class Scalar>
std::vector<Scalar> readTimeStepFile(const std::string& file);

} // namespace Opm

#endif
