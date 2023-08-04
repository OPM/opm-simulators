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

#ifndef PRESSURE_AVERAGE_HPP
#define PRESSURE_AVERAGE_HPP

#include <vector>

namespace Opm {
namespace detail {

//! \brief Calculates average pressure value.
template<class Scalar>
Scalar pressureAverage(const Scalar pressurePvHydrocarbon,
                       const Scalar pvHydrocarbon,
                       const Scalar pressurePv,
                       const Scalar pv,
                       const bool hydrocarbon);

//! \brief Calculates average pressure value for a vector.
template<class Scalar>
std::vector<Scalar>
pressureAverage(const std::vector<Scalar>& pressurePvHydrocarbon,
                const std::vector<Scalar>& pvHydrocarbon,
                const std::vector<Scalar>& pressurePv,
                const std::vector<Scalar>& pv,
                const bool hydrocarbon);

} // namespace detail
} // namespace Opm

#endif // PRESSURE_AVERAGE_HPP
