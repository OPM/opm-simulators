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

#include <config.h>
#include <opm/simulators/utils/PressureAverage.hpp>

#include <cassert>
#include <cstddef>

namespace Opm {
namespace detail {

template<class Scalar>
Scalar
pressureAverage(const Scalar pressurePvHydrocarbon,
                const Scalar pvHydrocarbon,
                const Scalar pressurePv,
                const Scalar pv,
                const bool hydrocarbon)
{
    if (hydrocarbon && (pvHydrocarbon > 1e-10))
        return pressurePvHydrocarbon / pvHydrocarbon;

    return pressurePv / pv;
}

template<class Scalar>
std::vector<Scalar>
pressureAverage(const std::vector<Scalar>& pressurePvHydrocarbon,
                const std::vector<Scalar>& pvHydrocarbon,
                const std::vector<Scalar>& pressurePv,
                const std::vector<Scalar>& pv,
                const bool hydrocarbon)
{
    const std::size_t size = pressurePvHydrocarbon.size();
    assert(pvHydrocarbon.size() == size);
    assert(pressurePv.size() == size);
    assert(pv.size() == size);

    std::vector<Scalar> fraction(size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        fraction[i] = pressureAverage(pressurePvHydrocarbon[i],
                                      pvHydrocarbon[i],
                                      pressurePv[i],
                                      pv[i],
                                      hydrocarbon);
    }

    return fraction;
}

template double
pressureAverage<double>(const double,
                        const double,
                        const double,
                        const double,
                        const bool);

template std::vector<double>
pressureAverage<double>(const std::vector<double>&,
                        const std::vector<double>&,
                        const std::vector<double>&,
                        const std::vector<double>&,
                        const bool);

} // namespace detail
} // namespace Opm
