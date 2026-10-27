/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_NUMERICAL_AQUIFER_EQUILIBRATION_CHECK_HPP
#define OPM_NUMERICAL_AQUIFER_EQUILIBRATION_CHECK_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <vector>

namespace Opm {

class EclipseState;

/// Connection of a numerical aquifer cell to a cell of another
/// equilibration region with a significant initial potential difference of
/// a phase which is mobile in the upstream cell.
struct NumericalAquiferEquilibrationMismatch
{
    std::size_t aquiferId{};
    int aquiferRegion{};          //!< EQLNUM region of the aquifer cell
    int cellRegion{};             //!< EQLNUM region of the connected cell
    double potentialDifference{}; //!< Phase potential, aquifer minus cell [Pa]
};

/// Initial potential differences [Pa] below this limit are not considered
/// to give an initial transient.
inline constexpr double numericalAquiferEquilibrationTolerance = 1.0e4;

/// Reports numerical aquifers which are not in hydrostatic equilibrium with
/// the cells they are connected to, because the cells belong to different
/// equilibration regions.
///
/// Collective.  The mismatches found on all processes are gathered, and the
/// I/O rank logs a description which includes the deck changes resolving
/// the issue, as an error, or as a warning if \p allowMismatch is true.
///
/// \return Total number of mismatched connections on all processes.
std::size_t
reportNumericalAquiferEquilibration(const std::vector<NumericalAquiferEquilibrationMismatch>& localMismatches,
                                    const EclipseState& eclState,
                                    const Parallel::Communication& comm,
                                    bool allowMismatch);

} // namespace Opm

#endif // OPM_NUMERICAL_AQUIFER_EQUILIBRATION_CHECK_HPP
