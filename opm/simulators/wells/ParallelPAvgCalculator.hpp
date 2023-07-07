/*
  Copyright 2023 Equinor ASA.

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

#ifndef PARALLEL_PAVG_CALCULATOR_HPP
#define PARALLEL_PAVG_CALCULATOR_HPP

#include <opm/input/eclipse/Schedule/Well/PAvgCalculator.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <functional>

namespace Opm {
    class GridDims;
    class WellConnections;
} // namespace Opm

namespace Opm {

/// Facility for deriving well-level pressure values from selected
/// block-averaging procedures.  Applicable to stopped wells which don't
/// have a flowing bottom-hole pressure.  Mainly useful for reporting.
///
/// Parallel edition.  Handles distributed wells.
class ParallelPAvgCalculator : public PAvgCalculator
{
public:
    /// Constructor
    ///
    /// \param[in] comm MPI communication object.  Typically \code
    ///   ParallelWellInfo::communication() \endcode.
    ///
    /// \param[in] cellIndexMap Cell index triple map ((I,J,K) <-> global).
    ///
    /// \param[in] connections List of reservoir connections for single
    ///   well.
    ParallelPAvgCalculator(const Parallel::Communication& comm,
                           const GridDims&                cellIndexMap,
                           const WellConnections&         connections);

private:
    /// MPI communication object.
    std::reference_wrapper<const Parallel::Communication> comm_;

    /// Communicate local contributions and collect global (off-rank)
    /// contributions.
    ///
    /// Reads from and writes to base class data members \c accumCTF_ and \c
    /// accumPV_.
    void collectGlobalContributions() override;
};

} // namespace Opm

#endif // PARALLEL_PAVG_CALCULATOR_HPP
