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
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_RFT_CONTAINER_HPP
#define OPM_RFT_CONTAINER_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Opm {

namespace data { class Wells; }
class EclipseState;
class Schedule;

/// Collection of cell-level RFT data--i.e., pressures and saturations--in
/// cells intersected by wells.  Will collect dynamic state values only for
/// those cells for which RFT data has been requested through the Schedule's
/// RFT configuration (member function ScheduleState::rft_config()).
///
/// \tparam FluidSystem Run's fluid system.  Needed, in particular, to infer
/// the simulator's \c Scalar type and its active phases.
template<class FluidSystem>
class RFTContainer {
    /// Simulator's floating-point element type.
    using Scalar = typename FluidSystem::Scalar;

    /// Convenience type alias for a linear sequence of floating-point
    /// values.
    using ScalarBuffer = std::vector<Scalar>;

    /// Phase index for gas.
    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;

    /// Phase index for oil.
    static constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;

    /// Phase index for water.
    static constexpr auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

public:
    /// Call-back function type for collecting cell level dynamic state
    /// values.
    using AssignmentFunc = std::function<Scalar()>;

    /// Call-back predicate type for inferring MPI characteristics for a
    /// named well.  Common examples of such characteristics are whether or
    /// not the current rank "owns" the well or if the well is intersected
    /// on the current rank.
    using WellQueryFunc = std::function<bool(const std::string&)>;

    /// Constructor.
    ///
    /// \param[in] eclState Run's static properties and configurations.
    ///
    /// \param[in] schedule Run's dynamic input objects such as its wells
    /// and RFT configuration.
    ///
    /// \param[in] wellIsOwnedByCurrent Predicate for whether or not a
    /// particular named well object is owned by the current rank.  The RFT
    /// state data for a particular well will be collected on the owning
    /// rank for that well.
    ///
    /// \param[in] wellOnCurrent Predicate for whether or not a particular
    /// named well object is intersected on the current rank.  This function
    /// will typically return \c true for those ranks that have local
    /// connections for the well and false otherwise.
    RFTContainer(const EclipseState& eclState,
                 const Schedule& schedule,
                 const WellQueryFunc& wellIsOwnedByCurrent,
                 const WellQueryFunc& wellOnCurrent)
        : eclState_(eclState)
        , schedule_(schedule)
        , wellIsOwnedByCurrentRank_(wellIsOwnedByCurrent)
        , wellOnCurrentRank_(wellOnCurrent)
    {}

    /// Export RFT cell level state data to requisite connections.
    ///
    /// Will populate the \code cell_pressure \endcode, the \code
    /// cell_saturation_water \endcode, and the \code cell_saturation_gas
    /// \endcode data members of the pertinent \code data::Connection
    /// \endcode objects depending on the run's active phases.
    ///
    /// \param[in,out] wellDatas Well and connection level dynamic report
    /// data at the current report step.  On exit, also contains all
    /// relevant and available RFT state data for the wells owned by the
    /// current rank.
    ///
    /// \param[in] reportStepNum Report step.  RFT data will be exported
    /// only for those wells that request such values at this time.
    ///
    /// \param[in] comm MPI communication object.  Needed to collect RFT
    /// state data from all ranks that share a well.
    void addToWells(data::Wells& wellDatas,
                    const std::size_t reportStepNum,
                    const Parallel::Communication& comm);

    /// Prepare internal data structures to collect RFT state.
    ///
    /// \param[in] reportStepNum Report step.  RFT data will be exported
    /// only for those wells that request such values at this time.
    void allocate(const std::size_t reportStepNum);

    /// Collect cell level RFT state into internal data structures.
    ///
    /// Does nothing if the cell is not among those for which RFT state is
    /// requested.
    ///
    /// \param[in] cartesianIndex Linearised global cell ID.
    ///
    /// \param[in] oil Call-back for transferring cell level pressure values
    /// in cell \p cartesianIndex from the simulator and into the container.
    /// Will be invoked only if oil is active in the current run.
    ///
    /// \param[in] water Call-back for transferring cell level water
    /// saturation values in cell \p cartesianIndex from the simulator and
    /// into the container.  Will be invoked only if water is active in the
    /// current run.
    ///
    /// \param[in] gas Call-back for transferring cell level gas saturation
    /// values in cell \p cartesianIndex from the simulator and into the
    /// container.  Will be invoked only if gas is active in the current
    /// run.
    void assign(const unsigned cartesianIndex,
                const AssignmentFunc& oil,
                const AssignmentFunc& water,
                const AssignmentFunc& gas);

private:
    /// Run's static properties and configurations.
    const EclipseState& eclState_;

    /// Run's dynamic input objects, e.g., wells and RFT configuration.
    const Schedule& schedule_;

    /// Predicate function for whether or not a particular named well is
    /// owned by the current rank.
    ///
    /// RFT state data will be collected on the well's owning rank.
    WellQueryFunc wellIsOwnedByCurrentRank_;

    /// Predicate function for whether or not a particular named well is
    /// intersected on the current rank.
    ///
    /// Needed to properly allocate the internal buffers.
    WellQueryFunc wellOnCurrentRank_;

    /// Cell level oil pressure values for all pertinent well connections.
    ///
    /// Keyed by the linearised global cell ID.  Will be populated only if
    /// oil is active in the current run.
    std::map<std::size_t, Scalar> oilConnectionPressures_;

    /// Cell level water saturation values for all pertinent well
    /// connections.
    ///
    /// Keyed by the linearised global cell ID.  Will be populated only if
    /// water is active in the current run.
    std::map<std::size_t, Scalar> waterConnectionSaturations_;

    /// Cell level gas saturation values for all pertinent well connections.
    ///
    /// Keyed by the linearised global cell ID.  Will be populated only if
    /// gas is active in the current run.
    std::map<std::size_t, Scalar> gasConnectionSaturations_;
};

} // namespace Opm

#endif // OPM_RFT_CONTAINER_HPP
