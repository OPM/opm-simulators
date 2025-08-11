/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2017 IRIS AS

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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>

#include <opm/output/data/Wells.hpp>

#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellContainer.hpp>

#include <opm/simulators/utils/BlackoilPhases.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

template<class Scalar> class ParallelWellInfo;
template<class Scalar> struct PerforationData;
template<class Scalar> class ConnFracStatistics;
class Schedule;
enum class WellStatus;

/// The state of a set of wells, tailored for use by the fully
/// implicit blackoil simulator.
template<typename FluidSystem, typename Indices>
class WellState
{
public:
    using Scalar = typename FluidSystem::Scalar;

    static const std::uint64_t event_mask = ScheduleEvents::WELL_STATUS_CHANGE
        | ScheduleEvents::PRODUCTION_UPDATE
        | ScheduleEvents::INJECTION_UPDATE;

    // TODO: same definition with WellInterface, eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;

    // Only usable for testing purposes
    explicit WellState(const ParallelWellInfo<Scalar>& pinfo);

    WellState() = default;

    // explicit WellState(const PhaseUsage& pu)
    //     : phase_usage_(pu)
    // {}

    static WellState serializationTestObject(const ParallelWellInfo<Scalar>& pinfo);

    std::size_t size() const
    {
        return this->wells_.size();
    }

    std::vector<std::string> wells() const
    {
        return this->wells_.wells();
    }

    int numWells() const
    {
        return this->size();
    }

    const ParallelWellInfo<Scalar>& parallelWellInfo(std::size_t well_index) const;

    /// Allocate and initialize if wells is non-null.  Also tries
    /// to give useful initial values to the bhp(), wellRates()
    /// and perfPhaseRatesORG() fields, depending on controls
    void init(const std::vector<Scalar>& cellPressures,
              const std::vector<Scalar>& cellTemperatures,
              const Schedule& schedule,
              const std::vector<Well>& wells_ecl,
              const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
              const int report_step,
              const WellState* prevState,
              const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
              const SummaryState& summary_state,
              const bool enableDistributedWells);

    void resize(const std::vector<Well>& wells_ecl,
                const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
                const Schedule& schedule,
                const bool handle_ms_well,
                const std::size_t numCells,
                const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
                const SummaryState& summary_state,
                const bool enable_distributed_wells);

    void setCurrentWellRates(const std::string& wellName,
                             const std::vector<Scalar>& new_rates)
    {
        auto& [owner, rates] = this->well_rates.at(wellName);
        if (owner)
            rates = new_rates;
    }

    const std::vector<Scalar>& currentWellRates(const std::string& wellName) const;

    bool hasWellRates(const std::string& wellName) const
    {
        return this->well_rates.find(wellName) != this->well_rates.end();
    }

    void clearWellRates()
    {
        this->well_rates.clear();
    }

    void gatherVectorsOnRoot(const std::vector<data::Connection>& from_connections,
                             std::vector<data::Connection>& to_connections,
                             const Parallel::Communication& comm) const;

    data::Wells
    report(const int* globalCellIdxMap,
           const std::function<bool(const int)>& wasDynamicallyClosed) const;

    void reportConnections(std::vector<data::Connection>& connections,
                           const PhaseUsage &pu,
                           std::size_t well_index,
                           const int* globalCellIdxMap) const;

    /// init the MS well related.
    void initWellStateMSWell(const std::vector<Well>& wells_ecl,
                             const WellState* prev_well_state);

    static void calculateSegmentRates(const ParallelWellInfo<Scalar>&      pw_info,
                                      const std::vector<std::vector<int>>& segment_inlets,
                                      const std::vector<std::vector<int>>& segment_perforations,
                                      const std::vector<Scalar>&           perforation_rates,
                                      const int                            np,
                                      const int                            segment,
                                      std::vector<Scalar>&                 segment_rates);


    void communicateGroupRates(const Parallel::Communication& comm);

    void updateGlobalIsGrup(const Parallel::Communication& comm);
    void updateEfficiencyScalingFactor(const std::string& wellName,
                                       const Scalar value);

    bool isInjectionGrup(const std::string& name) const
    {
        return this->global_well_info.value().in_injecting_group(name);
    }

    bool isProductionGrup(const std::string& name) const
    {
        return this->global_well_info.value().in_producing_group(name);
    }

    bool isOpen(const std::string& name) const
    {
        return this->global_well_info.value().is_open(name);
    }

    Scalar getGlobalEfficiencyScalingFactor(const std::string& name) const
    {
        return this->global_well_info.value().efficiency_scaling_factor(name);
    }

    // If the ALQ has changed since the previous time step,
    // reset current_alq and update default_alq. ALQ is used for
    // constant lift gas injection and for gas lift optimization
    // (THP controlled wells).
    void updateWellsDefaultALQ(const Schedule& schedule,
                              const int report_step,
                              const SummaryState& summary_state);

    void gliftTimeStepInit()
    {
        for (size_t i = 0; i < this->size(); ++i) {
            this->wells_[i].alq_state.reset_count();
        }
    }


    int wellNameToGlobalIdx(const std::string& name)
    {
        return this->global_well_info.value().well_index(name);
    }

    std::string globalIdxToWellName(const int index)
    {
        return this->global_well_info.value().well_name(index);
    }

    bool wellIsOwned(std::size_t well_index,
                     const std::string& wellName) const;

    bool wellIsOwned(const std::string& wellName) const;

    bool isRank0() const {
        return this->global_well_info.value().isRank0();
    }

    void updateStatus(int well_index, WellStatus status);

    void openWell(int well_index);
    void shutWell(int well_index);
    void stopWell(int well_index);

    /// The number of phases present.
    constexpr int numPhases() const
    {
        return Indices::numPhases;
    }

    // const PhaseUsage& phaseUsage() const
    // {
    //     return this->phase_usage_;
    // }

    /// One rate per well and phase.
    std::vector<Scalar>& wellRates(std::size_t well_index)
    { return this->wells_[well_index].surface_rates; }
    const std::vector<Scalar>& wellRates(std::size_t well_index) const
    { return this->wells_[well_index].surface_rates; }

    const std::string& name(std::size_t well_index) const
    {
        return this->wells_.well_name(well_index);
    }

    std::optional<std::size_t> index(const std::string& well_name) const
    {
        return this->wells_.well_index(well_name);
    }

    const SingleWellState<FluidSystem, Indices>& operator[](std::size_t well_index) const
    {
        return this->wells_[well_index];
    }

    const SingleWellState<FluidSystem, Indices>& operator[](const std::string& well_name) const
    {
        return this->wells_[well_name];
    }

    SingleWellState<FluidSystem, Indices>& operator[](std::size_t well_index)
    {
        return this->wells_[well_index];
    }

    SingleWellState<FluidSystem, Indices>& operator[](const std::string& well_name)
    {
        return this->wells_[well_name];
    }

    const SingleWellState<FluidSystem, Indices>& well(std::size_t well_index) const
    {
        return this->operator[](well_index);
    }

    const SingleWellState<FluidSystem, Indices>& well(const std::string& well_name) const
    {
        return this->operator[](well_name);
    }

    SingleWellState<FluidSystem, Indices>& well(std::size_t well_index)
    {
        return this->operator[](well_index);
    }

    SingleWellState<FluidSystem, Indices>& well(const std::string& well_name)
    {
        return this->operator[](well_name);
    }

    bool has(const std::string& well_name) const
    {
        return this->wells_.has(well_name);
    }

    bool operator==(const WellState&) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(well_rates);
        if (serializer.isSerializing()) {
            serializer(wells_.size());
        } else {
            std::size_t size = 0;
            serializer(size);
            if (size != wells_.size()) {
                OPM_THROW(std::runtime_error, "Error deserializing WellState: size mismatch");
            }
        }
        for (auto& w : wells_) {
            serializer(w);
        }
        serializer(permanently_inactive_well_names_);
    }

    bool is_permanently_inactive_well(const std::string& wname) const {
        return std::find(this->permanently_inactive_well_names_.begin(), this->permanently_inactive_well_names_.end(), wname) != this->permanently_inactive_well_names_.end();
    }

private:
    bool enableDistributedWells_ = false;

//    PhaseUsage phase_usage_;

    // The wells_ variable is essentially a map of all the wells on the current
    // process. Observe that since a well can be split over several processes a
    // well might appear in the WellContainer on different processes.
    WellContainer<SingleWellState<FluidSystem, Indices>> wells_;

    // The members global_well_info and well_rates are map like
    // structures which will have entries for *all* the wells in the system.

    // Use of std::optional<> here is a technical crutch, the
    // WellStateFullyImplicitBlackoil class should be default constructible,
    // whereas the GlobalWellInfo is not.
    std::optional<GlobalWellInfo<Scalar>> global_well_info;

    // The well_rates variable is defined for all wells on all processors. The
    // bool in the value pair is whether the current process owns the well or
    // not.
    std::map<std::string, std::pair<bool, std::vector<Scalar>>> well_rates;

    // Keep track of permanently inactive well names
    std::vector<std::string> permanently_inactive_well_names_;

    data::Segment
    reportSegmentResults(const int         well_id,
                         const int         seg_ix,
                         const int         seg_no) const;


    /// Allocate and initialize if wells is non-null.
    /// Also tries to give useful initial values to the bhp() and
    /// wellRates() fields, depending on controls.  The
    /// perfRates() field is filled with zero, and perfPress()
    /// with -1e100.
    void base_init(const std::vector<Scalar>& cellPressures,
                   const std::vector<Scalar>& cellTemperatures,
                   const std::vector<Well>& wells_ecl,
                   const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
                   const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
                   const SummaryState& summary_state);

    void initSingleWell(const std::vector<Scalar>& cellPressures,
                        const std::vector<Scalar>& cellTemperatures,
                        const Well& well,
                        const std::vector<PerforationData<Scalar>>& well_perf_data,
                        const ParallelWellInfo<Scalar>& well_info,
                        const SummaryState& summary_state);

    void initSingleProducer(const Well& well,
                            const ParallelWellInfo<Scalar>& well_info,
                            Scalar pressure_first_connection,
                            const std::vector<PerforationData<Scalar>>& well_perf_data,
                            const SummaryState& summary_state);

    void initSingleInjector(const Well& well,
                            const ParallelWellInfo<Scalar>& well_info,
                            Scalar pressure_first_connection,
                            Scalar temperature_first_connection,
                            const std::vector<PerforationData<Scalar>>& well_perf_data,
                            const SummaryState& summary_state);

    static void calculateSegmentRatesBeforeSum(const ParallelWellInfo<Scalar>&      pw_info,
                                               const std::vector<std::vector<int>>& segment_inlets,
                                               const std::vector<std::vector<int>>& segment_perforations,
                                               const std::vector<Scalar>&           perforation_rates,
                                               const int                            np,
                                               const int                            segment,
                                               std::vector<Scalar>&                 segment_rates);

    void reportConnectionFactors(const std::size_t well_index,
                                 std::vector<data::Connection>& connections) const;

    void reportConnectionPressuresAndRates(const std::size_t well_index,
                                           const PhaseUsage& pu,
                                           std::vector<data::Connection>& connections) const;

    void reportConnectionFilterCake(const std::size_t well_index,
                                    std::vector<data::Connection>& connections) const;

    void reportFractureStatistics(const std::vector<ConnFracStatistics<Scalar>>& stats,
                                  std::vector<data::Connection>& connections) const;
};

} // namespace Opm

#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
