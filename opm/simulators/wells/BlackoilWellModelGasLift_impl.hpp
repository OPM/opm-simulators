/*
  Copyright 2016 - 2019 SINTEF Digital, Mathematics & Cybernetics.
  Copyright 2016 - 2018 Equinor ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 Norce AS

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
#define OPM_BLACKOILWELLMODEL_GASLIFT_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_GASLIFT_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_BLACKOILWELLMODEL_GASLIFT_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelGasLift.hpp>
#endif
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

namespace Opm {

template<typename TypeTag>
bool
BlackoilWellModelGasLift<TypeTag>::
maybeDoGasLiftOptimize(const Simulator& simulator,
                       const std::vector<WellInterfacePtr>& well_container,
                       const std::map<std::string, Scalar>& node_pressures,
                       const bool updatePotentials,
                       WellStateType& wellState,
                       GroupState<Scalar>& groupState,
                       DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const auto& glo = simulator.vanguard().schedule().glo(simulator.episodeIndex());
    if(!glo.active()) {
        return false;
    }
    bool do_glift_optimization = false;
    int num_wells_changed = 0;
    const double simulation_time = simulator.time();
    const Scalar min_wait = glo.min_wait();
    // We only optimize if a min_wait time has past.
    // If all_newton is true we still want to optimize several times pr timestep
    // i.e. we also optimize if check simulation_time == last_glift_opt_time_
    // that is when the last_glift_opt_time is already updated with the current time step
    if (simulation_time == this->last_glift_opt_time_  ||
        simulation_time >= (this->last_glift_opt_time_ + min_wait))
    {
        do_glift_optimization = true;
        this->last_glift_opt_time_ = simulation_time;
    }

    // update potentials if changed by network (needed by the gaslift algorithm)
    if (updatePotentials) {
        updateWellPotentials(simulator, well_container, node_pressures, wellState, deferred_logger);
    }

    if (do_glift_optimization) {
        GLiftOptWells glift_wells;
        GLiftProdWells prod_wells;
        GLiftWellStateMap state_map;
        // NOTE: To make GasLiftGroupInfo (see below) independent of the TypeTag
        //  associated with *this (i.e. BlackoilWellModel<TypeTag>) we observe
        //  that GasLiftGroupInfo's only dependence on *this is that it needs to
        //  access the eclipse Wells in the well container (the eclipse Wells
        //  themselves are independent of the TypeTag).
        //  Hence, we extract them from the well container such that we can pass
        //  them to the GasLiftGroupInfo constructor.
        GLiftEclWells ecl_well_map;
        initGliftEclWellMap(well_container, ecl_well_map);
        GasLiftGroupInfo group_info {
            ecl_well_map,
            simulator.vanguard().schedule(),
            simulator.vanguard().summaryState(),
            simulator.episodeIndex(),
            simulator.model().newtonMethod().numIterations(),
            deferred_logger,
            wellState,
            groupState,
            simulator.vanguard().grid().comm(),
            this->glift_debug
        };
        group_info.initialize();

        gasLiftOptimizationStage1(simulator,
                                  well_container,
                                  wellState,
                                  groupState,
                                  prod_wells,
                                  glift_wells,
                                  group_info,
                                  state_map,
                                  deferred_logger);

        this->gasLiftOptimizationStage2(simulator.vanguard().gridView().comm(),
                                        simulator.vanguard().schedule(),
                                        simulator.vanguard().summaryState(),
                                        wellState,
                                        groupState,
                                        prod_wells,
                                        glift_wells,
                                        group_info,
                                        state_map,
                                        simulator.episodeIndex(),
                                        deferred_logger);

        if constexpr (glift_debug) {
            std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*> wc;
            wc.reserve(well_container.size());
            std::transform(well_container.begin(), well_container.end(),
                           std::back_inserter(wc),
                           [](const auto& w)
                           { return static_cast<WellInterfaceGeneric<Scalar, IndexTraits>*>(w.get()); });
            this->gliftDebugShowALQ(wc,
                                    wellState,
                                    deferred_logger);
        }
        num_wells_changed = glift_wells.size();
    }
    num_wells_changed = simulator.vanguard().gridView().comm().sum(num_wells_changed);
    return num_wells_changed > 0;
}

template<typename TypeTag>
void
BlackoilWellModelGasLift<TypeTag>::
gasLiftOptimizationStage1(const Simulator& simulator,
                          const std::vector<WellInterfacePtr>& well_container,
                          WellStateType& wellState,
                          GroupState<Scalar>& groupState,
                          GLiftProdWells& prod_wells,
                          GLiftOptWells &glift_wells,
                          GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                          GLiftWellStateMap& state_map,
                          DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    auto comm = simulator.vanguard().grid().comm();
    int num_procs = comm.size();
    // NOTE: Gas lift optimization stage 1 seems to be difficult
    //  to do in parallel since the wells are optimized on different
    //  processes and each process needs to know the current ALQ allocated
    //  to each group it is a memeber of in order to check group limits and avoid
    //  allocating more ALQ than necessary.  (Surplus ALQ is removed in
    //  stage 2). In stage1, as each well is adding ALQ, the current group ALQ needs
    //  to be communicated to the other processes.  But there is no common
    //  synchronization point that all process will reach in the
    //  runOptimizeLoop_() in GasLiftSingleWell.cpp.
    //
    //  TODO: Maybe a better solution could be invented by distributing
    //    wells according to certain parent groups. Then updated group rates
    //    might not have to be communicated to the other processors.

    //  Currently, the best option seems to be to run this part sequentially
    //    (not in parallel).
    //
    //  TODO: The simplest approach seems to be if a) one process could take
    //    ownership of all the wells (the union of all the wells in the
    //    well_container_ of each process) then this process could do the
    //    optimization, while the other processes could wait for it to
    //    finish (e.g. comm.barrier()), or alternatively, b) if all
    //    processes could take ownership of all the wells.  Then there
    //    would be no need for synchronization here..
    //
    for (int i = 0; i< num_procs; i++) {
        int num_rates_to_sync = 0;  // communication variable
        GLiftSyncGroups groups_to_sync;
        if (comm.rank() ==  i) {
            // Run stage1: Optimize single wells while also checking group limits
            for (const auto& well : well_container) {
                // NOTE: Only the wells in "group_info" needs to be optimized
                if (group_info.hasWell(well->name())) {
                    gasLiftOptimizationStage1SingleWell(well.get(),
                                                        simulator,
                                                        wellState,
                                                        groupState,
                                                        prod_wells,
                                                        glift_wells,
                                                        group_info,
                                                        state_map,
                                                        groups_to_sync,
                                                        deferred_logger);
                }
            }
            num_rates_to_sync = groups_to_sync.size();
        }
        {
            OPM_TIMEBLOCK(WaitForGasLiftSyncGroups);
            num_rates_to_sync = comm.sum(num_rates_to_sync);
        }
        if (num_rates_to_sync > 0) {
            OPM_TIMEBLOCK(GasLiftSyncGroups);
            std::vector<int> group_indexes;
            group_indexes.reserve(num_rates_to_sync);
            std::vector<Scalar> group_alq_rates;
            group_alq_rates.reserve(num_rates_to_sync);
            std::vector<Scalar> group_oil_rates;
            group_oil_rates.reserve(num_rates_to_sync);
            std::vector<Scalar> group_gas_rates;
            group_gas_rates.reserve(num_rates_to_sync);
            std::vector<Scalar> group_water_rates;
            group_water_rates.reserve(num_rates_to_sync);
            if (comm.rank() == i) {
                for (auto idx : groups_to_sync) {
                    auto [oil_rate, gas_rate, water_rate, alq] = group_info.getRates(idx);
                    group_indexes.push_back(idx);
                    group_oil_rates.push_back(oil_rate);
                    group_gas_rates.push_back(gas_rate);
                    group_water_rates.push_back(water_rate);
                    group_alq_rates.push_back(alq);
                }
            } else {
                group_indexes.resize(num_rates_to_sync);
                group_oil_rates.resize(num_rates_to_sync);
                group_gas_rates.resize(num_rates_to_sync);
                group_water_rates.resize(num_rates_to_sync);
                group_alq_rates.resize(num_rates_to_sync);
            }
#if HAVE_MPI
            Parallel::MpiSerializer ser(comm);
            ser.broadcast(Parallel::RootRank{i}, group_indexes, group_oil_rates,
                          group_gas_rates, group_water_rates, group_alq_rates);
#endif
            if (comm.rank() != i) {
                for (int j = 0; j < num_rates_to_sync; ++j) {
                    group_info.updateRate(group_indexes[j],
                                          group_oil_rates[j],
                                          group_gas_rates[j],
                                          group_water_rates[j],
                                          group_alq_rates[j]);
                }
            }
            if constexpr (glift_debug) {
                int counter = 0;
                if (comm.rank() == i) {
                    counter = wellState.gliftGetDebugCounter();
                }
                counter = comm.sum(counter);
                if (comm.rank() != i) {
                    wellState.gliftSetDebugCounter(counter);
                }
            }
        }
    }
}

// NOTE: this method cannot be const since it passes this->wellState()
//   (see below) to the GasLiftSingleWell constructor which accepts WellState
//   as a non-const reference..
template<typename TypeTag>
void
BlackoilWellModelGasLift<TypeTag>::
gasLiftOptimizationStage1SingleWell(WellInterface<TypeTag>* well,
                                    const Simulator& simulator,
                                    WellStateType& wellState,
                                    GroupState<Scalar>& groupState,
                                    GLiftProdWells& prod_wells,
                                    GLiftOptWells& glift_wells,
                                    GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                                    GLiftWellStateMap& state_map,
                                    GLiftSyncGroups& sync_groups,
                                    DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const auto& summary_state = simulator.vanguard().summaryState();
    auto glift = std::make_unique<GasLiftSingleWell<TypeTag>>(*well,
                                                              simulator,
                                                              summary_state,
                                                              deferred_logger,
                                                              wellState,
                                                              groupState,
                                                              group_info,
                                                              sync_groups,
                                                              simulator.vanguard().gridView().comm(),
                                                              this->glift_debug);
    auto state = glift->runOptimize(simulator.model().newtonMethod().numIterations());
    if (state) {
        state_map.emplace(well->name(), std::move(state));
        glift_wells.emplace(well->name(), std::move(glift));
        return;
    }
    prod_wells.insert({well->name(), well});
}

template<typename TypeTag>
void
BlackoilWellModelGasLift<TypeTag>::
initGliftEclWellMap(const std::vector<WellInterfacePtr>& well_container,
                    GLiftEclWells& ecl_well_map)
{
    for (const auto& well : well_container) {
        ecl_well_map.try_emplace(well->name(), &well->wellEcl(), well->indexOfWell());
    }
}

template<typename TypeTag>
void
BlackoilWellModelGasLift<TypeTag>::
updateWellPotentials(const Simulator& simulator,
                     const std::vector<WellInterfacePtr>& well_container,
                     const std::map<std::string, Scalar>& node_pressures,
                     WellStateType& wellState,
                     DeferredLogger& deferred_logger)
{
    auto well_state_copy = wellState;
    const int np = wellState.numPhases();
    std::string exc_msg;
    auto exc_type = ExceptionType::NONE;
    for (const auto& well : well_container) {
        if (well->isInjector() || !well->wellEcl().predictionMode())
            continue;

        const auto it = node_pressures.find(well->wellEcl().groupName());
        if (it != node_pressures.end()) {
            std::string cur_exc_msg;
            auto cur_exc_type = ExceptionType::NONE;
            try {
                std::vector<Scalar> potentials;
                well->computeWellPotentials(simulator, well_state_copy, potentials, deferred_logger);
                auto& ws = wellState.well(well->indexOfWell());
                for (int p = 0; p < np; ++p) {
                    // make sure the potentials are positive
                    ws.well_potentials[p] = std::max(Scalar{0.0}, potentials[p]);
                }
            }
            // catch all possible exception and store type and message.
            OPM_PARALLEL_CATCH_CLAUSE(cur_exc_type, cur_exc_msg);
            if (cur_exc_type != ExceptionType::NONE) {
                exc_msg += fmt::format("\nFor well {}: {}", well->name(), cur_exc_msg);
            }
            exc_type = std::max(exc_type, cur_exc_type);
        }
    }
    if (exc_type != ExceptionType::NONE) {
        const std::string msg = "Updating well potentials after network balancing failed. Continue with current potentials";
        deferred_logger.warning("WELL_POT_SOLVE_AFTER_NETWORK_FAILED", msg + exc_msg);
    }
}

} // namespace Opm
