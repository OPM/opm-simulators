/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#ifndef OPM_BLACKOILWELLMODEL_GASLIFT_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_GASLIFT_HEADER_INCLUDED

#include "opm/models/utils/basicproperties.hh"
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>

#include <memory>
#include <map>
#include <string>

namespace Opm {

class DeferredLogger;
template<class Scalar> class GroupState;
template<typename Scalar, typename IndexTraits> class WellState;
template<class TypeTag> class WellInterface;

template<typename Scalar, typename IndexTraits>
class BlackoilWellModelGasLiftGeneric
{
public:
    using GLiftOptWells = std::map<std::string, std::unique_ptr<GasLiftSingleWellGeneric<Scalar, IndexTraits>>>;
    using GLiftProdWells = std::map<std::string, const WellInterfaceGeneric<Scalar, IndexTraits>*>;
    using GLiftWellStateMap = std::map<std::string, std::unique_ptr<GasLiftWellState<Scalar>>>;
    using GLiftEclWells = typename GasLiftGroupInfo<Scalar, IndexTraits>::GLiftEclWells;
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::GLiftSyncGroups;

    explicit BlackoilWellModelGasLiftGeneric(bool terminal_output)
        : terminal_output_(terminal_output)
    {}

    static constexpr bool glift_debug = false;

    void gliftDebug(const std::string& msg,
                    DeferredLogger& deferred_logger) const;

    bool terminalOutput() const { return terminal_output_; }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(last_glift_opt_time_);
    }

    bool operator==(const BlackoilWellModelGasLiftGeneric& that) const
    { return this->last_glift_opt_time_ == that.last_glift_opt_time_; }

protected:
    void gliftDebugShowALQ(const std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*>& well_container,
                           const WellState<Scalar, IndexTraits>& wellState,
                           DeferredLogger& deferred_logger);

    void gasLiftOptimizationStage2(const Parallel::Communication& comm,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   WellState<Scalar, IndexTraits>& wellState,
                                   GroupState<Scalar>& groupState,
                                   GLiftProdWells& prod_wells,
                                   GLiftOptWells& glift_wells,
                                   GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                                   GLiftWellStateMap& map,
                                   const int episodeIndex,
                                   DeferredLogger& deferred_logger);

    bool terminal_output_;
    double last_glift_opt_time_ = -1.0;
};

/// Class for handling the gaslift in the blackoil well model.
template<typename TypeTag>
class BlackoilWellModelGasLift :
    public BlackoilWellModelGasLiftGeneric<GetPropType<TypeTag, Properties::Scalar>,
                                           typename GetPropType<TypeTag, Properties::FluidSystem>::IndexTraitsType>
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using Base = BlackoilWellModelGasLiftGeneric<Scalar, IndexTraits>;
    using Base::glift_debug;
    using GLiftEclWells = typename GasLiftGroupInfo<Scalar, IndexTraits>::GLiftEclWells;
    using GLiftOptWells = typename Base::GLiftOptWells;
    using GLiftProdWells = typename Base::GLiftProdWells;
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::GLiftSyncGroups;
    using GLiftWellStateMap =  typename Base::GLiftWellStateMap;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag>>;
    using WellStateType = WellState<Scalar, IndexTraits>;

    explicit BlackoilWellModelGasLift(bool terminal_output)
        : Base(terminal_output)
    {}

    static void initGliftEclWellMap(const std::vector<WellInterfacePtr>& well_container,
                                    GLiftEclWells& ecl_well_map);

    bool maybeDoGasLiftOptimize(const Simulator& simulator,
                                const std::vector<WellInterfacePtr>& well_container,
                                const std::map<std::string, Scalar>& node_pressures,
                                const bool updatePotentials,
                                WellStateType& wellState,
                                GroupState<Scalar>& groupState,
                                DeferredLogger& deferred_logger);

private:
    void gasLiftOptimizationStage1(const Simulator& simulator,
                                   const std::vector<WellInterfacePtr>& well_container,
                                   WellStateType& wellState,
                                   GroupState<Scalar>& groupState,
                                   GLiftProdWells& prod_wells,
                                   GLiftOptWells& glift_wells,
                                   GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                                   GLiftWellStateMap& state_map,
                                   DeferredLogger& deferred_logger);

    // cannot be const since it accesses the non-const WellState
    void gasLiftOptimizationStage1SingleWell(WellInterface<TypeTag>* well,
                                             const Simulator& simulator,
                                             WellStateType& wellState,
                                             GroupState<Scalar>& groupState,
                                             GLiftProdWells& prod_wells,
                                             GLiftOptWells& glift_wells,
                                             GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                                             GLiftWellStateMap& state_map,
                                             GLiftSyncGroups& groups_to_sync,
                                             DeferredLogger& deferred_logger);

    void updateWellPotentials(const Simulator& simulator,
                              const std::vector<WellInterfacePtr>& well_container,
                              const std::map<std::string, Scalar>& node_pressures,
                              WellState<Scalar, IndexTraits>& wellState,
                              DeferredLogger& deferred_logger);

};

} // namespace Opm

#include "BlackoilWellModelGasLift_impl.hpp"

#endif
