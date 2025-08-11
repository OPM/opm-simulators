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

#ifndef OPM_GROUP_ECONOMIC_LIMITS_CHECKER_HEADER_INCLUDED
#define OPM_GROUP_ECONOMIC_LIMITS_CHECKER_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <array>
#include <map>
#include <string>

namespace Opm
{

template<typename FluidSystem, typename Indices> class BlackoilWellModelGeneric;
class DeferredLogger;
class Group;
template<typename FluidSystem, typename Indices> class WellState;
class WellTestState;

template<typename FluidSystem, typename Indices>
class GroupEconomicLimitsChecker
{
public:
    using Scalar = typename FluidSystem::Scalar;

    GroupEconomicLimitsChecker(const BlackoilWellModelGeneric<FluidSystem, Indices>& well_model,
                               WellTestState& well_test_state,
                               const Group& group,
                               const double simulation_time,
                               const int report_step_idx,
                               DeferredLogger& deferred_logger);
    void closeWells();
    bool minGasRate();
    bool minOilRate();
    bool waterCut();
    bool GOR();
    bool WGR();
    void doWorkOver();
    bool endRun();
    int numProducersOpenInitially();
    int numProducersOpen();
    void activateEndRun();
    std::string message_separator(const char sep_char = '*',
                                  const size_t sep_length = 110) const
    { return std::string(sep_length, sep_char); }

    static constexpr int NUM_PHASES = 3;

private:
    void displayDebugMessage(const std::string& msg) const;
    void addPrintMessage(const std::string& msg,
                         const Scalar value,
                         const Scalar limit,
                         const UnitSystem::measure measure);
    bool closeWellsRecursive(const Group& group, int level = 0);
    void throwNotImplementedError(const std::string& error) const;

    const BlackoilWellModelGeneric<FluidSystem, Indices>& well_model_;
    const Group& group_;
    const double simulation_time_;
    const int report_step_idx_;
    DeferredLogger& deferred_logger_;
    const std::string date_string_;
    const UnitSystem& unit_system_;
    const WellState<FluidSystem, Indices>& well_state_;
    WellTestState& well_test_state_;
    const Schedule& schedule_;
    GroupEconProductionLimits::GEconGroupProp gecon_props_;
    bool debug_ = true;
    std::array<Scalar,NUM_PHASES> production_rates_;
    std::map<int, unsigned> phase_idx_map_ = {
        {0, FluidSystem::oilPhaseIdx},
        {1, FluidSystem::gasPhaseIdx},
        {2, FluidSystem::waterPhaseIdx}
    };
    std::map<unsigned, int> phase_idx_reverse_map_;
    std::string message_;
};

} // namespace Opm

#endif // OPM_GROUP_ECONOMIC_LIMITS_CHECKER_HEADER_INCLUDED
