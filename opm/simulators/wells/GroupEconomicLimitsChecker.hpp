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

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <map>
#include <string>

namespace Opm
{

class BlackoilWellModelGeneric;
class DeferredLogger;
class Group;
class WellState;
class WellTestState;

    class GroupEconomicLimitsChecker
    {
    public:
        GroupEconomicLimitsChecker(
            const BlackoilWellModelGeneric &well_model,
            WellTestState &well_test_state,
            const Group &group,
            const double simulation_time,
            const int report_step_idx,
            DeferredLogger &deferred_logger
        );
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
        std::string message_separator(const char sep_char='*', const size_t sep_length=110) const { return std::string(sep_length, sep_char); }

        static constexpr int NUM_PHASES = 3;
    private:
        void displayDebugMessage(const std::string &msg) const;
        void addPrintMessage(const std::string &msg, const double value, const double limit, const UnitSystem::measure measure);
        bool closeWellsRecursive(const Group& group, int level = 0);
        void throwNotImplementedError(const std::string &error) const;
        const BlackoilWellModelGeneric &well_model_;
        const Group &group_;
        const double simulation_time_;
        const int report_step_idx_;
        DeferredLogger &deferred_logger_;
        const std::string date_string_;
        const UnitSystem& unit_system_;
        const WellState &well_state_;
        WellTestState &well_test_state_;
        const Schedule &schedule_;
        GroupEconProductionLimits::GEconGroupProp gecon_props_;
        bool debug_ = true;
        double production_rates_[NUM_PHASES];
        std::map<int, BlackoilPhases::PhaseIndex> phase_idx_map_ = {
            {0, BlackoilPhases::Liquid},
            {1, BlackoilPhases::Vapour},
            {2, BlackoilPhases::Aqua}};
        std::map<BlackoilPhases::PhaseIndex, int> phase_idx_reverse_map_;
        std::string message_;
    };

} // namespace Opm

#endif // OPM_GROUP_ECONOMIC_LIMITS_CHECKER_HEADER_INCLUDED
