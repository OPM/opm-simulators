/*
  Copyright 2019 Norce.

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


#ifndef OPM_FRACTION_CALCULATOR_HEADER_INCLUDED
#define OPM_FRACTION_CALCULATOR_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>

#include <string>

namespace Opm {
class GroupState;
struct PhaseUsage;
class Schedule;
class WellState;
}

namespace Opm::WellGroupHelpers {

class FractionCalculator
{
public:
    FractionCalculator(const Schedule& schedule,
                       const WellState& well_state,
                       const GroupState& group_state,
                       const int report_step,
                       const GuideRate* guide_rate,
                       const GuideRateModel::Target target,
                       const PhaseUsage& pu,
                       const bool is_producer,
                       const Phase injection_phase);
    double fraction(const std::string& name,
                    const std::string& control_group_name,
                    const bool always_include_this);
    double localFraction(const std::string& name,
                         const std::string& always_included_child);

private:
    std::string parent(const std::string& name);
    double guideRateSum(const Group& group,
                        const std::string& always_included_child);
    double guideRate(const std::string& name,
                     const std::string& always_included_child);
    int groupControlledWells(const std::string& group_name,
                             const std::string& always_included_child);
    GuideRate::RateVector getGroupRateVector(const std::string& group_name);
    const Schedule& schedule_;
    const WellState& well_state_;
    const GroupState& group_state_;
    int report_step_;
    const GuideRate* guide_rate_;
    GuideRateModel::Target target_;
    const PhaseUsage& pu_;
    bool is_producer_;
    Phase injection_phase_;
};

} // namespace Opm::WellGroupHelpers

#endif // OPM_FRACTION_CALCULATOR_HEADER_INCLUDED
