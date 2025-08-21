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
template<class Scalar> class GroupState;
class Schedule;
template<typename Scalar, typename IndexTraits> class WellState;
}

namespace Opm::WGHelpers {

template<typename Scalar, typename IndexTraits>
class FractionCalculator
{
public:
    FractionCalculator(const Schedule& schedule,
                       const WellState<Scalar, IndexTraits>& well_state,
                       const GroupState<Scalar>& group_state,
                       const SummaryState& summary_state,
                       const int report_step,
                       const GuideRate* guide_rate,
                       const GuideRateModel::Target target,
                       const bool is_producer,
                       const Phase injection_phase);
    Scalar fraction(const std::string& name,
                    const std::string& control_group_name,
                    const bool always_include_this);
    Scalar localFraction(const std::string& name,
                         const std::string& always_included_child);

private:
    std::string parent(const std::string& name);

    // returns the sum of the guiderates of the given group
    // and the number of sub-groups/wells that contributed to the sum
    std::pair<Scalar,int> guideRateSum(const Group& group,
                        const std::string& always_included_child,
                        const bool always_use_potentials);
    Scalar guideRate(const std::string& name,
                     const std::string& always_included_child,
                     const bool always_use_potentials);
    int groupControlledWells(const std::string& group_name,
                             const std::string& always_included_child);
    GuideRate::RateVector getGroupRateVector(const std::string& group_name);
    const Schedule& schedule_;
    const WellState<Scalar, IndexTraits>& well_state_;
    const GroupState<Scalar>& group_state_;
    const SummaryState& summary_state_;
    int report_step_;
    const GuideRate* guide_rate_;
    GuideRateModel::Target target_;
    bool is_producer_;
    Phase injection_phase_;
};

} // namespace Opm::WGHelpers

#endif // OPM_FRACTION_CALCULATOR_HEADER_INCLUDED
