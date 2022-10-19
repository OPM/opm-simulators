/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELLINTERFACE_EVAL_HEADER_INCLUDED
#define OPM_WELLINTERFACE_EVAL_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <functional>

namespace Opm
{

class DeferredLogger;
class Group;
class GroupState;
class Schedule;
class SummaryState;
template<class FluidSystem> class WellInterfaceFluidSystem;
class WellState;

template<class FluidSystem>
class WellInterfaceEval {
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

public:
    template <class EvalWell>
    EvalWell calculateBhpFromThp(const WellState& well_state,
                                 const std::vector<EvalWell>& rates,
                                 const Well& well,
                                 const SummaryState& summaryState,
                                 const double rho,
                                 DeferredLogger& deferred_logger) const;

    template<class EvalWell, class BhpFromThpFunc>
    void assembleControlEqProd(const WellState& well_state,
                               const GroupState& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const Well::ProductionControls& controls,
                               const EvalWell& bhp,
                               const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                               BhpFromThpFunc bhp_from_thp,
                               EvalWell& control_eq,
                               DeferredLogger& deferred_logger) const
    {
        std::function<EvalWell()> eval = [&bhp_from_thp]() { return bhp_from_thp(); };
        assembleControlEqProd_(well_state,
                               group_state,
                               schedule,
                               summaryState,
                               controls,
                               bhp,
                               rates,
                               eval,
                               control_eq,
                               deferred_logger);
    }

    template<class EvalWell>
    void assembleControlEqProd_(const WellState& well_state,
                                const GroupState& group_state,
                                const Schedule& schedule,
                                const SummaryState& summaryState,
                                const Well::ProductionControls& controls,
                                const EvalWell& bhp,
                                const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                const std::function<EvalWell()>& bhp_from_thp,
                                EvalWell& control_eq,
                                DeferredLogger& deferred_logger) const;

    template<class EvalWell, class BhpFromThpFunc>
    void assembleControlEqInj(const WellState& well_state,
                              const GroupState& group_state,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              const Well::InjectionControls& controls,
                              const EvalWell& bhp,
                              const EvalWell& injection_rate,
                              BhpFromThpFunc bhp_from_thp,
                              EvalWell& control_eq,
                              DeferredLogger& deferred_logger) const
    {
        std::function<EvalWell()> eval = [&bhp_from_thp]() { return bhp_from_thp(); };
        assembleControlEqInj_(well_state,
                              group_state,
                              schedule,
                              summaryState,
                              controls,
                              bhp,
                              injection_rate,
                              eval,
                              control_eq,
                              deferred_logger);
    }

    template<class EvalWell>
    void assembleControlEqInj_(const WellState& well_state,
                               const GroupState& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const Well::InjectionControls& controls,
                               const EvalWell& bhp,
                               const EvalWell& injection_rate,
                               const std::function<EvalWell()>& bhp_from_thp,
                               EvalWell& control_eq,
                               DeferredLogger& deferred_logger) const;

protected:
    WellInterfaceEval(const WellInterfaceFluidSystem<FluidSystem>& baseif);

    const WellInterfaceFluidSystem<FluidSystem>& baseif_;
};

}

#endif // OPM_WELLINTERFACE_EVAL_HEADER_INCLUDED
