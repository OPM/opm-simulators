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


#ifndef OPM_WELL_ASSEMBLE_HEADER_INCLUDED
#define OPM_WELL_ASSEMBLE_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

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
struct WellInjectionControls;
struct WellProductionControls;

template<class FluidSystem>
class WellAssemble {
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

public:
    WellAssemble(const WellInterfaceFluidSystem<FluidSystem>& well);

    template<class EvalWell>
    void assembleControlEqProd(const WellState& well_state,
                               const GroupState& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const WellProductionControls& controls,
                               const EvalWell& bhp,
                               const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                               const std::function<EvalWell()>& bhp_from_thp,
                               EvalWell& control_eq,
                               DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    void assembleControlEqInj(const WellState& well_state,
                              const GroupState& group_state,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              const WellInjectionControls& controls,
                              const EvalWell& bhp,
                              const EvalWell& injection_rate,
                              const std::function<EvalWell()>& bhp_from_thp,
                              EvalWell& control_eq,
                              DeferredLogger& deferred_logger) const;

private:
    const WellInterfaceFluidSystem<FluidSystem>& well_;
};



}

#endif // OPM_WELL_ASSEMBLE_HEADER_INCLUDED
