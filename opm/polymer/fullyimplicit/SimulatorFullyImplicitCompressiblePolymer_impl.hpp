/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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
#ifndef OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED

namespace Opm
{

/// Class collecting all necessary components for a two-phase simulation.
template <class GridT>
SimulatorFullyImplicitCompressiblePolymer<GridT>::
SimulatorFullyImplicitCompressiblePolymer(const ParameterGroup& param,
                                          const GridT& grid,
                                          DerivedGeology& geo,
                                          BlackoilPropsAdFromDeck& props,
                                          const PolymerPropsAd&    polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          const Deck& deck,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity)
: BaseType(param,
           grid,
           geo,
           props,
           rock_comp_props,
           linsolver,
           gravity,
           /*disgas=*/false,
           /*vapoil=*/false,
           eclipse_state,
           output_writer,
           /*threshold_pressures_by_face=*/std::vector<double>(),
           // names of deactivated wells in parallel run
           std::unordered_set<std::string>())
    , deck_(deck)
    , polymer_props_(polymer_props)

{
}

template <class GridT>
auto SimulatorFullyImplicitCompressiblePolymer<GridT>::
createSolver(const WellModel& well_model)
    -> std::unique_ptr<Solver>
{
    return std::unique_ptr<Solver>(new Solver(BaseType::grid_,
                                              BaseType::props_,
                                              BaseType::geo_,
                                              BaseType::rock_comp_props_,
                                              polymer_props_,
                                              // *wells,
                                              // TODO: it is resulted from refactoring of other simulators.
                                              well_model.wells(),
                                              BaseType::solver_));
}

template <class GridT>
void SimulatorFullyImplicitCompressiblePolymer<GridT>::
handleAdditionalWellInflow(SimulatorTimer& timer,
                           WellsManager& wells_manager,
                           typename BaseType::WellState& well_state,
                           const Wells* wells)
{
    // compute polymer inflow
    std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
    if (deck_.hasKeyword("WPOLYMER")) {
        if (wells_manager.c_wells() == 0) {
            OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
        }
        polymer_inflow_ptr.reset(new PolymerInflowFromDeck( *BaseType::eclipse_state_, *wells, Opm::UgGridHelpers::numCells(BaseType::grid_), timer.currentStepNum()));
    } else {
        polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                        1.0*Opm::unit::day,
                                                        0.0));
    }
    std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(BaseType::grid_));
    polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(),
                                        timer.simulationTimeElapsed() + timer.currentStepLength(),
                                        polymer_inflow_c);
    well_state.polymerInflow() = polymer_inflow_c;
}





template <class GridT>
void
SimulatorFullyImplicitCompressiblePolymer<GridT>::
updateListEconLimited(const std::unique_ptr<Solver>& /*solver*/,
                      const Schedule& /*schedule*/,
                      const int /*current_step*/,
                      const Wells* /*wells*/,
                      const WellState& /*well_state*/,
                      DynamicListEconLimited& /*list_econ_limited*/) const
{

}


} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED
