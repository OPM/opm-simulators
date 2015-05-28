/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS
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

namespace Opm
{
    template <class GridT>
    SimulatorFullyImplicitBlackoilPolymer<GridT>::
    SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
                                          const GridT& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const PolymerPropsAd& polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity,
                                          const bool has_disgas,
                                          const bool has_vapoil,
                                          const bool has_polymer,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
                                          const std::vector<double>& threshold_pressures_by_face)
    : BaseType(param,
               grid,
               geo,
               props,
               rock_comp_props,
               linsolver,
               gravity,
               has_disgas,
               has_vapoil,
               eclipse_state,
               output_writer,
               threshold_pressures_by_face)
        , polymer_props_(polymer_props)
        , has_polymer_(has_polymer)
        , deck_(deck)
    {
    }

    template <class GridT>
    auto SimulatorFullyImplicitBlackoilPolymer<GridT>::
    createSolver(const Wells* wells)
        -> std::shared_ptr<Solver>
    {
        typedef typename Traits::Model Model;
        typedef typename Model::ModelParameters ModelParams;
        ModelParams modelParams( BaseType::param_ );
        typedef NewtonSolver<Model> Solver;

        auto model = std::make_shared<Model>(modelParams,
                                             BaseType::grid_,
                                             BaseType::props_,
                                             BaseType::geo_,
                                             BaseType::rock_comp_props_,
                                             polymer_props_,
                                             wells,
                                             BaseType::solver_,
                                             BaseType::has_disgas_,
                                             BaseType::has_vapoil_,
                                             has_polymer_,
                                             BaseType::terminal_output_);

        if (!BaseType::threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(BaseType::threshold_pressures_by_face_);
        }

        typedef typename Solver::SolverParameters SolverParams;
        SolverParams solverParams( BaseType::param_ );
        return std::make_shared<Solver>(solverParams, model);
    }

    template <class GridT>
    void SimulatorFullyImplicitBlackoilPolymer<GridT>::
    handleAdditionalWellInflow(SimulatorTimer& timer,
                                    WellsManager& wells_manager,
                                    typename BaseType::WellState& well_state,
                                    const Wells* wells)
    {
        // compute polymer inflow
        std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
        if (deck_->hasKeyword("WPOLYMER")) {
            if (wells_manager.c_wells() == 0) {
                OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
            }
            polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, BaseType::eclipse_state_, *wells, Opm::UgGridHelpers::numCells(BaseType::grid_), timer.currentStepNum()));
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

} // namespace Opm
