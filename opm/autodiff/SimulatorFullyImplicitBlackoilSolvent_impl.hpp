/*
  Copyright 2015 IRIS AS

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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_IMPL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_IMPL_HEADER_INCLUDED

namespace Opm
{
    template <class GridT>
    SimulatorFullyImplicitBlackoilSolvent<GridT>::
    SimulatorFullyImplicitBlackoilSolvent(const ParameterGroup& param,
                                          const GridT& grid,
                                          DerivedGeology& geo,
                                          BlackoilPropsAdFromDeck& props,
                                          const SolventPropsAdFromDeck& solvent_props,
                                          const RockCompressibility* rock_comp_props,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity,
                                          const bool has_disgas,
                                          const bool has_vapoil,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          std::shared_ptr<Schedule> schedule,
                                          std::shared_ptr<SummaryConfig> summary_config,
                                          BlackoilOutputWriter& output_writer,
                                          std::shared_ptr< Deck > deck,
                                          const std::vector<double>& threshold_pressures_by_face,
                                          const bool has_solvent)
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
               schedule,
               summary_config,
               output_writer,
               threshold_pressures_by_face,
               // names of deactivated wells in parallel run
               std::unordered_set<std::string>())
    , has_solvent_(has_solvent)
    , deck_(deck)
    , solvent_props_(solvent_props)
    , is_miscible_(false)
    {
        if(deck->hasKeyword("MISCIBLE")) {
            is_miscible_ = true;
        }
    }

    template <class GridT>
    auto SimulatorFullyImplicitBlackoilSolvent<GridT>::
    createSolver(const WellModel& well_model)
        -> std::unique_ptr<Solver>
    {
        typedef typename Traits::Model Model;


        auto model = std::unique_ptr<Model>(new Model(BaseType::model_param_,
                                                      BaseType::grid_,
                                                      BaseType::props_,
                                                      BaseType::geo_,
                                                      BaseType::rock_comp_props_,
                                                      solvent_props_,
                                                      well_model,
                                                      BaseType::solver_,
                                                      BaseType::eclipse_state_,
                                                      BaseType::schedule_,
                                                      BaseType::summary_config_,
                                                      BaseType::has_disgas_,
                                                      BaseType::has_vapoil_,
                                                      BaseType::terminal_output_,
                                                      has_solvent_,
                                                      is_miscible_));

        if (!BaseType::threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(BaseType::threshold_pressures_by_face_);
        }

        return std::unique_ptr<Solver>(new Solver(BaseType::solver_param_, std::move(model)));
    }

    template <class GridT>
    void SimulatorFullyImplicitBlackoilSolvent<GridT>::
    handleAdditionalWellInflow(SimulatorTimer& timer,
			       WellsManager& /*wells_manager*/,
			       typename BaseType::WellState& well_state,
			       const Wells* wells)
    {
        // compute solvent inflow
        const int nw = wells->number_of_wells;
        std::vector<double> perfcells_fraction(wells->well_connpos[nw], 0.0);

        size_t currentStep = timer.currentStepNum();

        for (const auto&  well_solvent : BaseType::schedule_->getWells( currentStep )) {
            if (well_solvent->getStatus( currentStep ) == WellCommon::SHUT) {
                continue;
            }

            WellInjectionProperties injection = well_solvent->getInjectionProperties(currentStep);
            if (injection.injectorType == WellInjector::GAS) {

                double solventFraction = well_solvent->getSolventFraction(currentStep);

                // Find the solvent well in the well list and add properties to it
                int wix = 0;
                for (; wix < nw; ++wix) {
                    if (well_solvent->name() == wells->name[wix]) {
                        break;
                    }
                }
                if (wix == wells->number_of_wells) {
                    OPM_THROW(std::runtime_error, "Could not find a match for well "
                              << well_solvent->name()
                              << " from WSOLVENT.");
                }
                for (int j = wells->well_connpos[wix]; j < wells->well_connpos[wix+1]; ++j) {
                    perfcells_fraction[j] = solventFraction;
                }
            }
        }
        well_state.solventFraction() = perfcells_fraction;
    }

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_IMPL_HEADER_INCLUDED
