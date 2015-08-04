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

#ifndef OPM_SIMULATORFULLYIMPLICITEXTENDEDBLACKOIL_IMPL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITEXTENDEDBLACKOIL_IMPL_HEADER_INCLUDED

namespace Opm
{
    template <class GridT>
    SimulatorFullyImplicitExtendedBlackoil<GridT>::
    SimulatorFullyImplicitExtendedBlackoil(const parameter::ParameterGroup& param,
                                          const GridT& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const RockCompressibility* rock_comp_props,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity,
                                          const bool has_disgas,
                                          const bool has_vapoil,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
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
               output_writer,
               threshold_pressures_by_face)
    , has_solvent_(has_solvent)
    , deck_(deck)
    {
        if(deck->hasKeyword("MISCIBLE")) {
            std::cerr << "MISICIBLE keyword is present. Mixing is not currently supported" << std::endl;
        }
    }

    template <class GridT>
    auto SimulatorFullyImplicitExtendedBlackoil<GridT>::
    createSolver(const Wells* wells)
        -> std::unique_ptr<Solver>
    {
        typedef typename Traits::Model Model;


        auto model = std::unique_ptr<Model>(new Model(BaseType::model_param_,
                                                      BaseType::grid_,
                                                      BaseType::props_,
                                                      BaseType::geo_,
                                                      BaseType::rock_comp_props_,
                                                      wells,
                                                      BaseType::solver_,
                                                      BaseType::eclipse_state_,
                                                      BaseType::has_disgas_,
                                                      BaseType::has_vapoil_,
                                                      BaseType::terminal_output_,
                                                      has_solvent_));

        if (!BaseType::threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(BaseType::threshold_pressures_by_face_);
        }

        return std::unique_ptr<Solver>(new Solver(BaseType::solver_param_, std::move(model)));
    }

    template <class GridT>
    void SimulatorFullyImplicitExtendedBlackoil<GridT>::
    handleAdditionalWellInflow(SimulatorTimer& timer,
                                    WellsManager& wells_manager,
                                    typename BaseType::WellState& well_state,
                                    const Wells* wells)
    {
        // compute solvent inflow
        // TODO: Allow for changing WSOLVENT during simulation.
        if (deck_->hasKeyword("WSOLVENT")) {
            const int nw = wells->number_of_wells;
            std::vector<double> perfcells_fraction(wells->well_connpos[nw]);
            //perfcells_fraction.resize(wells->well_connpos[nw]);
            Opm::DeckKeywordConstPtr keyword = deck_->getKeyword("WSOLVENT");
            const int num_keywords = keyword->size();
            for (int i = 0; i < num_keywords; ++i) {
                int wix = 0;
                for (; wix < nw; ++wix) {
                    if (keyword->getRecord(i)->getItem("WELL")->getString(0) == wells->name[wix]) {
                        break;
                    }
                }
                if (wix == nw ) {
                    OPM_THROW(std::runtime_error, "Could not find a match for well "
                              << keyword->getRecord(i)->getItem("WELL")->getString(0)
                              << " from WSOLVENT.");
                }
                for (int j = wells->well_connpos[wix]; j < wells->well_connpos[wix+1]; ++j) {
                    perfcells_fraction[j] =
                        keyword->getRecord(i)->getItem("SOLVENT_FRACTION")->getSIDouble(0);
                }
            }
            well_state.solventFraction() = perfcells_fraction;
        }


    }


} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITEXTENDEDBLACKOIL_IMPL_HEADER_INCLUDED
