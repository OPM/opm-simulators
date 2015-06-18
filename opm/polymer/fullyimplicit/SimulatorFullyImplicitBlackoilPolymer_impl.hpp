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
                                          const bool has_plyshlog,
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
        , has_plyshlog_(has_plyshlog)
        , deck_(deck)
    {
    }

    template <class GridT>
    auto SimulatorFullyImplicitBlackoilPolymer<GridT>::
    createSolver(const Wells* wells)
        -> std::unique_ptr<Solver>
    {
        typedef typename Traits::Model Model;


        auto model = std::unique_ptr<Model>(new Model(BaseType::model_param_,
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
                                                      has_plyshlog_,
                                                      wells_rep_radius_,
                                                      wells_perf_length_,
                                                      BaseType::terminal_output_));

        if (!BaseType::threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(BaseType::threshold_pressures_by_face_);
        }

        return std::unique_ptr<Solver>(new Solver(BaseType::solver_param_, std::move(model)));
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

        computeRepRadiusPerfLength(BaseType::eclipse_state_, timer.currentStepNum(), BaseType::grid_, wells_rep_radius_, wells_perf_length_);
    }


    template <class GridT>
    void SimulatorFullyImplicitBlackoilPolymer<GridT>::
    setupCompressedToCartesian(const int* global_cell, int number_of_cells,
                               std::map<int,int>& cartesian_to_compressed )
    {
        if (global_cell) {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
            }
        }
        else {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(i, i));
            }
        }

    }


    template <class GridT>
    void SimulatorFullyImplicitBlackoilPolymer<GridT>::
    computeRepRadiusPerfLength(const Opm::EclipseStateConstPtr eclipseState,
                               const size_t                    timeStep,
                               const GridT&                    grid,
                               std::vector<double>&            wells_rep_radius,
                               std::vector<double>&            wells_perf_length)
    {

        // TODO, the function does not work for parallel running
        // to be fixed later.
        int number_of_cells = Opm::UgGridHelpers::numCells(grid);
        const int* global_cell = Opm::UgGridHelpers::globalCell(grid);
        const int* cart_dims = Opm::UgGridHelpers::cartDims(grid);
        auto cell_to_faces = Opm::UgGridHelpers::cell2Faces(grid);
        auto begin_face_centroids = Opm::UgGridHelpers::beginFaceCentroids(grid);

        if (eclipseState->getSchedule()->numWells() == 0) {
            OPM_MESSAGE("No wells specified in Schedule section, "
                        "initializing no wells");
            return;
        }

        const size_t n_perf = wells_rep_radius.size();

        wells_rep_radius.clear();
        wells_perf_length.clear();

        wells_rep_radius.reserve(n_perf);
        wells_perf_length.reserve(n_perf);

        std::map<int,int> cartesian_to_compressed;

        setupCompressedToCartesian(global_cell, number_of_cells,
                                   cartesian_to_compressed);

        ScheduleConstPtr          schedule = eclipseState->getSchedule();
        std::vector<WellConstPtr> wells    = schedule->getWells(timeStep);

        int well_index = 0;

        for (auto wellIter= wells.begin(); wellIter != wells.end(); ++wellIter) {
             WellConstPtr well = (*wellIter);

             if (well->getStatus(timeStep) == WellCommon::SHUT) {
                 continue;
             }
             {   // COMPDAT handling
                 CompletionSetConstPtr completionSet = well->getCompletions(timeStep);
                 for (size_t c=0; c<completionSet->size(); c++) {
                     CompletionConstPtr completion = completionSet->get(c);
                     if (completion->getState() == WellCompletion::OPEN) {
                         int i = completion->getI();
                         int j = completion->getJ();
                         int k = completion->getK();

                         const int* cpgdim = cart_dims;
                         int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                         std::map<int, int>::const_iterator cgit = cartesian_to_compressed.find(cart_grid_indx);
                         if (cgit == cartesian_to_compressed.end()) {
                             OPM_THROW(std::runtime_error, "Cell with i,j,k indices " << i << ' ' << j << ' '
                                       << k << " not found in grid (well = " << well->name() << ')');
                         }
                         int cell = cgit->second;

                         {
                             double radius = 0.5*completion->getDiameter();
                             if (radius <= 0.0) {
                                 radius = 0.5*unit::feet;
                                 OPM_MESSAGE("**** Warning: Well bore internal radius set to " << radius);
                             }

                             const std::array<double, 3> cubical =
                             WellsManagerDetail::getCubeDim<3>(cell_to_faces, begin_face_centroids, cell);

                             WellCompletion::DirectionEnum direction = completion->getDirection();

                             double re; // area equivalent radius of the grid block
                             double perf_length; // the length of the well perforation

                             switch (direction) {
                                 case Opm::WellCompletion::DirectionEnum::X:
                                     re = std::sqrt(cubical[1] * cubical[2] / M_PI);
                                     perf_length = cubical[0];
                                     break;
                                 case Opm::WellCompletion::DirectionEnum::Y:
                                     re = std::sqrt(cubical[0] * cubical[2] / M_PI);
                                     perf_length = cubical[1];
                                     break;
                                 case Opm::WellCompletion::DirectionEnum::Z:
                                     re = std::sqrt(cubical[0] * cubical[1] / M_PI);
                                     perf_length = cubical[2];
                                     break;
                                 default:
                                     OPM_THROW(std::runtime_error, " Dirtecion of well is not supported ");
                             }

                             double repR = std::sqrt(re * radius);
                             wells_rep_radius.push_back(repR);
                             wells_perf_length.push_back(perf_length);
                         }
                     } else {
                         if (completion->getState() != WellCompletion::SHUT) {
                             OPM_THROW(std::runtime_error, "Completion state: " << WellCompletion::StateEnum2String( completion->getState() ) << " not handled");
                         }
                     }

                 }
            }
            well_index++;
        }
    }

} // namespace Opm
