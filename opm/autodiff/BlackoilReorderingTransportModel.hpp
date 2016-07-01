/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED
#define OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/core/grid.h>
#include <opm/autodiff/DebugTimeReport.hpp>
#include <opm/core/transport/reorder/reordersequence.h>

namespace Opm {

    /// A model implementation for the transport equation in three-phase black oil.
    template<class Grid, class WellModel>
    class BlackoilReorderingTransportModel
        : public BlackoilModelBase<Grid, WellModel, BlackoilReorderingTransportModel<Grid, WellModel> >
    {
    public:
        typedef BlackoilModelBase<Grid, WellModel, BlackoilReorderingTransportModel<Grid, WellModel> > Base;
        friend Base;

        typedef typename Base::ReservoirState ReservoirState;
        typedef typename Base::WellState WellState;
        typedef typename Base::SolutionState SolutionState;
        typedef typename Base::V V;


        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells_arg        well structure
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilReorderingTransportModel(const typename Base::ModelParameters&   param,
                               const Grid&                             grid,
                               const BlackoilPropsAdInterface&         fluid,
                               const DerivedGeology&                   geo,
                               const RockCompressibility*              rock_comp_props,
                               const StandardWells&                    std_wells,
                               const NewtonIterationBlackoilInterface& linsolver,
                               Opm::EclipseStateConstPtr               eclState,
                               const bool                              has_disgas,
                               const bool                              has_vapoil,
                               const bool                              terminal_output)
            : Base(param, grid, fluid, geo, rock_comp_props, std_wells, linsolver,
                   eclState, has_disgas, has_vapoil, terminal_output)
            , reservoir_state0_(0, 0, 0)
            , well_state0_()
        {
        }





        void prepareStep(const double dt,
                         const ReservoirState& reservoir_state,
                         const WellState& well_state)
        {
            Base::prepareStep(dt, reservoir_state, well_state);
            Base::param_.solve_welleq_initially_ = false;
            reservoir_state0_ = reservoir_state;
            well_state0_ = well_state;
        }





        template <class NonlinearSolverType>
        IterationReport nonlinearIteration(const int /* iteration */,
                                           const double /* dt */,
                                           NonlinearSolverType& /* nonlinear_solver */,
                                           ReservoirState& reservoir_state,
                                           const WellState& well_state)
        {
            // Extract reservoir and well fluxes.
            {
                DebugTimeReport tr("Extracting fluxes");
                extractFluxes(reservoir_state, well_state);
            }

            // Compute cell ordering based on total flux.
            {
                DebugTimeReport tr("Topological sort");
                computeOrdering();
            }

            // Solve in every component (cell or block of cells), in order.
            {
                DebugTimeReport tr("Solving all components");
                solveComponents();
            }

            // Update states for output.

            // Create report and exit.
            const bool failed = false;
            const bool converged = true;
            const int linear_iterations = 0;
            const int well_iterations = std::numeric_limits<int>::min();
            return IterationReport{failed, converged, linear_iterations, well_iterations};
        }





        void afterStep(const double /* dt */,
                       const ReservoirState& /* reservoir_state */,
                       const WellState& /* well_state */)
        {
            // Does nothing in this model.
        }





        using Base::numPhases;


    protected:

        // ============  Data members  ============
        using Base::grid_;
        using Base::ops_;

        ReservoirState reservoir_state0_;
        WellState well_state0_;
        V total_flux_;
        V total_wellperf_flux_;
        DataBlock comp_wellperf_flux_;
        std::vector<int> sequence_;
        std::vector<int> components_;


        // ============  Member functions  ============


        void extractFluxes(const ReservoirState& reservoir_state,
                           const WellState& well_state)
        {
            // Input face fluxes are for interior faces only, while rest of code deals with all faces.
            const V face_flux = Eigen::Map<const V>(reservoir_state.faceflux().data(),
                                                    reservoir_state.faceflux().size());
            using namespace Opm::AutoDiffGrid;
            const int num_faces = numFaces(grid_);
            assert(face_flux.size() < num_faces); // Expected to be internal only.
            total_flux_ = superset(face_flux, ops_.internal_faces, num_faces);
            total_wellperf_flux_ = Eigen::Map<const V>(well_state.perfRates().data(),
                                                       well_state.perfRates().size());
            comp_wellperf_flux_ = Eigen::Map<const DataBlock>(well_state.perfPhaseRates().data(),
                                                              well_state.perfRates().size(),
                                                              numPhases());
            assert(numPhases() * well_state.perfRates().size() == well_state.perfPhaseRates().size());
        }





        void computeOrdering()
        {
            static_assert(std::is_same<Grid, UnstructuredGrid>::value,
                          "compute_sequence() is written in C and therefore requires an UnstructuredGrid, "
                          "it must be rewritten to use other grid classes such as CpGrid");
            using namespace Opm::AutoDiffGrid;
            const int num_cells = numCells(grid_);
            sequence_.resize(num_cells);
            components_.resize(num_cells + 1); // max possible size
            int num_components = -1;
            compute_sequence(&grid_, total_flux_.data(), sequence_.data(), components_.data(), &num_components);
            OpmLog::debug(std::string("Number of components: ") + std::to_string(num_components));
            components_.resize(num_components + 1); // resize to fit actually used part
        }




        void solveComponents()
        {
            const int num_components = components_.size() - 1;
            for (int comp = 0; comp < num_components; ++comp) {
                const int comp_size = components_[comp + 1] - components_[comp];
                if (comp_size == 1) {
                    solveSingleCell(sequence_[components_[comp]]);
                } else {
                    solveMultiCell(comp_size, &sequence_[components_[comp]]);
                }
            }
        }





        void solveSingleCell(const int cell)
        {
        }





        void solveMultiCell(const int comp_size, const int* cell_array)
        {
        }
    };








    /// Providing types by template specialisation of ModelTraits for BlackoilReorderingTransportModel.
    template <class Grid, class WellModel>
    struct ModelTraits< BlackoilReorderingTransportModel<Grid, WellModel> >
    {
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;
    };

} // namespace Opm




#endif // OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED
