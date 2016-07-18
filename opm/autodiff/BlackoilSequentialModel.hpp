/*
  Copyright 2015, 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil AS.

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

#ifndef OPM_BLACKOILSEQUENTIALMODEL_HEADER_INCLUDED
#define OPM_BLACKOILSEQUENTIALMODEL_HEADER_INCLUDED


#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/autodiff/BlackoilPressureModel.hpp>
#include <opm/autodiff/BlackoilTransportModel.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>

namespace Opm {

    struct BlackoilSequentialModelParameters : public BlackoilModelParameters
    {
        bool iterate_to_fully_implicit;
        explicit BlackoilSequentialModelParameters( const parameter::ParameterGroup& param )
            : BlackoilModelParameters(param),
              iterate_to_fully_implicit(param.getDefault("iterate_to_fully_implicit", false))
        {
        }
    };


    /// A sequential splitting model implementation for three-phase black oil.
    template<class Grid, class WellModel>
    class BlackoilSequentialModel
    {
    public:
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilSequentialModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilSequentialModel(const ModelParameters&          param,
                                const Grid&                     grid ,
                                const BlackoilPropsAdInterface& fluid,
                                const DerivedGeology&           geo  ,
                                const RockCompressibility*      rock_comp_props,
                                const WellModel                 well_model,
                                const NewtonIterationBlackoilInterface& linsolver,
                                Opm::EclipseStateConstPtr eclState,
                                const bool has_disgas,
                                const bool has_vapoil,
                                const bool terminal_output)
        : pressure_model_(new PressureModel(param, grid, fluid, geo, rock_comp_props, well_model,
                                            linsolver, eclState, has_disgas, has_vapoil, terminal_output)),
          transport_model_(new TransportModel(param, grid, fluid, geo, rock_comp_props, well_model,
                                              linsolver, eclState, has_disgas, has_vapoil, terminal_output)),
          pressure_solver_(typename PressureSolver::SolverParameters(), std::move(pressure_model_)),
          transport_solver_(typename TransportSolver::SolverParameters(), std::move(transport_model_)),
          initial_reservoir_state_(0, 0, 0), // will be overwritten
          iterate_to_fully_implicit_(param.iterate_to_fully_implicit)
        {
        }





        /// Called once before each time step.
        /// \param[in] timer             simulation timer
        /// \param[in] reservoir_state   reservoir state variables
        /// \param[in] well_state        well state variables
        void prepareStep(const SimulatorTimerInterface& /*timer*/,
                         const ReservoirState& reservoir_state,
                         const WellState& well_state)
        {
            initial_reservoir_state_ = reservoir_state;
            initial_well_state_ = well_state;
        }





        /// Called once per nonlinear iteration.
        /// This model will first solve the pressure model to convergence, then the
        /// transport model.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer             simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        template <class NonlinearSolverType>
        IterationReport nonlinearIteration(const int iteration,
                                           const SimulatorTimerInterface& timer,
                                           NonlinearSolverType& /* nonlinear_solver */,
                                           ReservoirState& reservoir_state,
                                           WellState& well_state)
        {
            if (!iterate_to_fully_implicit_) {
                // Do a single pressure solve, followed by a single transport solve.
                if (terminalOutputEnabled()) {
                    OpmLog::info("Using sequential model.");
                }

                // Pressure solve.
                if (terminalOutputEnabled()) {
                    OpmLog::info("Solving the pressure equation.");
                }
                ReservoirState initial_state = reservoir_state;
                const int pressure_liniter = pressure_solver_.step(timer, reservoir_state, well_state);
                if (pressure_liniter == -1) {
                    OPM_THROW(std::runtime_error, "Pressure solver failed to converge.");
                }

                // Transport solve.
                if (terminalOutputEnabled()) {
                    OpmLog::info("Solving the transport equations.");
                }
                const int transport_liniter = transport_solver_.step(timer, initial_state, well_state, reservoir_state, well_state);
                if (transport_liniter == -1) {
                    OPM_THROW(std::runtime_error, "Transport solver failed to converge.");
                }

                // Report and return.
                return IterationReport { false, true, pressure_liniter + transport_liniter, 0 };
            } else {
                // Iterate to fully implicit solution.
                // This call is just for a single iteration (one pressure and one transport solve),
                // we return a 'false' converged status if more are needed
                if (terminalOutputEnabled()) {
                    OpmLog::info("Using sequential model in iterative mode, outer iteration " + std::to_string(iteration));
                }

                // Pressure solve.
                if (terminalOutputEnabled()) {
                    OpmLog::info("Solving the pressure equation.");
                }
                const int pressure_liniter = pressure_solver_.step(timer, initial_reservoir_state_, initial_well_state_, reservoir_state, well_state);
                if (pressure_liniter == -1) {
                    OPM_THROW(std::runtime_error, "Pressure solver failed to converge.");
                }

                // Transport solve.
                if (terminalOutputEnabled()) {
                    OpmLog::info("Solving the transport equations.");
                }
                const int transport_liniter = transport_solver_.step(timer, initial_reservoir_state_, initial_well_state_, reservoir_state, well_state);
                if (transport_liniter == -1) {
                    OPM_THROW(std::runtime_error, "Transport solver failed to converge.");
                }

                // Report and return.
                const bool converged = iteration >= 3; // TODO: replace this with a proper convergence check
                return IterationReport { false, converged, pressure_liniter + transport_liniter, 0 };
            }
        }





        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void afterStep(const SimulatorTimerInterface& /* timer */,
                       ReservoirState& /* reservoir_state */,
                       WellState& /* well_state */)
        {
        }





        /// \brief Set threshold pressures that prevent or reduce flow.
        /// This prevents flow across faces if the potential
        /// difference is less than the threshold. If the potential
        /// difference is greater, the threshold value is subtracted
        /// before calculating flow. This is treated symmetrically, so
        /// flow is prevented or reduced in both directions equally.
        /// \param[in]  threshold_pressures_by_face   array of size equal to the number of faces
        ///                                   of the grid passed in the constructor.
        void setThresholdPressures(const std::vector<double>& threshold_pressures_by_face)
        {
            pressure_solver_.model().setThresholdPressures(threshold_pressures_by_face);
            transport_solver_.model().setThresholdPressures(threshold_pressures_by_face);
        }





        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return pressure_solver_.model().terminalOutputEnabled();
        }





        /// Return the relative change in variables relevant to this model.
        double relativeChange(const SimulationDataContainer& previous,
                              const SimulationDataContainer& current ) const
        {
            // TODO: this is a quick stopgap implementation, and should be evaluated more carefully.
            return std::max(pressure_solver_.model().relativeChange(previous, current),
                            transport_solver_.model().relativeChange(previous, current));
        }

        /// Return the well model
        const WellModel& wellModel() const
        {
            return pressure_model_->wellModel();
        }


        /// Compute fluid in place.
        V computeFluidInPlace(const ReservoirState& x,
                              const WellState& xw) const
        {
            return transport_solver_.computeFluidInPlace(x, xw);
        }



    protected:
        typedef BlackoilPressureModel<Grid, WellModel> PressureModel;
        typedef BlackoilTransportModel<Grid, WellModel> TransportModel;
        typedef NonlinearSolver<PressureModel> PressureSolver;
        typedef NonlinearSolver<TransportModel> TransportSolver;

        std::unique_ptr<PressureModel> pressure_model_;
        std::unique_ptr<TransportModel> transport_model_;
        PressureSolver pressure_solver_;
        TransportSolver transport_solver_;

        ReservoirState initial_reservoir_state_;
        WellState initial_well_state_;

        bool iterate_to_fully_implicit_;
    };

} // namespace Opm



#endif // OPM_BLACKOILSEQUENTIALMODEL_HEADER_INCLUDED
