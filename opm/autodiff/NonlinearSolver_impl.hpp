/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
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

#ifndef OPM_NONLINEARSOLVER_IMPL_HEADER_INCLUDED
#define OPM_NONLINEARSOLVER_IMPL_HEADER_INCLUDED

#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

namespace Opm
{
    template <class PhysicalModel>
    NonlinearSolver<PhysicalModel>::NonlinearSolver(const SolverParameters& param,
                                                    std::unique_ptr<PhysicalModel> model_arg)
        : param_(param),
          model_(std::move(model_arg)),
          linearizations_(0),
          nonlinearIterations_(0),
          linearIterations_(0),
          wellIterations_(0),
          nonlinearIterationsLast_(0),
          linearIterationsLast_(0),
          wellIterationsLast_(0)
    {
        if (!model_) {
            OPM_THROW(std::logic_error, "Must provide a non-null model argument for NonlinearSolver.");
        }
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::linearizations() const
    {
        return linearizations_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::nonlinearIterations() const
    {
        return nonlinearIterations_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::linearIterations() const
    {
        return linearIterations_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::wellIterations() const
    {
        return wellIterations_;
    }

    template <class PhysicalModel>
    const PhysicalModel& NonlinearSolver<PhysicalModel>::model() const
    {
        return *model_;
    }

    template <class PhysicalModel>
    PhysicalModel& NonlinearSolver<PhysicalModel>::model()
    {
        return *model_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::nonlinearIterationsLastStep() const
    {
        return nonlinearIterationsLast_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::linearIterationsLastStep() const
    {
        return linearIterationsLast_;
    }

    template <class PhysicalModel>
    int NonlinearSolver<PhysicalModel>::wellIterationsLastStep() const
    {
        return wellIterationsLast_;
    }

    template <class PhysicalModel>
    SimulatorReport
    NonlinearSolver<PhysicalModel>::
    step(const SimulatorTimerInterface& timer,
         ReservoirState& reservoir_state,
         WellState& well_state)
    {
        return step(timer, reservoir_state, well_state, reservoir_state, well_state);
    }



    template <class PhysicalModel>
    SimulatorReport
    NonlinearSolver<PhysicalModel>::
    step(const SimulatorTimerInterface& timer,
         const ReservoirState& initial_reservoir_state,
         const WellState& initial_well_state,
         ReservoirState& reservoir_state,
         WellState& well_state)
    {
        SimulatorReport iterReport;
        SimulatorReport report;

        // Do model-specific once-per-step calculations.
        model_->prepareStep(timer, initial_reservoir_state, initial_well_state);

        int iteration = 0;

        // Let the model do one nonlinear iteration.

        // Set up for main solver loop.
        bool converged = false;

        // ----------  Main nonlinear solver loop  ----------
        do {
            // Do the nonlinear step. If we are in a converged state, the
            // model will usually do an early return without an expensive
            // solve, unless the minIter() count has not been reached yet.
            iterReport = model_->nonlinearIteration(iteration, timer, *this, reservoir_state, well_state);

            report += iterReport;
            report.converged = iterReport.converged;

            converged = report.converged;
            iteration += 1;
        } while ( (!converged && (iteration <= maxIter())) || (iteration <= minIter()));

        if (!converged) {
            std::string msg = "Failed to complete a time step within " + std::to_string(maxIter()) + " iterations.";
            if (model_->terminalOutputEnabled()) {
                OpmLog::problem(msg);
            }
            OPM_THROW_NOLOG(Opm::NumericalProblem, msg);
        }

        // Do model-specific post-step actions.
        model_->afterStep(timer, reservoir_state, well_state);
        report.converged = true;

        return report;
    }



    template <class PhysicalModel>
    void NonlinearSolver<PhysicalModel>::SolverParameters::
    reset()
    {
        // default values for the solver parameters
        relax_type_      = DAMPEN;
        relax_max_       = 0.5;
        relax_increment_ = 0.1;
        relax_rel_tol_   = 0.2;
        max_iter_        = 15;
        min_iter_        = 1;
    }

    template <class PhysicalModel>
    NonlinearSolver<PhysicalModel>::SolverParameters::
    SolverParameters()
    {
        // set default values
        reset();
    }

    template <class PhysicalModel>
    NonlinearSolver<PhysicalModel>::SolverParameters::
    SolverParameters( const parameter::ParameterGroup& param )
    {
        // set default values
        reset();

        // overload with given parameters
        relax_max_   = param.getDefault("relax_max", relax_max_);
        max_iter_    = param.getDefault("max_iter", max_iter_);
        min_iter_    = param.getDefault("min_iter", min_iter_);

        std::string relaxation_type = param.getDefault("relax_type", std::string("dampen"));
        if (relaxation_type == "dampen") {
            relax_type_ = DAMPEN;
        } else if (relaxation_type == "sor") {
            relax_type_ = SOR;
        } else {
            OPM_THROW(std::runtime_error, "Unknown Relaxtion Type " << relaxation_type);
        }
    }

    template <class PhysicalModel>
    void
    NonlinearSolver<PhysicalModel>::detectOscillations(const std::vector<std::vector<double>>& residual_history,
                                                       const int it,
                                                       bool& oscillate, bool& stagnate) const
    {
        // The detection of oscillation in two primary variable results in the report of the detection
        // of oscillation for the solver.
        // Only the saturations are used for oscillation detection for the black oil model.
        // Stagnate is not used for any treatment here.

        if ( it < 2 ) {
            oscillate = false;
            stagnate = false;
            return;
        }

        stagnate = true;
        int oscillatePhase = 0;
        const std::vector<double>& F0 = residual_history[it];
        const std::vector<double>& F1 = residual_history[it - 1];
        const std::vector<double>& F2 = residual_history[it - 2];
        for (int p= 0; p < model_->numPhases(); ++p){
            const double d1 = std::abs((F0[p] - F2[p]) / F0[p]);
            const double d2 = std::abs((F0[p] - F1[p]) / F0[p]);

            oscillatePhase += (d1 < relaxRelTol()) && (relaxRelTol() < d2);

            // Process is 'stagnate' unless at least one phase
            // exhibits significant residual change.
            stagnate = (stagnate && !(std::abs((F1[p] - F2[p]) / F2[p]) > 1.0e-3));
        }

        oscillate = (oscillatePhase > 1);
    }


    template <class PhysicalModel>
    template <class BVector>
    void
    NonlinearSolver<PhysicalModel>::stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld, const double omega) const
    {
        // The dxOld is updated with dx.
        // If omega is equal to 1., no relaxtion will be appiled.

        BVector tempDxOld = dxOld;
        dxOld = dx;

        switch (relaxType()) {
            case DAMPEN: {
                if (omega == 1.) {
                    return;
                }
                auto i = dx.size();
                for (i = 0; i < dx.size(); ++i) {
                    dx[i] *= omega;
                }
                return;
            }
            case SOR: {
                if (omega == 1.) {
                    return;
                }
                auto i = dx.size();
                for (i = 0; i < dx.size(); ++i) {
                    dx[i] *= omega;
                    tempDxOld[i] *= (1.-omega);
                    dx[i] += tempDxOld[i];
                }
                return;
            }
            default:
                OPM_THROW(std::runtime_error, "Can only handle DAMPEN and SOR relaxation type.");
        }

        return;
    }
} // namespace Opm


#endif // OPM_FULLYIMPLICITSOLVER_IMPL_HEADER_INCLUDED
