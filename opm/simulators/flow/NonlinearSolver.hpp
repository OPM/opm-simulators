/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_NON_LINEAR_SOLVER_HPP
#define OPM_NON_LINEAR_SOLVER_HPP

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/nonlinear/newtonmethodproperties.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/basicproperties.hh>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>

#include <memory>

namespace Opm::Parameters {

template<class Scalar>
struct NewtonMaxRelax { static constexpr Scalar value = 0.5; };

struct NewtonMinIterations { static constexpr int value = 2; };
struct NewtonRelaxationType { static constexpr auto value = "dampen"; };

} // namespace Opm::Parameters

namespace Opm {

// Available relaxation scheme types.
enum class NonlinearRelaxType {
    Dampen,
    SOR
};

namespace detail {

/// Detect oscillation or stagnation in a given residual history.
template<class Scalar>
void detectOscillations(const std::vector<std::vector<Scalar>>& residualHistory,
                        const int it, const int numPhases, const Scalar relaxRelTol,
                        const int minimumOscillatingPhases,
                        bool& oscillate, bool& stagnate);

/// Apply a stabilization to dx, depending on dxOld and relaxation parameters.
/// Implemention for Dune block vectors.
template <class BVector, class Scalar>
void stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld,
                              const Scalar omega, NonlinearRelaxType relaxType);

}

// Solver parameters controlling nonlinear process.
template<class Scalar>
struct NonlinearSolverParameters
{
    NonlinearRelaxType relaxType_;
    Scalar relaxMax_;
    Scalar relaxIncrement_;
    Scalar relaxRelTol_;
    int maxIter_; // max nonlinear iterations
    int minIter_; // min nonlinear iterations

    NonlinearSolverParameters();

    static void registerParameters();

    void reset();
};

    /// A nonlinear solver class suitable for general fully-implicit models,
    /// as well as pressure, transport and sequential models.
    template <class TypeTag, class PhysicalModel>
    class NonlinearSolver
    {
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    public:
        using SolverParameters = NonlinearSolverParameters<Scalar>;

        // ---------  Public methods  ---------

        /// Construct solver for a given model.
        ///
        /// The model is a std::unique_ptr because the object to which model points to is
        /// not allowed to be deleted as long as the NonlinearSolver object exists.
        ///
        /// \param[in]      param   parameters controlling nonlinear process
        /// \param[in, out] model   physical simulation model.
        NonlinearSolver(const SolverParameters& param,
                        std::unique_ptr<PhysicalModel> model)
            : param_(param)
            , model_(std::move(model))
            , linearizations_(0)
            , nonlinearIterations_(0)
            , linearIterations_(0)
            , wellIterations_(0)
            , nonlinearIterationsLast_(0)
            , linearIterationsLast_(0)
            , wellIterationsLast_(0)
        {
            if (!model_) {
                OPM_THROW(std::logic_error, "Must provide a non-null model argument for NonlinearSolver.");
            }
        }


        SimulatorReportSingle step(const SimulatorTimerInterface& timer, const TimeStepControlInterface *timeStepControl)
        {
            SimulatorReportSingle report;
            report.global_time = timer.simulationTimeElapsed();
            report.timestep_length = timer.currentStepLength();

            // Do model-specific once-per-step calculations.
            report += model_->prepareStep(timer);

            int iteration = 0;

            // Let the model do one nonlinear iteration.

            // Set up for main solver loop.
            bool converged = false;

            // ----------  Main nonlinear solver loop  ----------
            do {
                try {
                    // Do the nonlinear step. If we are in a converged state, the
                    // model will usually do an early return without an expensive
                    // solve, unless the minIter() count has not been reached yet.
                    auto iterReport = model_->nonlinearIteration(iteration, timer, *this);
                    iterReport.global_time = timer.simulationTimeElapsed();
                    report += iterReport;
                    report.converged = iterReport.converged;

                    converged = report.converged;
                    iteration += 1;
                }
                catch (...) {
                    // if an iteration fails during a time step, all previous iterations
                    // count as a failure as well
                    failureReport_ = report;
                    failureReport_ += model_->failureReport();
                    throw;
                }
            }
            while ( (!converged && (iteration <= maxIter())) || (iteration <= minIter()));

            if (!converged) {
                failureReport_ = report;

                std::string msg = "Solver convergence failure - Failed to complete a time step within " + std::to_string(maxIter()) + " iterations.";
                OPM_THROW_NOLOG(TooManyIterations, msg);
            }
            auto relativeChange = model_->relativeChange();
            if (timeStepControl && !timeStepControl->timeStepAccepted(relativeChange, timer.currentStepLength())) {
                report.converged = false;
                report.time_step_rejected = true;
                failureReport_ = report;

                std::string msg = "Relative change in solution for time step was " + std::to_string(relativeChange) + ", which is larger than the tolerance accepted by the timestepping algorithm.";
                OPM_THROW_NOLOG(TimeSteppingBreakdown, msg);
            }

            // Do model-specific post-step actions.
            report += model_->afterStep(timer);
            report.converged = true;
            return report;
        }

        /// return the statistics if the step() method failed
        const SimulatorReportSingle& failureReport() const
        { return failureReport_; }

        /// Number of linearizations used in all calls to step().
        int linearizations() const
        { return linearizations_; }

        /// Number of full nonlinear solver iterations used in all calls to step().
        int nonlinearIterations() const
        { return nonlinearIterations_; }

        /// Number of linear solver iterations used in all calls to step().
        int linearIterations() const
        { return linearIterations_; }

        /// Number of well iterations used in all calls to step().
        int wellIterations() const
        { return wellIterations_; }

        /// Number of nonlinear solver iterations used in the last call to step().
        int nonlinearIterationsLastStep() const
        { return nonlinearIterationsLast_; }

        /// Number of linear solver iterations used in the last call to step().
        int linearIterationsLastStep() const
        { return linearIterationsLast_; }

        /// Number of well iterations used in all calls to step().
        int wellIterationsLastStep() const
        { return wellIterationsLast_; }

        std::vector<std::vector<Scalar> >
        computeFluidInPlace(const std::vector<int>& fipnum) const
        { return model_->computeFluidInPlace(fipnum); }

        /// Reference to physical model.
        const PhysicalModel& model() const
        { return *model_; }

        /// Mutable reference to physical model.
        PhysicalModel& model()
        { return *model_; }

        /// Detect oscillation or stagnation in a given residual history.
        void detectOscillations(const std::vector<std::vector<Scalar>>& residualHistory,
                                const int it, bool& oscillate, bool& stagnate) const
        {
            detail::detectOscillations(residualHistory, it, model_->numPhases(),
                                       this->relaxRelTol(), 2, oscillate, stagnate);
        }


        /// Apply a stabilization to dx, depending on dxOld and relaxation parameters.
        /// Implemention for Dune block vectors.
        template <class BVector>
        void stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld, const Scalar omega) const
        {
            detail::stabilizeNonlinearUpdate(dx, dxOld, omega, this->relaxType());
        }

        /// The greatest relaxation factor (i.e. smallest factor) allowed.
        Scalar relaxMax() const
        { return param_.relaxMax_; }

        /// The step-change size for the relaxation factor.
        Scalar relaxIncrement() const
        { return param_.relaxIncrement_; }

        /// The relaxation type (Dampen or SOR).
        NonlinearRelaxType relaxType() const
        { return param_.relaxType_; }

        /// The relaxation relative tolerance.
        Scalar relaxRelTol() const
        { return param_.relaxRelTol_; }

        /// The maximum number of nonlinear iterations allowed.
        int maxIter() const
        { return param_.maxIter_; }

        /// The minimum number of nonlinear iterations allowed.
        int minIter() const
        { return param_.minIter_; }

        /// Set parameters to override those given at construction time.
        void setParameters(const SolverParameters& param)
        { param_ = param; }

    private:
        // ---------  Data members  ---------
        SimulatorReportSingle failureReport_;
        SolverParameters param_;
        std::unique_ptr<PhysicalModel> model_;
        int linearizations_;
        int nonlinearIterations_;
        int linearIterations_;
        int wellIterations_;
        int nonlinearIterationsLast_;
        int linearIterationsLast_;
        int wellIterationsLast_;
    };

} // namespace Opm

#endif // OPM_NON_LINEAR_SOLVER_HPP
