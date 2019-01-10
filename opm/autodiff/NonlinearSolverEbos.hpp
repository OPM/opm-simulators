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

#ifndef OPM_NON_LINEAR_SOLVER_EBOS_HPP
#define OPM_NON_LINEAR_SOLVER_EBOS_HPP

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>

#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/propertysystem.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <memory>

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowNonLinearSolver);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NewtonMaxRelax);
NEW_PROP_TAG(FlowNewtonMaxIterations);
NEW_PROP_TAG(FlowNewtonMinIterations);
NEW_PROP_TAG(NewtonRelaxationType);

SET_SCALAR_PROP(FlowNonLinearSolver, NewtonMaxRelax, 0.5);
SET_INT_PROP(FlowNonLinearSolver, FlowNewtonMaxIterations, 20);
SET_INT_PROP(FlowNonLinearSolver, FlowNewtonMinIterations, 1);
SET_STRING_PROP(FlowNonLinearSolver, NewtonRelaxationType, "dampen");

END_PROPERTIES

namespace Opm {


    /// A nonlinear solver class suitable for general fully-implicit models,
    /// as well as pressure, transport and sequential models.
    template <class TypeTag, class PhysicalModel>
    class NonlinearSolverEbos
    {
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    public:
        // Available relaxation scheme types.
        enum RelaxType {
            Dampen,
            SOR
        };

        // Solver parameters controlling nonlinear process.
        struct SolverParameters
        {
            RelaxType relaxType_;
            double relaxMax_;
            double relaxIncrement_;
            double relaxRelTol_;
            int maxIter_; // max nonlinear iterations
            int minIter_; // min nonlinear iterations

            SolverParameters()
            {
                // set default values
                reset();

                // overload with given parameters
                relaxMax_ = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxRelax);
                maxIter_ = EWOMS_GET_PARAM(TypeTag, int, FlowNewtonMaxIterations);
                minIter_ = EWOMS_GET_PARAM(TypeTag, int, FlowNewtonMinIterations);

                const auto& relaxationTypeString = EWOMS_GET_PARAM(TypeTag, std::string, NewtonRelaxationType);
                if (relaxationTypeString == "dampen") {
                    relaxType_ = Dampen;
                } else if (relaxationTypeString == "sor") {
                    relaxType_ = SOR;
                } else {
                    OPM_THROW(std::runtime_error, "Unknown Relaxtion Type " << relaxationTypeString);
                }
            }

            static void registerParameters()
            {
                EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonMaxRelax, "The maximum relaxation factor of a Newton iteration used by flow");
                EWOMS_REGISTER_PARAM(TypeTag, int, FlowNewtonMaxIterations, "The maximum number of Newton iterations per time step used by flow");
                EWOMS_REGISTER_PARAM(TypeTag, int, FlowNewtonMinIterations, "The minimum number of Newton iterations per time step used by flow");
                EWOMS_REGISTER_PARAM(TypeTag, std::string, NewtonRelaxationType, "The type of relaxation used by flow's Newton method");
            }

            void reset()
            {
                // default values for the solver parameters
                relaxType_ = Dampen;
                relaxMax_ = 0.5;
                relaxIncrement_ = 0.1;
                relaxRelTol_ = 0.2;
                maxIter_ = 10;
                minIter_ = 1;
            }

        };

        // Forwarding types from PhysicalModel.
        typedef typename PhysicalModel::WellState WellState;

        // ---------  Public methods  ---------

        /// Construct solver for a given model.
        ///
        /// The model is a std::unique_ptr because the object to which model points to is
        /// not allowed to be deleted as long as the NonlinearSolver object exists.
        ///
        /// \param[in]      param   parameters controlling nonlinear process
        /// \param[in, out] model   physical simulation model.
        NonlinearSolverEbos(const SolverParameters& param,
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


        SimulatorReport step(const SimulatorTimerInterface& timer)
        {
            SimulatorReport iterReport;
            SimulatorReport report;
            failureReport_ = SimulatorReport();

            // Do model-specific once-per-step calculations.
            model_->prepareStep(timer);
            if (timer.initialStep()) {
                model_->adjoint_serialize();
            }
	    //Hopefully this is when wells is fully initialized
	    //if(param_.use_adjoint_){
	      model_->serialize_well(true);
	      //}


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
                    iterReport = model_->nonlinearIteration(iteration, timer, *this);

                    report += iterReport;
                    report.converged = iterReport.converged;

                    converged = report.converged;
                    iteration += 1;
                }
                catch (...) {
                    // if an iteration fails during a time step, all previous iterations
                    // count as a failure as well
                    failureReport_ += report;
                    failureReport_ += model_->failureReport();
                    throw;
                }
            }
            while ( (!converged && (iteration <= maxIter())) || (iteration <= minIter()));

            if (!converged) {
                failureReport_ += report;

                std::string msg = "Solver convergence failure - Failed to complete a time step within " + std::to_string(maxIter()) + " iterations.";
                OPM_THROW_NOLOG(Opm::TooManyIterations, msg);
            }

            // Do model-specific post-step actions.
            model_->afterStep(timer);
            //model_->adjoint_serialize();
            report.converged = true;

            return report;
        }

        /// return the statistics if the step() method failed
        const SimulatorReport& failureReport() const
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

        std::vector<std::vector<double> >
        computeFluidInPlace(const std::vector<int>& fipnum) const
        { return model_->computeFluidInPlace(fipnum); }

        /// Reference to physical model.
        const PhysicalModel& model() const
        { return *model_; }

        /// Mutable reference to physical model.
        PhysicalModel& model()
        { return *model_; }

        /// Detect oscillation or stagnation in a given residual history.
        void detectOscillations(const std::vector<std::vector<double>>& residualHistory,
                                const int it, bool& oscillate, bool& stagnate) const
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
            const std::vector<double>& F0 = residualHistory[it];
            const std::vector<double>& F1 = residualHistory[it - 1];
            const std::vector<double>& F2 = residualHistory[it - 2];
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


        /// Apply a stabilization to dx, depending on dxOld and relaxation parameters.
        /// Implemention for Dune block vectors.
        template <class BVector>
        void stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld, const double omega) const
        {
            // The dxOld is updated with dx.
            // If omega is equal to 1., no relaxtion will be appiled.

            BVector tempDxOld = dxOld;
            dxOld = dx;

            switch (relaxType()) {
            case Dampen: {
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
                OPM_THROW(std::runtime_error, "Can only handle Dampen and SOR relaxation type.");
            }

            return;
        }

        /// The greatest relaxation factor (i.e. smallest factor) allowed.
        double relaxMax() const
        { return param_.relaxMax_; }

        /// The step-change size for the relaxation factor.
        double relaxIncrement() const
        { return param_.relaxIncrement_; }

        /// The relaxation type (Dampen or SOR).
        enum RelaxType relaxType() const
        { return param_.relaxType_; }

        /// The relaxation relative tolerance.
        double relaxRelTol() const
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
        SimulatorReport failureReport_;
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

#endif // OPM_NON_LINEAR_SOLVER_EBOS_HPP
