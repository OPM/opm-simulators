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
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <memory>

namespace Opm::Properties {

NEW_TYPE_TAG(FlowNonLinearSolver);
NEW_PROP_TAG(NewtonMaxRelax);
NEW_PROP_TAG(FlowNewtonMaxIterations);
NEW_PROP_TAG(FlowNewtonMinIterations);
NEW_PROP_TAG(NewtonRelaxationType);

SET_SCALAR_PROP(FlowNonLinearSolver, NewtonMaxRelax, 0.5);
SET_INT_PROP(FlowNonLinearSolver, FlowNewtonMaxIterations, 20);
SET_INT_PROP(FlowNonLinearSolver, FlowNewtonMinIterations, 1);
SET_STRING_PROP(FlowNonLinearSolver, NewtonRelaxationType, "dampen");

} // namespace Opm::Properties

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


        SimulatorReportSingle stepSequential(const SimulatorTimerInterface& timer,bool implicit){
            LinearizationType linearizationType;
            bool converged = false;
            SimulatorReportSingle reportpre;
            SimulatorReportSingle reportseq;
            int seqiterations  = 0;
            int maxseqiterations = 10;
            // should probalby store for trying onse more
            auto solutionOld = model_->ebosSimulator().model().solution(/*timeIdx=*/0);
            auto oldTotalSaturation = model_->ebosSimulator().problem().getTotalSaturation();
            bool first= true;
            while(not(converged) && (seqiterations < maxseqiterations) ){
                // do pressure step
                if(seqiterations>0){
                    auto& prevsol = model_->ebosSimulator().model().solution(/*timeIdx=*/1);
                    model_->ebosSimulator().model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/1);
                    //NB is the storag cach allso invalidated??
                    prevsol=solutionOld;
                }
                linearizationType.type = Opm::LinearizationType::pressure;                
                model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);  
                model_->updateSolution();                               
                SimulatorReportSingle lreportpre = this->stepDefault(timer, /*next*/first, false);
                first=false;
                // do transport step
                model_->ebosSimulator().problem().wellModel().endTimeStep();// shouldupdate previous state NB do wee need to store it if transportstep fails?
                model_->ebosSimulator().problem().updatePressureAndFluxes();//update pressure ans ST 
                linearizationType.type = Opm::LinearizationType::seqtransport;
                model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);

                model_->updateSolution();
                reportpre += lreportpre;
                try{
                    //model_->ebosSimulator().model().advanceTimeLevel();
                    //model_->ebosSimulator().problem().beginTimeStep();
                    SimulatorReportSingle lreportseq = this->stepDefault(timer,/*next*/false, false);
                    // hopefully do not change any thing only update well quantities
                    //model_->wellModel().solveWells(model_->ebosSimulator().timeStepSize());
                    reportseq += lreportseq;
                }catch (...){
                    //set the prevois value to the staring point
                    auto& prevsol = model_->ebosSimulator().model().solution(/*timeIdx=*/1);
                    prevsol = solutionOld;
                    auto& currsol =  model_->ebosSimulator().model().solution(/*timeIdx=*/0);
                    currsol = solutionOld;
                    // set back totalSaturation pressure and fluxes
                    model_->ebosSimulator().problem().setTotalSaturation(oldTotalSaturation);
                    // total pressure  fluxes should be updated in the pressure solve anyway
                    
                    throw;
                }
                
                // for no not seq implicit
                if(implicit){
                    // for now do full implicit solve
                    //NB storagecache need to be invalidated
                    //bool storagecache = EWOMS_GET_PARAM(TypeTag, bool, EnableStorageCache);
                    model_->ebosSimulator().model().setEnableStorageCache(false);
                    auto& prevsol = model_->ebosSimulator().model().solution(/*timeIdx=*/1);
                    prevsol=solutionOld;
                    model_->ebosSimulator().model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/1);
                    model_->ebosSimulator().model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
                    linearizationType.type = Opm::LinearizationType::implicit;
                    model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);                   
                    model_->updateSolution();
                    
                    //NB her only the convergenceshould be checked
                    //SimulatorReportSingle lreportsim = this->stepDefault(timer,/*next*/false, false);
                    std::vector<double> residual_norms;
                    auto asreport = model_->assembleReservoir(timer, seqiterations);
                    auto convrep = model_->getConvergence(timer, seqiterations,residual_norms);
                    converged = convrep.converged();
                    // if this used one probably have to invalidate the storage cache
                    //model_->ebosSimulator().model().setEnableStorageCache(storagecache);
                    model_->ebosSimulator().problem().totalSaturationOne();
                }else{
                    converged = true;
                }
                
                std::cout << "Sequantial fullimplicit iteration " << seqiterations << std::endl;
                seqiterations += 1;
            }
            SimulatorReportSingle report;
            if(not(converged)){
                report.converged = false;
            }else{
                // model_->updateSolution();
                model_->ebosSimulator().problem().updateTotalSaturation();  
                model_->afterStep(timer);
                // NOTE: updating solution again for information used in the relativeChange function.
                linearizationType.type = Opm::LinearizationType::pressure;
                model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);
                model_->updateSolution();
                report.converged = true;
            }            
            //todo fill this report correctly
            return report;
                
        }
        SimulatorReportSingle step(const SimulatorTimerInterface& timer){
            std::string simulationtype  = EWOMS_GET_PARAM(TypeTag, std::string, SimulationType);
            SimulatorReportSingle report;
            if(simulationtype == "implicit"){
                report = this->stepFull(timer);
            }else if(simulationtype == "pressure"){
                report = this->stepPressure(timer);
            }else if(simulationtype == "seq"){
                report = this->stepSequential(timer,false);
            }else if(simulationtype == "fimseq"){
                report = this->stepSequential(timer,true);
            }else{
                OPM_THROW(std::logic_error, "No such simulation type valid is: implicit, pressure, seq,  fimseq ");
            }
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
        SimulatorReportSingle stepFull(const SimulatorTimerInterface& timer){
            LinearizationType linearizationType;// use default
            model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);
            // incase previous step was not fully implicit
            model_->ebosSimulator().problem().updatePressureAndFluxes();//update pressure ans ST                               
            model_->updateSolution();
            SimulatorReportSingle report = this->stepDefault(timer,true,true);
            return report;
        }

        SimulatorReportSingle stepPressure(const SimulatorTimerInterface& timer){
            LinearizationType linearizationType;
            SimulatorReportSingle reportpre;
            // do pressure step
            linearizationType.type = Opm::LinearizationType::pressure;
            model_->ebosSimulator().problem().updatePressureAndFluxes();//update pressure ans ST
            model_->ebosSimulator().model().linearizer().setLinearizationType(linearizationType);  
            model_->updateSolution();// should conver to the new solution time         
            SimulatorReportSingle lreportpre = this->stepDefault(timer,true,true);
            SimulatorReportSingle report = lreportpre;
            //todo fill this report correctly
            return report;
                
        }
                
        //SimulatorReportSingle stepDefault(const SimulatorTimerInterface& timer,bool next = true, bool end)
        SimulatorReportSingle stepDefault(const SimulatorTimerInterface& timer,bool next, bool endstep)    
        {
            SimulatorReportSingle report;
            report.global_time = timer.simulationTimeElapsed();
            report.timestep_length = timer.currentStepLength();

            // Do model-specific once-per-step calculations.
            
            model_->prepareStep(timer, next);
           

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
                OPM_THROW_NOLOG(Opm::TooManyIterations, msg);
            }

            // Do model-specific post-step actions.
            if(endstep){
                model_->afterStep(timer);
            }
            report.converged = true;
            return report;
        }

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

#endif // OPM_NON_LINEAR_SOLVER_EBOS_HPP
