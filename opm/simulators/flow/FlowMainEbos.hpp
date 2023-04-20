/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
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
#ifndef OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
#define OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED

#include <opm/simulators/flow/Banners.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <memory>
#include <string_view>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableDryRun {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct OutputInterval {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableLoggingFalloutWarning {
    using type = UndefinedProperty;
};

// TODO: enumeration parameters. we use strings for now.
template<class TypeTag>
struct EnableDryRun<TypeTag, TTag::EclFlowProblem> {
    static constexpr auto value = "auto";
};
// Do not merge parallel output files or warn about them
template<class TypeTag>
struct EnableLoggingFalloutWarning<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct OutputInterval<TypeTag, TTag::EclFlowProblem> {
    static constexpr int value = 1;
};

} // namespace Opm::Properties

namespace Opm {
namespace detail {

void mergeParallelLogFiles(std::string_view output_dir,
                           std::string_view deck_filename,
                           bool enableLoggingFalloutWarning);

void handleExtraConvergenceOutput(SimulatorReport& report,
                                  std::string_view option,
                                  std::string_view optionName,
                                  std::string_view output_dir,
                                  std::string_view base_name);

}

    class Deck;

    // The FlowMain class is the ebos based black-oil simulator.
    template <class TypeTag>
    class FlowMainEbos
    {
    public:
        using MaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
        using EbosSimulator = GetPropType<TypeTag, Properties::Simulator>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        using Simulator = SimulatorFullyImplicitBlackoilEbos<TypeTag>;

        FlowMainEbos(int argc, char **argv, bool output_cout, bool output_files )
            : argc_{argc}, argv_{argv},
              output_cout_{output_cout}, output_files_{output_files}
        {

        }

        // Read the command line parameters. Throws an exception if something goes wrong.
        static int setupParameters_(int argc, char** argv, Parallel::Communication comm)
        {
            using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
            if (!ParamsMeta::registrationOpen()) {
                // We have already successfully run setupParameters_().
                // For the dynamically chosen runs (as from the main flow
                // executable) we must run this function again with the
                // real typetag to be used, as the first time was with the
                // "FlowEarlyBird" typetag. However, for the static ones (such
                // as 'flow_onephase_energy') it has already been run with the
                // correct typetag.
                return EXIT_SUCCESS;
            }
            // register the flow specific parameters
            EWOMS_REGISTER_PARAM(TypeTag, std::string, EnableDryRun,
                                 "Specify if the simulation ought to be actually run, or just pretended to be");
            EWOMS_REGISTER_PARAM(TypeTag, int, OutputInterval,
                                 "Specify the number of report steps between two consecutive writes of restart data");
            EWOMS_REGISTER_PARAM(TypeTag, bool, EnableLoggingFalloutWarning,
                                 "Developer option to see whether logging was on non-root processors. In that case it will be appended to the *.DBG or *.PRT files");

            Simulator::registerParameters();

            // register the parameters inherited from ebos
            registerAllParameters_<TypeTag>(/*finalizeRegistration=*/false);

            // hide the parameters unused by flow. TODO: this is a pain to maintain
            EWOMS_HIDE_PARAM(TypeTag, EnableGravity);
            EWOMS_HIDE_PARAM(TypeTag, EnableGridAdaptation);

            // this parameter is actually used in eWoms, but the flow well model
            // hard-codes the assumption that the intensive quantities cache is enabled,
            // so flow crashes. Let's hide the parameter for that reason.
            EWOMS_HIDE_PARAM(TypeTag, EnableIntensiveQuantityCache);

            // thermodynamic hints are not implemented/required by the eWoms blackoil
            // model
            EWOMS_HIDE_PARAM(TypeTag, EnableThermodynamicHints);

            // in flow only the deck file determines the end time of the simulation
            EWOMS_HIDE_PARAM(TypeTag, EndTime);

            // time stepping is not done by the eWoms code in flow
            EWOMS_HIDE_PARAM(TypeTag, InitialTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepDivisions);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MinTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, PredeterminedTimeStepsFile);

            EWOMS_HIDE_PARAM(TypeTag, EclMaxTimeStepSizeAfterWellEvent);
            EWOMS_HIDE_PARAM(TypeTag, EclRestartShrinkFactor);
            EWOMS_HIDE_PARAM(TypeTag, EclEnableTuning);

            // flow also does not use the eWoms Newton method
            EWOMS_HIDE_PARAM(TypeTag, NewtonMaxError);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTolerance);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTargetIterations);
            EWOMS_HIDE_PARAM(TypeTag, NewtonVerbose);
            EWOMS_HIDE_PARAM(TypeTag, NewtonWriteConvergence);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumTolerance);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumToleranceExponent);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonStrictIterations);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedVolumeFraction);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedTolerance);

            // the default eWoms checkpoint/restart mechanism does not work with flow
            EWOMS_HIDE_PARAM(TypeTag, RestartTime);
            EWOMS_HIDE_PARAM(TypeTag, RestartWritingInterval);
            // hide all vtk related it is not currently possible to do this dependet on if the vtk writing is used
            //if(not(EWOMS_GET_PARAM(TypeTag,bool,EnableVtkOutput))){
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteOilFormationVolumeFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteOilSaturationPressure);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteOilVaporizationFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWritePorosity);
                EWOMS_HIDE_PARAM(TypeTag, VtkWritePotentialGradients);
                EWOMS_HIDE_PARAM(TypeTag, VtkWritePressures);
                EWOMS_HIDE_PARAM(TypeTag, VtkWritePrimaryVars);
                EWOMS_HIDE_PARAM(TypeTag, VtkWritePrimaryVarsMeaning);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteProcessRank);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteRelativePermeabilities);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteSaturatedGasOilVaporizationFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteSaturatedOilGasDissolutionFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteSaturationRatios);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteSaturations);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteTemperature);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteViscosities);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteWaterFormationVolumeFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteGasDissolutionFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteGasFormationVolumeFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteGasSaturationPressure);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteIntrinsicPermeabilities);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteEclTracerConcentration);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteExtrusionFactor);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteFilterVelocities);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteDensities);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteDofIndex);
                EWOMS_HIDE_PARAM(TypeTag, VtkWriteMobilities);
                //}
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteAverageMolarMasses);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteFugacities);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteFugacityCoeffs);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteMassFractions);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteMolarities);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteMoleFractions);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteTotalMassFractions);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteTotalMoleFractions);

            EWOMS_HIDE_PARAM(TypeTag, VtkWriteTortuosities);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteDiffusionCoefficients);
            EWOMS_HIDE_PARAM(TypeTag, VtkWriteEffectiveDiffusionCoefficients);

            EWOMS_END_PARAM_REGISTRATION(TypeTag);

            int mpiRank = comm.rank();

            // read in the command line parameters
            int status = ::Opm::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*doRegistration=*/false, /*allowUnused=*/true, /*handleHelp=*/(mpiRank==0));
            if (status == 0) {

                // deal with unknown parameters.

                int unknownKeyWords = 0;
                if (mpiRank == 0) {
                    unknownKeyWords = Parameters::printUnused<TypeTag>(std::cerr);
                }
                int globalUnknownKeyWords = comm.sum(unknownKeyWords);
                unknownKeyWords = globalUnknownKeyWords;
                if ( unknownKeyWords )
                {
                    if ( mpiRank == 0 )
                    {
                        std::string msg = "Aborting simulation due to unknown "
                            "parameters. Please query \"flow --help\" for "
                            "supported command line parameters.";
                        if (OpmLog::hasBackend("STREAMLOG"))
                        {
                            OpmLog::error(msg);
                        }
                        else {
                            std::cerr << msg << std::endl;
                        }
                    }
                    return EXIT_FAILURE;
                }

                // deal with --print-properties and --print-parameters and unknown parameters.

                bool doExit = false;

                if (EWOMS_GET_PARAM(TypeTag, int, PrintProperties) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        Properties::printValues<TypeTag>(std::cout);
                }

                if (EWOMS_GET_PARAM(TypeTag, int, PrintParameters) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        Parameters::printValues<TypeTag>();
                }

                if (doExit)
                    return -1;
            }

            return status;
        }

        /// This is the main function of Flow.  It runs a complete simulation with the
        /// given grid and simulator classes, based on the user-specified command-line
        /// input.
        int execute()
        {
            return execute_(&FlowMainEbos::runSimulator, /*cleanup=*/true);
        }

        int executeInitStep()
        {
            return execute_(&FlowMainEbos::runSimulatorInit, /*cleanup=*/false);
        }

        // Returns true unless "EXIT" was encountered in the schedule
        //   section of the input datafile.
        int executeStep()
        {
            return simulator_->runStep(*simtimer_);
        }

        // Called from Python to cleanup after having executed the last
        // executeStep()
        int executeStepsCleanup()
        {
            SimulatorReport report = simulator_->finalize();
            runSimulatorAfterSim_(report);
            return report.success.exit_status;
        }

        EbosSimulator *getSimulatorPtr() {
            return ebosSimulator_.get();
        }

        SimulatorTimer* getSimTimer() {
            return simtimer_.get();
        }

    private:
        // called by execute() or executeInitStep()
        int execute_(int (FlowMainEbos::* runOrInitFunc)(), bool cleanup)
        {
            try {
                // deal with some administrative boilerplate

                int status = setupParameters_(this->argc_, this->argv_, EclGenericVanguard::comm());
                if (status)
                    return status;

                setupParallelism();
                setupEbosSimulator();
                createSimulator();

                // if run, do the actual work, else just initialize
                int exitCode = (this->*runOrInitFunc)();
                if (cleanup) {
                    executeCleanup_();
                }
                return exitCode;
            }
            catch (const std::exception& e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if (this->output_cout_) {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }
#if HAVE_MPI
                if (this->mpi_size_ > 1)
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
                return EXIT_FAILURE;
            }
        }

        void executeCleanup_() {
            // clean up
            mergeParallelLogFiles();
        }

    protected:
        void setupParallelism()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized
            // here. (yes, the name of this method is misleading.)
            auto comm = EclGenericVanguard::comm();
            mpi_rank_ = comm.rank();
            mpi_size_ = comm.size();

#if _OPENMP
            // if openMP is available, default to 2 threads per process.
            if (!getenv("OMP_NUM_THREADS"))
                omp_set_num_threads(std::min(2, omp_get_num_procs()));
#endif

            using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
            ThreadManager::init();
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if (mpi_rank_ != 0 || mpi_size_ < 2 || !this->output_files_) {
                return;
            }

            detail::mergeParallelLogFiles(eclState().getIOConfig().getOutputDir(),
                                          EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName),
                                          EWOMS_GET_PARAM(TypeTag, bool, EnableLoggingFalloutWarning));
        }

        void setupEbosSimulator()
        {
            ebosSimulator_ = std::make_unique<EbosSimulator>(EclGenericVanguard::comm(), /*verbose=*/false);
            ebosSimulator_->executionTimer().start();
            ebosSimulator_->model().applyInitialSolution();

            try {
                // Possible to force initialization only behavior (NOSIM).
                const std::string& dryRunString = EWOMS_GET_PARAM(TypeTag, std::string, EnableDryRun);
                if (dryRunString != "" && dryRunString != "auto") {
                    bool yesno;
                    if (dryRunString == "true"
                        || dryRunString == "t"
                        || dryRunString == "1")
                        yesno = true;
                    else if (dryRunString == "false"
                             || dryRunString == "f"
                             || dryRunString == "0")
                        yesno = false;
                    else
                        throw std::invalid_argument("Invalid value for parameter EnableDryRun: '"
                                                    +dryRunString+"'");
                    auto& ioConfig = eclState().getIOConfig();
                    ioConfig.overrideNOSIM(yesno);
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object" << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }
        }

        const Deck& deck() const
        { return ebosSimulator_->vanguard().deck(); }

        Deck& deck()
        { return ebosSimulator_->vanguard().deck(); }

        const EclipseState& eclState() const
        { return ebosSimulator_->vanguard().eclState(); }

        EclipseState& eclState()
        { return ebosSimulator_->vanguard().eclState(); }

        const Schedule& schedule() const
        { return ebosSimulator_->vanguard().schedule(); }

        // Run the simulator.
        int runSimulator()
        {
            return runSimulatorInitOrRun_(&FlowMainEbos::runSimulatorRunCallback_);
        }

        int runSimulatorInit()
        {
            return runSimulatorInitOrRun_(&FlowMainEbos::runSimulatorInitCallback_);
        }

    private:
        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorRunCallback_()
        {
            SimulatorReport report = simulator_->run(*simtimer_);
            runSimulatorAfterSim_(report);
            return report.success.exit_status;
        }

        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorInitCallback_()
        {
            simulator_->init(*simtimer_);
            return EXIT_SUCCESS;
        }

        // Output summary after simulation has completed
        void runSimulatorAfterSim_(SimulatorReport &report)
        {
            if (! this->output_cout_) {
                return;
            }

            const int threads
#if !defined(_OPENMP) || !_OPENMP
                = 1;
#else
                = omp_get_max_threads();
#endif

            printFlowTrailer(mpi_size_, threads, report);

            detail::handleExtraConvergenceOutput(report,
                                                 EWOMS_GET_PARAM(TypeTag, std::string, OutputExtraConvergenceInfo),
                                                 R"(OutputExtraConvergenceInfo (--output-extra-convergence-info))",
                                                 eclState().getIOConfig().getOutputDir(),
                                                 eclState().getIOConfig().getBaseName());
        }

        // Run the simulator.
        int runSimulatorInitOrRun_(int (FlowMainEbos::* initOrRunFunc)())
        {

            const auto& schedule = this->schedule();
            auto& ioConfig = eclState().getIOConfig();
            simtimer_ = std::make_unique<SimulatorTimer>();

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer_->init(schedule, (size_t)initConfig.getRestartStep());

            if (this->output_cout_) {
                std::ostringstream oss;

                // This allows a user to catch typos and misunderstandings in the
                // use of simulator parameters.
                if (Parameters::printUnused<TypeTag>(oss)) {
                    std::cout << "-----------------   Unrecognized parameters:   -----------------\n";
                    std::cout << oss.str();
                    std::cout << "----------------------------------------------------------------" << std::endl;
                }
            }

            if (!ioConfig.initOnly()) {
                if (this->output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                return (this->*initOrRunFunc)();
            }
            else {
                if (this->output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }
                return EXIT_SUCCESS;
            }
        }

    protected:

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_ = std::make_unique<Simulator>(*ebosSimulator_);
        }

        Grid& grid()
        { return ebosSimulator_->vanguard().grid(); }

    private:
        std::unique_ptr<EbosSimulator> ebosSimulator_;
        int  mpi_rank_ = 0;
        int  mpi_size_ = 1;
        std::any parallel_information_;
        std::unique_ptr<Simulator> simulator_;
        std::unique_ptr<SimulatorTimer> simtimer_;
        int argc_;
        char **argv_;
        bool output_cout_;
        bool output_files_;
    };

} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
