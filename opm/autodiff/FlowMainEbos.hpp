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


#include <sys/utsname.h>


#include <opm/simulators/ParallelFileMerger.hpp>
#include <opm/simulators/ensureDirectoryExists.hpp>

#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>
#include <opm/autodiff/MissingFeatures.hpp>
#include <opm/autodiff/moduleVersion.hpp>
#include <opm/autodiff/ExtractParallelGridInformationToISTL.hpp>
#include <opm/autodiff/RedistributeDataHandles.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Opm
{
    // The FlowMain class is the ebos based black-oil simulator.
    template <class TypeTag>
    class FlowMainEbos
    {
        enum FileOutputValue{
            //! \brief No output to files.
            OUTPUT_NONE = 0,
            //! \brief Output only to log files, no eclipse output.
            OUTPUT_LOG_ONLY = 1,
            //! \brief Output to all files.
            OUTPUT_ALL = 3
        };

    public:
        typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager MaterialLawManager;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) EbosSimulator;
        typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) EbosThreadManager;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        typedef Opm::SimulatorFullyImplicitBlackoilEbos<TypeTag> Simulator;
        typedef typename Simulator::ReservoirState ReservoirState;
        typedef typename Simulator::OutputWriter OutputWriter;

        /// This is the main function of Flow.
        /// It runs a complete simulation, with the given grid and
        /// simulator classes, based on user command-line input.  The
        /// content of this function used to be in the main() function of
        /// flow.cpp.
        int execute(int argc, char** argv)
        {
            try {
                setupParallelism();
                printStartupMessage();
                const bool ok = setupParameters(argc, argv);
                if (!ok) {
                    return EXIT_FAILURE;
                }

                setupEbosSimulator();
                setupOutput();
                setupLogging();
                printPRTHeader();
                extractMessages();
                runDiagnostics();
                writeInit();
                setupOutputWriter();
                setupLinearSolver();
                createSimulator();

                // Run.
                auto ret =  runSimulator();

                mergeParallelLogFiles();

                return ret;
            }
            catch (const std::exception &e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if( output_cout_ )
                {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }

                return EXIT_FAILURE;
            }
        }

    protected:
        void setupParallelism()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized here.
#if HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
            int mpi_size;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
            mpi_rank_ = 0;
            const int mpi_size = 1;
#endif
            output_cout_ = ( mpi_rank_ == 0 );
            must_distribute_ = ( mpi_size > 1 );

#ifdef _OPENMP
            // OpenMP setup.
            if (!getenv("OMP_NUM_THREADS")) {
                // Default to at most 4 threads, regardless of
                // number of cores (unless ENV(OMP_NUM_THREADS) is defined)
                int num_cores = omp_get_num_procs();
                int num_threads = std::min(4, num_cores);
                omp_set_num_threads(num_threads);
            }
            // omp_get_num_threads() only works as expected within a parallel region.
            const int num_omp_threads = omp_get_max_threads();
            if (mpi_size == 1) {
                std::cout << "OpenMP using " << num_omp_threads << " threads." << std::endl;
            } else {
                std::cout << "OpenMP using " << num_omp_threads << " threads on MPI rank " << mpi_rank_ << "." << std::endl;
            }
#endif
        }

        // Print startup message if on output rank.
        void printStartupMessage()
        {

            if (output_cout_) {
                const int lineLen = 70;
                const std::string version = moduleVersionName();
                const std::string banner = "This is flow "+version;
                const int bannerPreLen = (lineLen - 2 - banner.size())/2;
                const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
                std::cout << "*             including solvent and polymer capabilities.            *\n";
                std::cout << "*          For more information, see http://opm-project.org          *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";
            }
        }

        // Read parameters, see if a deck was specified on the command line, and if
        // it was, insert it into parameters.
        // Writes to:
        //   param_
        // Returns true if ok, false if not.
        bool setupParameters(int argc, char** argv)
        {
            param_ = ParameterGroup(argc, argv, false, output_cout_);

            // See if a deck was specified on the command line.
            if (!param_.unhandledArguments().empty()) {
                if (param_.unhandledArguments().size() != 1) {
                    std::cerr << "You can only specify a single input deck on the command line.\n";
                    return false;
                } else {
                    const auto casename = this->simulationCaseName( param_.unhandledArguments()[ 0 ] );
                    param_.insertParameter("deck_filename", casename.string() );
                }
            }

            // We must have an input deck. Grid and props will be read from that.
            if (!param_.has("deck_filename")) {
                std::cerr << "This program must be run with an input deck.\n"
                    "Specify the deck filename either\n"
                    "    a) as a command line argument by itself\n"
                    "    b) as a command line parameter with the syntax deck_filename=<path to your deck>, or\n"
                    "    c) as a parameter in a parameter file (.param or .xml) passed to the program.\n";
                return false;
            }
            return true;
        }

        // Set output_to_files_ and set/create output dir. Write parameter file.
        // Writes to:
        //   output_to_files_
        //   output_dir_
        // Throws std::runtime_error if failed to create (if requested) output dir.
        void setupOutput()
        {
            const std::string output = param_.getDefault("output", std::string("all"));
            static std::map<std::string, FileOutputValue> string2OutputEnum =
                { {"none", OUTPUT_NONE },
                  {"false", OUTPUT_LOG_ONLY },
                  {"log", OUTPUT_LOG_ONLY },
                  {"all" , OUTPUT_ALL },
                  {"true" , OUTPUT_ALL }};
            auto converted = string2OutputEnum.find(output);
            if ( converted != string2OutputEnum.end() )
            {
                output_ = string2OutputEnum[output];
            }
            else
            {
                std::cerr << "Value " << output <<
                    " passed to option output was invalid. Using \"all\" instead."
                          << std::endl;
            }

            output_to_files_ = output_cout_ && output_ > OUTPUT_NONE;

            // Setup output directory.
            auto& ioConfig = eclState().getIOConfig();
            // Default output directory is the directory where the deck is found.
            const std::string default_output_dir = ioConfig.getOutputDir();
            output_dir_ = param_.getDefault("output_dir", default_output_dir);
            // Override output directory if user specified.
            ioConfig.setOutputDir(output_dir_);

            // Write parameters used for later reference. (only if rank is zero)
            if (output_to_files_) {
                // Create output directory if needed.
                ensureDirectoryExists(output_dir_);
                // Write simulation parameters.
                param_.writeParam(output_dir_ + "/simulation.param");
            }
        }

        // Setup OpmLog backend with output_dir.
        void setupLogging()
        {
            std::string deck_filename = param_.get<std::string>("deck_filename");
            // create logFile
            using boost::filesystem::path;
            path fpath(deck_filename);
            std::string baseName;
            std::ostringstream debugFileStream;
            std::ostringstream logFileStream;

            if (boost::to_upper_copy(path(fpath.extension()).string()) == ".DATA") {
                baseName = path(fpath.stem()).string();
            } else {
                baseName = path(fpath.filename()).string();
            }

            logFileStream << output_dir_ << "/" << baseName;
            debugFileStream << output_dir_ << "/" << "." << baseName;

            if ( must_distribute_ && mpi_rank_ != 0 )
            {
                // Added rank to log file for non-zero ranks.
                // This prevents message loss.
                debugFileStream << "."<< mpi_rank_;
                // If the following file appears then there is a bug.
                logFileStream << "." << mpi_rank_;
            }
            logFileStream << ".PRT";
            debugFileStream << ".DEBUG";

            logFile_ = logFileStream.str();

            if( output_ > OUTPUT_NONE)
            {
                std::shared_ptr<EclipsePRTLog> prtLog = std::make_shared<EclipsePRTLog>(logFile_ , Log::NoDebugMessageTypes, false, output_cout_);
                OpmLog::addBackend( "ECLIPSEPRTLOG" , prtLog );
                prtLog->setMessageLimiter(std::make_shared<MessageLimiter>());
                prtLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(false));
            }

            if( output_ >= OUTPUT_LOG_ONLY && !param_.getDefault("no_debug_log", false) )
            {
                std::string debugFile = debugFileStream.str();
                std::shared_ptr<StreamLog> debugLog = std::make_shared<EclipsePRTLog>(debugFile, Log::DefaultMessageTypes, false, output_cout_);
                OpmLog::addBackend( "DEBUGLOG" ,  debugLog);
            }

            std::shared_ptr<StreamLog> streamLog = std::make_shared<StreamLog>(std::cout, Log::StdoutMessageTypes);
            OpmLog::addBackend( "STREAMLOG", streamLog);
            const auto& msgLimits = schedule().getMessageLimits();
            const std::map<int64_t, int> limits = {{Log::MessageType::Note, msgLimits.getCommentPrintLimit(0)},
                                                   {Log::MessageType::Info, msgLimits.getMessagePrintLimit(0)},
                                                   {Log::MessageType::Warning, msgLimits.getWarningPrintLimit(0)},
                                                   {Log::MessageType::Error, msgLimits.getErrorPrintLimit(0)},
                                                   {Log::MessageType::Problem, msgLimits.getProblemPrintLimit(0)},
                                                   {Log::MessageType::Bug, msgLimits.getBugPrintLimit(0)}};
            streamLog->setMessageLimiter(std::make_shared<MessageLimiter>(10, limits));
            streamLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(true));

            if ( output_cout_ )
            {
            // Read Parameters.
                OpmLog::debug("\n---------------    Reading parameters     ---------------\n");
            }
        }

        void printPRTHeader()
        {
          // Print header for PRT file.
          if ( output_cout_ ) {
              const std::string version = moduleVersionName();
              const double megabyte = 1024 * 1024;
              unsigned num_cpu = std::thread::hardware_concurrency();
              struct utsname arch;
              const char* user = getlogin();
              time_t now = std::time(0);
              struct tm  tstruct;
              char      tmstr[80];
              tstruct = *localtime(&now);
              strftime(tmstr, sizeof(tmstr), "%d-%m-%Y at %X", &tstruct);
              const double mem_size = getTotalSystemMemory() / megabyte;
              std::ostringstream ss;
              ss << "\n\n\n";
              ss << " ########  #          ######   #           #\n";
              ss << " #         #         #      #   #         # \n";
              ss << " #####     #         #      #    #   #   #  \n";
              ss << " #         #         #      #     # # # #   \n";
              ss << " #         #######    ######       #   #    \n\n";
              ss << "Flow is a simulator for fully implicit three-phase black-oil flow,";
              ss << " and is part of OPM.\nFor more information visit: http://opm-project.org \n\n";
              ss << "Flow Version  =  " + version + "\n";
              if (uname(&arch) == 0) {
                 ss << "System        =  " << arch.nodename << " (Number of cores: " << num_cpu;
                 ss << ", RAM: " << std::fixed << std::setprecision (2) << mem_size << " MB) \n";
                 ss << "Architecture  =  " << arch.sysname << " " << arch.machine << " (Release: " << arch.release;
                 ss << ", Version: " << arch.version << " )\n";
                 }
              if (user) {
                 ss << "User          =  " << user << std::endl;
                 }
              ss << "Simulation started on " << tmstr << " hrs\n";
              OpmLog::note(ss.str());
            }
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if( mpi_rank_ != 0 || !must_distribute_ || !output_to_files_ )
            {
                return;
            }

            namespace fs = boost::filesystem;
            fs::path output_path(".");
            if ( param_.has("output_dir") )
            {
                output_path = fs::path(output_dir_);
            }

            fs::path deck_filename(param_.get<std::string>("deck_filename"));

            std::for_each(fs::directory_iterator(output_path),
                          fs::directory_iterator(),
                          detail::ParallelFileMerger(output_path, deck_filename.stem().string()));
        }

        void setupEbosSimulator()
        {
            std::vector<const char*> argv;

            argv.push_back("flow_ebos");

            std::string deckFileParam("--ecl-deck-file-name=");
            deckFileParam += param_.get<std::string>("deck_filename");
            argv.push_back(deckFileParam.c_str());

#if defined(_OPENMP)
            std::string numThreadsParam("--threads-per-process=");
            int numThreads = omp_get_max_threads();

            numThreadsParam += std::to_string(numThreads);
            argv.push_back(numThreadsParam.c_str());
#endif // defined(_OPENMP)

            EbosSimulator::registerParameters();
            Ewoms::setupParameters_<TypeTag>(argv.size(), &argv[0]);
            EbosThreadManager::init();
            ebosSimulator_.reset(new EbosSimulator(/*verbose=*/false));
            ebosSimulator_->model().applyInitialSolution();

            // Create a grid with a global view.
            globalGrid_.reset(new Grid(grid()));
            globalGrid_->switchToGlobalView();

            try {
                if (output_cout_) {
                    MissingFeatures::checkKeywords(deck());
                }

                // Possible to force initialization only behavior (NOSIM).
                if (param_.has("nosim")) {
                    const bool nosim = param_.get<bool>("nosim");
                    auto& ioConfig = eclState().getIOConfig();
                    ioConfig.overrideNOSIM( nosim );
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object. See logfile: " << logFile_ << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }

            // Possibly override IOConfig setting (from deck) for how often RESTART files should get written to disk (every N report step)
            if (param_.has("output_interval")) {
                const int output_interval = param_.get<int>("output_interval");
                eclState().getRestartConfig().overrideRestartWriteInterval( size_t( output_interval ) );
            }
        }

        const Deck& deck() const
        { return ebosSimulator_->gridManager().deck(); }

        Deck& deck()
        { return ebosSimulator_->gridManager().deck(); }

        const EclipseState& eclState() const
        { return ebosSimulator_->gridManager().eclState(); }

        EclipseState& eclState()
        { return ebosSimulator_->gridManager().eclState(); }

        const Schedule& schedule() const
        { return ebosSimulator_->gridManager().schedule(); }

        const SummaryConfig& summaryConfig() const
        { return ebosSimulator_->gridManager().summaryConfig(); }
  
        // Extract messages from parser.
        // Writes to:
        //    OpmLog singleton.
        void extractMessages()
        {
            if ( !output_cout_ )
            {
                return;
            }

            auto extractMessage = [this](const Message& msg) {
                auto log_type = this->convertMessageType(msg.mtype);
                const auto& location = msg.location;
                if (location) {
                    OpmLog::addMessage(log_type, Log::fileMessage(location.filename, location.lineno, msg.message));
                } else {
                    OpmLog::addMessage(log_type, msg.message);
                }
            };

            // Extract messages from Deck.
            for(const auto& msg : deck().getMessageContainer()) {
                extractMessage(msg);
            }

            // Extract messages from EclipseState.
            for (const auto& msg : eclState().getMessageContainer()) {
                extractMessage(msg);
            }
        }

        // Run diagnostics.
        // Writes to:
        //   OpmLog singleton.
        void runDiagnostics()
        {
            if( ! output_cout_ )
            {
                return;
            }

            // Run relperm diagnostics
            RelpermDiagnostics diagnostic;
            diagnostic.diagnosis(eclState(), deck(), this->grid());
        }

        void writeInit()
        {
            bool output      = ( output_ > OUTPUT_LOG_ONLY );
            bool output_ecl  = param_.getDefault("output_ecl", true);
            auto int_vectors  = computeCellRanks(output, output_ecl);

            if( output && output_ecl && grid().comm().rank() == 0 )
            {
                exportNncStructure_();

                const EclipseGrid& inputGrid = eclState().getInputGrid();
                eclIO_.reset(new EclipseIO(eclState(),
                                           UgGridHelpers::createEclipseGrid( this->globalGrid() , inputGrid ),
                                           schedule(),
                                           summaryConfig()));
                eclIO_->writeInitial(computeLegacySimProps_(), int_vectors, nnc_);
                Problem& problem = ebosProblem();
                problem.setEclIO(std::move(eclIO_));
            }
        }

        // Setup output writer.
        // Writes to:
        //   output_writer_
        void setupOutputWriter()
        {
            // create output writer after grid is distributed, otherwise the parallel output
            // won't work correctly since we need to create a mapping from the distributed to
            // the global view

            output_writer_.reset(new OutputWriter(*ebosSimulator_,
                                                   param_));

        }

        // Run the simulator.
        // Returns EXIT_SUCCESS if it does not throw.
        int runSimulator()
        {
            const auto& schedule = this->schedule();
            const auto& timeMap = schedule.getTimeMap();
            auto& ioConfig = eclState().getIOConfig();
            SimulatorTimer simtimer;

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer.init(timeMap, (size_t)initConfig.getRestartStep());

            if (!ioConfig.initOnly()) {
                if (output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                SimulatorReport successReport = simulator_->run(simtimer);
                SimulatorReport failureReport = simulator_->failureReport();

                if (output_cout_) {
                    std::ostringstream ss;
                    ss << "\n\n================    End of simulation     ===============\n\n";
                    successReport.reportFullyImplicit(ss, &failureReport);
                    OpmLog::info(ss.str());
                    if (param_.anyUnused()) {
                        // This allows a user to catch typos and misunderstandings in the
                        // use of simulator parameters.
                        std::cout << "--------------------   Unused parameters:   --------------------\n";
                        param_.displayUsage();
                        std::cout << "----------------------------------------------------------------" << std::endl;
                    }
                }

            } else {
                if (output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }

            }
            return EXIT_SUCCESS;
        }

        // Setup linear solver.
        // Writes to:
        //   fis_solver_
        void setupLinearSolver()
        {
            typedef typename BlackoilModelEbos<TypeTag> :: ISTLSolverType ISTLSolverType;

            extractParallelGridInformationToISTL(grid(), parallel_information_);
            fis_solver_.reset( new ISTLSolverType( param_, parallel_information_ ) );
        }

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_.reset(new Simulator(*ebosSimulator_,
                                           param_,
                                           *fis_solver_,
                                           FluidSystem::enableDissolvedGas(),
                                           FluidSystem::enableVaporizedOil(),
                                           *output_writer_));
        }

    private:
        boost::filesystem::path simulationCaseName( const std::string& casename ) {
            namespace fs = boost::filesystem;

            const auto exists = []( const fs::path& f ) -> bool {
                if( !fs::exists( f ) ) return false;

                if( fs::is_regular_file( f ) ) return true;

                return fs::is_symlink( f )
                && fs::is_regular_file( fs::read_symlink( f ) );
            };

            auto simcase = fs::path( casename );

            if( exists( simcase ) ) {
                return simcase;
            }

            for( const auto& ext : { std::string("data"), std::string("DATA") } ) {
                if( exists( simcase.replace_extension( ext ) ) ) {
                    return simcase;
                }
            }

            throw std::invalid_argument( "Cannot find input case " + casename );
        }

        unsigned long long getTotalSystemMemory()
        {
            long pages = sysconf(_SC_PHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            return pages * page_size;
        }

        int64_t convertMessageType(const Message::type& mtype)
        {
            switch (mtype) {
            case Message::type::Debug:
                return Log::MessageType::Debug;
            case Message::type::Info:
                return Log::MessageType::Info;
            case Message::type::Warning:
                return Log::MessageType::Warning;
            case Message::type::Error:
                return Log::MessageType::Error;
            case Message::type::Problem:
                return Log::MessageType::Problem;
            case Message::type::Bug:
                return Log::MessageType::Bug;
            case Message::type::Note:
                return Log::MessageType::Note;
            }
            throw std::logic_error("Invalid messages type!\n");
        }

        Grid& grid()
        { return ebosSimulator_->gridManager().grid(); }

        const Grid& globalGrid()
        { return *globalGrid_; }

        Problem& ebosProblem()
        { return ebosSimulator_->problem(); }

        const Problem& ebosProblem() const
        { return ebosSimulator_->problem(); }

        std::shared_ptr<MaterialLawManager> materialLawManager()
        { return ebosProblem().materialLawManager(); }

        Scalar gravity() const
        { return ebosProblem().gravity()[2]; }

        std::map<std::string, std::vector<int> > computeCellRanks(bool output, bool output_ecl)
        {
            std::map<std::string, std::vector<int> > integerVectors;

            if(  output && output_ecl && grid().comm().size() > 1 )
            {
                typedef typename Grid::LeafGridView GridView;
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
                using ElementMapper =  Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
#else
                // Get the owner rank number for each cell
                using ElementMapper =  Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout>;
#endif
                using Handle = CellOwnerDataHandle<ElementMapper>;
                const Grid& globalGrid = this->globalGrid();
                const auto& globalGridView = globalGrid.leafGridView();
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
                ElementMapper globalMapper(globalGridView, Dune::mcmgElementLayout());
#else
                ElementMapper globalMapper(globalGridView);
#endif
                const auto globalSize = globalGrid.size(0);
                std::vector<int> ranks(globalSize, -1);
                Handle handle(globalMapper, ranks);
                this->grid().gatherData(handle);
                integerVectors.emplace("MPI_RANK", ranks);
            }

            return integerVectors;
        }

        data::Solution computeLegacySimProps_()
        {
            const int* dims = UgGridHelpers::cartDims(grid());
            const int globalSize = dims[0]*dims[1]*dims[2];

            data::CellData tranx = {UnitSystem::measure::transmissibility, std::vector<double>( globalSize ), data::TargetType::INIT};
            data::CellData trany = {UnitSystem::measure::transmissibility, std::vector<double>( globalSize ), data::TargetType::INIT};
            data::CellData tranz = {UnitSystem::measure::transmissibility, std::vector<double>( globalSize ), data::TargetType::INIT};

            for (size_t i = 0; i < tranx.data.size(); ++i) {
                tranx.data[0] = 0.0;
                trany.data[0] = 0.0;
                tranz.data[0] = 0.0;
            }

            const Grid& globalGrid = this->globalGrid();
            const auto& globalGridView = globalGrid.leafGridView();
            typedef typename Grid::LeafGridView GridView;
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
            ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());
#else
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
            ElementMapper globalElemMapper(globalGridView);
#endif
            const auto& cartesianCellIdx = globalGrid.globalCell();

            const auto* globalTrans = &(ebosSimulator_->gridManager().globalTransmissibility());
            if (grid().comm().size() < 2) {
                // in the sequential case we must use the transmissibilites defined by
                // the problem. (because in the sequential case, the grid manager does
                // not compute "global" transmissibilities for performance reasons. in
                // the parallel case, the problem's transmissibilities can't be used
                // because this object refers to the distributed grid and we need the
                // sequential version here.)
                globalTrans = &ebosSimulator_->problem().eclTransmissibilities();
            }

            auto elemIt = globalGridView.template begin</*codim=*/0>();
            const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++ elemIt) {
                const auto& elem = *elemIt;

                auto isIt = globalGridView.ibegin(elem);
                const auto& isEndIt = globalGridView.iend(elem);
                for (; isIt != isEndIt; ++ isIt) {
                    const auto& is = *isIt;

                    if (!is.neighbor())
                    {
                        continue; // intersection is on the domain boundary
                    }

                    unsigned c1 = globalElemMapper.index(is.inside());
                    unsigned c2 = globalElemMapper.index(is.outside());

                    if (c1 > c2)
                    {
                        continue; // we only need to handle each connection once, thank you.
                    }


                    int gc1 = std::min(cartesianCellIdx[c1], cartesianCellIdx[c2]);
                    int gc2 = std::max(cartesianCellIdx[c1], cartesianCellIdx[c2]);
                    if (gc2 - gc1 == 1) {
                        tranx.data[gc1] = globalTrans->transmissibility(c1, c2);
                    }

                    if (gc2 - gc1 == dims[0]) {
                        trany.data[gc1] = globalTrans->transmissibility(c1, c2);
                    }

                    if (gc2 - gc1 == dims[0]*dims[1]) {
                        tranz.data[gc1] = globalTrans->transmissibility(c1, c2);
                    }
                }
            }

            return {{"TRANX" , tranx},
                    {"TRANY" , trany} ,
                    {"TRANZ" , tranz}};
        }

        void exportNncStructure_()
        {
            nnc_ = eclState().getInputNNC();
            int nx = eclState().getInputGrid().getNX();
            int ny = eclState().getInputGrid().getNY();
            //int nz = eclState().getInputGrid().getNZ()

            const Grid& globalGrid = this->globalGrid();
            const auto& globalGridView = globalGrid.leafGridView();
            typedef typename Grid::LeafGridView GridView;
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
            ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());
#else
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
            ElementMapper globalElemMapper(globalGridView);
#endif

            const auto* globalTrans = &(ebosSimulator_->gridManager().globalTransmissibility());
            if (grid().comm().size() < 2) {
                // in the sequential case we must use the transmissibilites defined by
                // the problem. (because in the sequential case, the grid manager does
                // not compute "global" transmissibilities for performance reasons. in
                // the parallel case, the problem's transmissibilities can't be used
                // because this object refers to the distributed grid and we need the
                // sequential version here.)
                globalTrans = &ebosSimulator_->problem().eclTransmissibilities();
            }

            auto elemIt = globalGridView.template begin</*codim=*/0>();
            const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++ elemIt) {
                const auto& elem = *elemIt;

                auto isIt = globalGridView.ibegin(elem);
                const auto& isEndIt = globalGridView.iend(elem);
                for (; isIt != isEndIt; ++ isIt) {
                    const auto& is = *isIt;

                    if (!is.neighbor())
                    {
                        continue; // intersection is on the domain boundary
                    }

                    unsigned c1 = globalElemMapper.index(is.inside());
                    unsigned c2 = globalElemMapper.index(is.outside());

                    if (c1 > c2)
                    {
                        continue; // we only need to handle each connection once, thank you.
                    }

                    // TODO (?): use the cartesian index mapper to make this code work
                    // with grids other than Dune::CpGrid. The problem is that we need
                    // the a mapper for the sequential grid, not for the distributed one.
                    int cc1 = globalGrid.globalCell()[c1];
                    int cc2 = globalGrid.globalCell()[c2];

                    if (std::abs(cc1 - cc2) != 1 &&
                        std::abs(cc1 - cc2) != nx &&
                        std::abs(cc1 - cc2) != nx*ny)
                    {
                        nnc_.addNNC(cc1, cc2, globalTrans->transmissibility(c1, c2));
                    }
                }
            }
        }

        /// Convert saturations from a vector of individual phase saturation vectors
        /// to an interleaved format where all values for a given cell come before all
        /// values for the next cell, all in a single vector.
        template <class FluidSystem>
        void convertSats(std::vector<double>& sat_interleaved, const std::vector< std::vector<double> >& sat, const PhaseUsage& pu)
        {
            assert(sat.size() == 3);
            const auto nc = sat[0].size();
            const auto np = sat_interleaved.size() / nc;
            for (size_t c = 0; c < nc; ++c) {
                if ( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    const int opos = pu.phase_pos[BlackoilPhases::Liquid];
                    const std::vector<double>& sat_p = sat[ FluidSystem::oilPhaseIdx];
                    sat_interleaved[np*c + opos] = sat_p[c];
                }
                if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
                    const std::vector<double>& sat_p = sat[ FluidSystem::waterPhaseIdx];
                    sat_interleaved[np*c + wpos] = sat_p[c];
                }
                if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                    const std::vector<double>& sat_p = sat[ FluidSystem::gasPhaseIdx];
                    sat_interleaved[np*c + gpos] = sat_p[c];
                }
            }
        }


        std::unique_ptr<EbosSimulator> ebosSimulator_;
        int  mpi_rank_ = 0;
        bool output_cout_ = false;
        FileOutputValue output_ = OUTPUT_ALL;
        bool must_distribute_ = false;
        ParameterGroup param_;
        bool output_to_files_ = false;
        std::string output_dir_ = std::string(".");
        NNC nnc_;
        std::unique_ptr<EclipseIO> eclIO_;
        std::unique_ptr<OutputWriter> output_writer_;
        boost::any parallel_information_;
        std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver_;
        std::unique_ptr<Simulator> simulator_;
        std::string logFile_;
        // Needs to be shared pointer because it gets initialzed before MPI_Init.
        std::shared_ptr<Grid> globalGrid_;
    };
} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
