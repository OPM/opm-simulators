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

#ifndef OPM_FLOWMAIN_HEADER_INCLUDED
#define OPM_FLOWMAIN_HEADER_INCLUDED


#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/common/utility/platform_dependent/reenable_warnings.h>


#include <opm/core/grid/GridManager.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/createGlobalCellArray.hpp>
#include <opm/autodiff/GridInit.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/initStateEquil.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/thresholdPressures.hpp> // Note: the GridHelpers must be included before this (to make overloads available). \TODO: Fix.

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/RedistributeDataHandles.hpp>
#include <opm/autodiff/moduleVersion.hpp>

#include <opm/core/utility/share_obj.hpp>

#include <opm/parser/eclipse/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/OpmLog/StreamLog.hpp>
#include <opm/parser/eclipse/OpmLog/CounterLog.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <memory>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <numeric>
#include <cstdlib>




namespace Opm
{

    /// Calling this will print the unused parameters, if any.
    /// This allows a user to catch typos and misunderstandings in the
    /// use of simulator parameters.
    inline void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }


    /// This is the main function of Flow.
    /// It runs a complete simulation, with the given grid and
    /// simulator classes, based on user command-line input.  The
    /// content of this function used to be in the main() function of
    /// flow.cpp.
    template <class Grid, class Simulator>
    class FlowMain
    {
    public:
        int execute(int argc, char** argv)
        try {
            setupParallelism(argc, argv);
            printStartupMessage();
            const bool ok = setupParameters(argc, argv);
            if (!ok) {
                return EXIT_FAILURE;
            }
            setupOutput();

            std::string logFile = output_dir_ + "/LOGFILE.txt";
            Opm::ParserPtr parser(new Opm::Parser());
            {
                std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(logFile , Opm::Log::DefaultMessageTypes);
                std::shared_ptr<Opm::CounterLog> counterLog = std::make_shared<Opm::CounterLog>(Opm::Log::DefaultMessageTypes);

                Opm::OpmLog::addBackend( "STREAM" , streamLog );
                Opm::OpmLog::addBackend( "COUNTER" , counterLog );
            }

            Opm::ParseMode parseMode({{ ParseMode::PARSE_RANDOM_SLASH , InputError::IGNORE }});
            Opm::DeckConstPtr deck;
            std::shared_ptr<EclipseState> eclipseState;
            std::string deck_filename = param_.get<std::string>("deck_filename");
            try {
                deck = parser->parseFile(deck_filename, parseMode);
                Opm::checkDeck(deck);
                eclipseState.reset(new Opm::EclipseState(deck , parseMode));
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid ECLIPSESTATE object. See logfile: " << logFile << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                return EXIT_FAILURE;
            }

            std::vector<double> porv = eclipseState->getDoubleGridProperty("PORV")->getData();
            GridInit<Grid> grid_init(deck, eclipseState, porv);
            auto&& grid = grid_init.grid();

            // Possibly override IOConfig setting (from deck) for how often RESTART files should get written to disk (every N report step)
            if (param_.has("output_interval")) {
                int output_interval = param_.get<int>("output_interval");
                IOConfigPtr ioConfig = eclipseState->getIOConfig();
                ioConfig->overrideRestartWriteInterval((size_t)output_interval);
            }

            const PhaseUsage pu = Opm::phaseUsageFromDeck(deck);

            std::vector<int> compressedToCartesianIdx;
            Opm::createGlobalCellArray(grid, compressedToCartesianIdx);

            typedef BlackoilPropsAdFromDeck::MaterialLawManager MaterialLawManager;
            auto materialLawManager = std::make_shared<MaterialLawManager>();
            materialLawManager->initFromDeck(deck, eclipseState, compressedToCartesianIdx);

            // Rock and fluid init
            BlackoilPropertiesFromDeck props( deck, eclipseState, materialLawManager,
                                              Opm::UgGridHelpers::numCells(grid),
                                              Opm::UgGridHelpers::globalCell(grid),
                                              Opm::UgGridHelpers::cartDims(grid),
                                              param_);

            BlackoilPropsAdFromDeck new_props( deck, eclipseState, materialLawManager, grid );
            // check_well_controls = param.getDefault("check_well_controls", false);
            // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
            // Rock compressibility.

            RockCompressibility rock_comp(deck, eclipseState);

            // Gravity.
            double gravity[3] = { 0.0 };
            gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

            typename Simulator::ReservoirState state;
            // Init state variables (saturation and pressure).
            if (param_.has("init_saturation")) {
                initStateBasic(Opm::UgGridHelpers::numCells(grid),
                               Opm::UgGridHelpers::globalCell(grid),
                               Opm::UgGridHelpers::cartDims(grid),
                               Opm::UgGridHelpers::numFaces(grid),
                               Opm::UgGridHelpers::faceCells(grid),
                               Opm::UgGridHelpers::beginFaceCentroids(grid),
                               Opm::UgGridHelpers::beginCellCentroids(grid),
                               Opm::UgGridHelpers::dimensions(grid),
                               props, param_, gravity[2], state);

                initBlackoilSurfvol(Opm::UgGridHelpers::numCells(grid), props, state);

                enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
                if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
                    const int numPhases = props.numPhases();
                    const int numCells  = Opm::UgGridHelpers::numCells(grid);
                    for (int c = 0; c < numCells; ++c) {
                        state.gasoilratio()[c] = state.surfacevol()[c*numPhases + pu.phase_pos[Gas]]
                            / state.surfacevol()[c*numPhases + pu.phase_pos[Oil]];
                    }
                }
            } else if (deck->hasKeyword("EQUIL") && props.numPhases() == 3) {
                state.init(Opm::UgGridHelpers::numCells(grid),
                           Opm::UgGridHelpers::numFaces(grid),
                           props.numPhases());
                const double grav = param_.getDefault("gravity", unit::gravity);
                initStateEquil(grid, props, deck, eclipseState, grav, state);
                state.faceflux().resize(Opm::UgGridHelpers::numFaces(grid), 0.0);
            } else {
                initBlackoilStateFromDeck(Opm::UgGridHelpers::numCells(grid),
                                          Opm::UgGridHelpers::globalCell(grid),
                                          Opm::UgGridHelpers::numFaces(grid),
                                          Opm::UgGridHelpers::faceCells(grid),
                                          Opm::UgGridHelpers::beginFaceCentroids(grid),
                                          Opm::UgGridHelpers::beginCellCentroids(grid),
                                          Opm::UgGridHelpers::dimensions(grid),
                                          props, deck, gravity[2], state);
            }


            // The capillary pressure is scaled in new_props to match the scaled capillary pressure in props.
            if (deck->hasKeyword("SWATINIT")) {
                const int numCells = Opm::UgGridHelpers::numCells(grid);
                std::vector<int> cells(numCells);
                for (int c = 0; c < numCells; ++c) { cells[c] = c; }
                std::vector<double> pc = state.saturation();
                props.capPress(numCells, state.saturation().data(), cells.data(), pc.data(),NULL);
                new_props.setSwatInitScaling(state.saturation(),pc);
            }

            bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
            const double *grav = use_gravity ? &gravity[0] : 0;

            const bool use_local_perm = param_.getDefault("use_local_perm", true);

            DerivedGeology geoprops(grid, new_props, eclipseState, use_local_perm, grav);
            boost::any parallel_information;

            // At this point all properties and state variables are correctly initialized
            // If there are more than one processors involved, we now repartition the grid
            // and initilialize new properties and states for it.
            if( must_distribute_ )
                {
                    Opm::distributeGridAndData( grid, deck, eclipseState, state, new_props, geoprops, materialLawManager, parallel_information, use_local_perm );
                }

            // create output writer after grid is distributed, otherwise the parallel output
            // won't work correctly since we need to create a mapping from the distributed to
            // the global view
            Opm::BlackoilOutputWriter outputWriter(grid, param_, eclipseState, pu, new_props.permeability() );

            // Solver for Newton iterations.
            std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver;
            {
                const std::string cprSolver = "cpr";
                const std::string interleavedSolver = "interleaved";
                const std::string directSolver = "direct";
                const std::string flowDefaultSolver = interleavedSolver;

                std::shared_ptr<const Opm::SimulationConfig> simCfg = eclipseState->getSimulationConfig();
                std::string solver_approach = flowDefaultSolver;

                if (param_.has("solver_approach")) {
                    solver_approach = param_.get<std::string>("solver_approach");
                }  else {
                    if (simCfg->useCPR()) {
                        solver_approach = cprSolver;
                    }
                }

                if (solver_approach == cprSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilCPR(param_, parallel_information));
                } else if (solver_approach == interleavedSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilInterleaved(param_, parallel_information));
                } else if (solver_approach == directSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilSimple(param_, parallel_information));
                } else {
                    OPM_THROW( std::runtime_error , "Internal error - solver approach " << solver_approach << " not recognized.");
                }

            }

            Opm::ScheduleConstPtr schedule = eclipseState->getSchedule();
            Opm::TimeMapConstPtr timeMap(schedule->getTimeMap());
            SimulatorTimer simtimer;

            // initialize variables
            simtimer.init(timeMap);

            std::map<std::pair<int, int>, double> maxDp;
            computeMaxDp(maxDp, deck, eclipseState, grid, state, props, gravity[2]);
            std::vector<double> threshold_pressures = thresholdPressures(deck, eclipseState, grid, maxDp);

            Simulator simulator(param_,
                                grid,
                                geoprops,
                                new_props,
                                rock_comp.isActive() ? &rock_comp : 0,
                                *fis_solver,
                                grav,
                                deck->hasKeyword("DISGAS"),
                                deck->hasKeyword("VAPOIL"),
                                eclipseState,
                                outputWriter,
                                threshold_pressures);

            if (!schedule->initOnly()){
                if( output_cout_ )
                    {
                        std::cout << "\n\n================ Starting main simulation loop ===============\n"
                                  << std::flush;
                    }

                SimulatorReport fullReport = simulator.run(simtimer, state);

                if( output_cout_ )
                    {
                        std::cout << "\n\n================    End of simulation     ===============\n\n";
                        fullReport.reportFullyImplicit(std::cout);
                    }

                if (output_to_files_) {
                    std::string filename = output_dir_ + "/walltime.txt";
                    std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
                    fullReport.reportParam(tot_os);
                    warnIfUnusedParams(param_);
                }
            } else {
                outputWriter.writeInit( simtimer );
                if ( output_cout_ )
                    {
                        std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                    }
            }
            return EXIT_SUCCESS;
        }
        catch (const std::exception &e) {
            std::cerr << "Program threw an exception: " << e.what() << "\n";
            return EXIT_FAILURE;
        }



    private:

        // ------------   Data members   ------------





        bool output_cout_ = false;
        bool must_distribute_ = false;
        parameter::ParameterGroup param_;
        bool output_to_files_ = false;
        std::string output_dir_ = "output";




        // ------------   Methods   ------------





        // Set up MPI and OpenMP.
        // Writes to:
        //   output_cout_
        //   must_distribute_
        void setupParallelism(int argc, char** argv)
        {
            // MPI setup.
            // Must ensure an instance of the helper is created to initialise MPI.
            // For a build without MPI the Dune::FakeMPIHelper is used, so rank will
            // be 0 and size 1.
            const Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc, argv);
            const int mpi_rank = mpi_helper.rank();
            const int mpi_size = mpi_helper.size();
            output_cout_ = ( mpi_rank == 0 );
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
#pragma omp parallel
            if (omp_get_thread_num() == 0) {
                // omp_get_num_threads() only works as expected within a parallel region.
                const int num_omp_threads = omp_get_num_threads();
                if (mpi_size == 1) {
                    std::cout << "OpenMP using " << num_omp_threads << " threads." << std::endl;
                } else {
                    std::cout << "OpenMP using " << num_omp_threads << " threads on MPI rank " << mpi_rank << "." << std::endl;
                }
            }
#endif
        }





        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (output_cout_) {
                const std::string version = moduleVersionName();
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*                   This is Flow (version " << version << ")"
                          << std::string(26 - version.size(), ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
                std::cout << "*            and is part of OPM. For more information see:           *\n";
                std::cout << "*                       http://opm-project.org                       *\n";
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
            // Read parameters.
            if ( output_cout_ ) {
                std::cout << "---------------    Reading parameters     ---------------" << std::endl;
            }
            param_ = parameter::ParameterGroup(argc, argv, false, output_cout_);

            // See if a deck was specified on the command line.
            if (!param_.unhandledArguments().empty()) {
                if (param_.unhandledArguments().size() != 1) {
                    std::cerr << "You can only specify a single input deck on the command line.\n";
                    return false;
                } else {
                    param_.insertParameter("deck_filename", param_.unhandledArguments()[0]);
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
            // Write parameters used for later reference. (only if rank is zero)
            output_to_files_ = output_cout_ && param_.getDefault("output", true);
            if (output_to_files_) {
                // Create output directory if needed.
                output_dir_ =
                    param_.getDefault("output_dir", std::string("output"));
                boost::filesystem::path fpath(output_dir_);
                try {
                    create_directories(fpath);
                }
                catch (...) {
                    OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
                }
                // Write simulation parameters.
                param_.writeParam(output_dir_ + "/simulation.param");
            }
        }


    }; // class FlowMain


} // namespace Opm

#endif // OPM_FLOWMAIN_HEADER_INCLUDED
