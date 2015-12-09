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
            using namespace Opm;

            // Must ensure an instance of the helper is created to initialise MPI.
            // For a build without MPI the Dune::FakeMPIHelper is used, so rank will
            // be 0 and size 1.
            const Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc, argv);
            const int mpi_rank = mpi_helper.rank();
            const int mpi_size = mpi_helper.size();

            // Write parameters used for later reference. (only if rank is zero)
            const bool output_cout = ( mpi_rank == 0 );

            if (output_cout)
                {
                    std::string version = moduleVersionName();
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

#ifdef _OPENMP
            if (!getenv("OMP_NUM_THREADS")) {
                //Default to at most 4 threads, regardless of 
                //number of cores (unless ENV(OMP_NUM_THREADS) is defined)
                int num_cores = omp_get_num_procs();
                int num_threads = std::min(4, num_cores);
                omp_set_num_threads(num_threads);
            }
#pragma omp parallel
            if (omp_get_thread_num() == 0){
                //opm_get_num_threads() only works as expected within a parallel region.
                std::cout << "OpenMP using " << omp_get_num_threads() << " threads." << std::endl;
            }
#endif

            // Read parameters, see if a deck was specified on the command line.
            if ( output_cout )
                {
                    std::cout << "---------------    Reading parameters     ---------------" << std::endl;
                }

            parameter::ParameterGroup param(argc, argv, false, output_cout);
            if( !output_cout )
                {
                    param.disableOutput();
                }

            if (!param.unhandledArguments().empty()) {
                if (param.unhandledArguments().size() != 1) {
                    std::cerr << "You can only specify a single input deck on the command line.\n";
                    return EXIT_FAILURE;
                } else {
                    param.insertParameter("deck_filename", param.unhandledArguments()[0]);
                }
            }

            // We must have an input deck. Grid and props will be read from that.
            if (!param.has("deck_filename")) {
                std::cerr << "This program must be run with an input deck.\n"
                    "Specify the deck filename either\n"
                    "    a) as a command line argument by itself\n"
                    "    b) as a command line parameter with the syntax deck_filename=<path to your deck>, or\n"
                    "    c) as a parameter in a parameter file (.param or .xml) passed to the program.\n";
                return EXIT_FAILURE;
            }

            // bool check_well_controls = false;
            // int max_well_control_iterations = 0;
            double gravity[3] = { 0.0 };
            std::string deck_filename = param.get<std::string>("deck_filename");

            // Write parameters used for later reference. (only if rank is zero)
            bool output = ( mpi_rank == 0 ) && param.getDefault("output", true);
            std::string output_dir;
            if (output) {
                // Create output directory if needed.
                output_dir =
                    param.getDefault("output_dir", std::string("output"));
                boost::filesystem::path fpath(output_dir);
                try {
                    create_directories(fpath);
                }
                catch (...) {
                    std::cerr << "Creating directories failed: " << fpath << std::endl;
                    return EXIT_FAILURE;
                }
                // Write simulation parameters.
                param.writeParam(output_dir + "/simulation.param");
            }

            std::string logFile = output_dir + "/LOGFILE.txt";
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
            if (param.has("output_interval")) {
                int output_interval = param.get<int>("output_interval");
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
                                              param);

            BlackoilPropsAdFromDeck new_props( deck, eclipseState, materialLawManager, grid );
            // check_well_controls = param.getDefault("check_well_controls", false);
            // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
            // Rock compressibility.

            RockCompressibility rock_comp(deck, eclipseState);

            // Gravity.
            gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

            typename Simulator::ReservoirState state;
            // Init state variables (saturation and pressure).
            if (param.has("init_saturation")) {
                initStateBasic(Opm::UgGridHelpers::numCells(grid),
                               Opm::UgGridHelpers::globalCell(grid),
                               Opm::UgGridHelpers::cartDims(grid),
                               Opm::UgGridHelpers::numFaces(grid),
                               Opm::UgGridHelpers::faceCells(grid),
                               Opm::UgGridHelpers::beginFaceCentroids(grid),
                               Opm::UgGridHelpers::beginCellCentroids(grid),
                               Opm::UgGridHelpers::dimensions(grid),
                               props, param, gravity[2], state);

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
                const double grav = param.getDefault("gravity", unit::gravity);
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

            const bool use_local_perm = param.getDefault("use_local_perm", true);

            DerivedGeology geoprops(grid, new_props, eclipseState, use_local_perm, grav);
            boost::any parallel_information;

            // At this point all properties and state variables are correctly initialized
            // If there are more than one processors involved, we now repartition the grid
            // and initilialize new properties and states for it.
            if( mpi_size > 1 )
                {
                    Opm::distributeGridAndData( grid, deck, eclipseState, state, new_props, geoprops, materialLawManager, parallel_information, use_local_perm );
                }

            // create output writer after grid is distributed, otherwise the parallel output
            // won't work correctly since we need to create a mapping from the distributed to
            // the global view
            Opm::BlackoilOutputWriter outputWriter(grid, param, eclipseState, pu, new_props.permeability() );

            // Solver for Newton iterations.
            std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver;
            {
                const std::string cprSolver = "cpr";
                const std::string interleavedSolver = "interleaved";
                const std::string directSolver = "direct";
                const std::string flowDefaultSolver = interleavedSolver;

                std::shared_ptr<const Opm::SimulationConfig> simCfg = eclipseState->getSimulationConfig();
                std::string solver_approach = flowDefaultSolver;

                if (param.has("solver_approach")) {
                    solver_approach = param.get<std::string>("solver_approach");
                }  else {
                    if (simCfg->useCPR()) {
                        solver_approach = cprSolver;
                    }
                }

                if (solver_approach == cprSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilCPR(param, parallel_information));
                } else if (solver_approach == interleavedSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilInterleaved(param, parallel_information));
                } else if (solver_approach == directSolver) {
                    fis_solver.reset(new NewtonIterationBlackoilSimple(param, parallel_information));
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

            Simulator simulator(param,
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
                if( output_cout )
                    {
                        std::cout << "\n\n================ Starting main simulation loop ===============\n"
                                  << std::flush;
                    }

                SimulatorReport fullReport = simulator.run(simtimer, state);

                if( output_cout )
                    {
                        std::cout << "\n\n================    End of simulation     ===============\n\n";
                        fullReport.reportFullyImplicit(std::cout);
                    }

                if (output) {
                    std::string filename = output_dir + "/walltime.txt";
                    std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
                    fullReport.reportParam(tot_os);
                    warnIfUnusedParams(param);
                }
            } else {
                outputWriter.writeInit( simtimer );
                if ( output_cout )
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

    }; // class FlowMain


} // namespace Opm

#endif // OPM_FLOWMAIN_HEADER_INCLUDED
