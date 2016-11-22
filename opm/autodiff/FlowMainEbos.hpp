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

#include <opm/autodiff/FlowMain.hpp>
#include <opm/autodiff/BlackoilModelEbos.hpp>

#include <ewoms/version.hh>

namespace Opm
{
    // The FlowMain class is the ebos based black-oil simulator.
    class FlowMainEbos : public FlowMainBase<FlowMainEbos, Dune::CpGrid, Opm::SimulatorFullyImplicitBlackoilEbos>
    {
    protected:
        typedef Opm::SimulatorFullyImplicitBlackoilEbos Simulator;
        typedef FlowMainBase<FlowMainEbos, Dune::CpGrid, Simulator> Base;
        friend Base;

        typedef typename TTAG(EclFlowProblem) TypeTag;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) EbosSimulator;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (output_cout_) {
                const int lineLen = 70;
                const std::string version = moduleVersionName();
                const std::string banner = "This is flow_ebos (version "+version+")";
                const std::string ewomsVersion = "(eWoms version: " + Ewoms::versionString() + ")";
                const int bannerPreLen = (lineLen - 2 - banner.size())/2;
                const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
                const int eVPreLen = (lineLen - 2 - ewomsVersion.size())/2;
                const int eVPostLen = eVPreLen + (lineLen - 2 - ewomsVersion.size())%2;
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
                std::cout << "*" << std::string(eVPreLen, ' ') << ewomsVersion << std::string(eVPostLen, ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
                std::cout << "*            and is part of OPM. For more information see:           *\n";
                std::cout << "*                       http://opm-project.org                       *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";
            }
        }

        // Parser the input and creates the Deck and EclipseState objects.
        // Writes to:
        //   deck_
        //   eclipse_state_
        // May throw if errors are encountered, here configured to be somewhat tolerant.
        void readDeckInput()
        {
            std::string progName("flow_ebos");
            std::string deckFile("--ecl-deck-file-name=");
            deckFile += param_.get<std::string>("deck_filename");
            char* ptr[2];
            ptr[ 0 ] = const_cast< char * > (progName.c_str());
            ptr[ 1 ] = const_cast< char * > (deckFile.c_str());
            EbosSimulator::registerParameters();
            Ewoms::setupParameters_< TypeTag > ( 2, ptr );
            ebosSimulator_.reset(new EbosSimulator(/*verbose=*/false));
            ebosSimulator_->model().applyInitialSolution();

            Base::deck_ = ebosSimulator_->gridManager().deck();
            Base::eclipse_state_ = ebosSimulator_->gridManager().eclState();

            try {
                if (Base::output_cout_) {
                    MissingFeatures::checkKeywords(*Base::deck_);
                }

                IOConfig& ioConfig = Base::eclipse_state_->getIOConfig();
                ioConfig.setOutputDir(Base::output_dir_);

                // Possible to force initialization only behavior (NOSIM).
                if (Base::param_.has("nosim")) {
                    const bool nosim = Base::param_.get<bool>("nosim");
                    ioConfig.overrideNOSIM( nosim );
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object. See logfile: " << logFile_ << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }

            // Possibly override IOConfig setting (from deck) for how often RESTART files should get written to disk (every N report step)
            if (Base::param_.has("output_interval")) {
                const int output_interval = Base::param_.get<int>("output_interval");
                eclipse_state_->getRestartConfig().overrideRestartWriteInterval( size_t( output_interval ) );

            }
        }

        // Create grid and property objects.
        // Writes to:
        //   grid_init_
        //   material_law_manager_
        //   fluidprops_
        //   rock_comp_
        //   gravity_
        //   use_local_perm_
        //   geoprops_
        void setupGridAndProps()
        {
            Dune::CpGrid& grid = ebosSimulator_->gridManager().grid();

            //Base::material_law_manager_ = ebosSimulator_->problem().materialLawManager();

            grid.switchToGlobalView();
            // Create material law manager.
            std::vector<int> compressedToCartesianIdx;
            Opm::createGlobalCellArray(grid, compressedToCartesianIdx);
            material_law_manager_.reset(new MaterialLawManager());
            material_law_manager_->initFromDeck(*deck_, *eclipse_state_, compressedToCartesianIdx);



            grid_init_.reset(new GridInit<Grid>());
            grid_init_->setGrid(grid);

            // create the legacy properties objects
            Base::fluidprops_.reset(new BlackoilPropsAdFromDeck(*deck_,
                                                                *eclipse_state_,
                                                                Base::material_law_manager_,
                                                                grid));
            // Gravity.
            assert(UgGridHelpers::dimensions(grid) == 3);
            Base::gravity_.fill(0.0);
            Base::gravity_[2] = Base::deck_->hasKeyword("NOGRAV")
                ? Base::param_.getDefault("gravity", 0.0)
                : Base::param_.getDefault("gravity", unit::gravity);

            // Geological properties
            Base::use_local_perm_ = true;
            geoprops_.reset(new DerivedGeology(grid,
                                               *Base::fluidprops_,
                                               *Base::eclipse_state_,
                                               Base::use_local_perm_,
                                               Base::gravity_.data()));

        }

        // Initialise the reservoir state. Updated fluid props for SWATINIT.
        // Writes to:
        //   state_
        //   threshold_pressures_
        //   fluidprops_ (if SWATINIT is used)
        void setupState()
        {
            const PhaseUsage pu = Opm::phaseUsageFromDeck(*deck_);
            const Grid& grid = Base::grid_init_->grid();

            // Need old-style fluid object for init purposes (only).
            BlackoilPropertiesFromDeck props(*deck_,
                                             *Base::eclipse_state_,
                                             Base::material_law_manager_,
                                             grid.size(/*codim=*/0),
                                             grid.globalCell().data(),
                                             grid.logicalCartesianSize().data(),
                                             param_);


            // Init state variables (saturation and pressure).
            if (param_.has("init_saturation")) {
                state_.reset(new ReservoirState(grid.size(/*codim=*/0),
                                                grid.numFaces(),
                                                props.numPhases()));

                initStateBasic(grid.size(/*codim=*/0),
                               grid.globalCell().data(),
                               grid.logicalCartesianSize().data(),
                               grid.numFaces(),
                               Opm::UgGridHelpers::faceCells(grid),
                               Opm::UgGridHelpers::beginFaceCentroids(grid),
                               Opm::UgGridHelpers::beginCellCentroids(grid),
                               Grid::dimension,
                               props, param_, gravity_[2], *state_);

                initBlackoilSurfvol(Opm::UgGridHelpers::numCells(grid), props, *state_);

                enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
                if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
                    const int numPhases = props.numPhases();
                    const int numCells  = Opm::UgGridHelpers::numCells(grid);

                    // Uglyness 1: The state is a templated type, here we however make explicit use BlackoilState.
                    auto& gor = state_->getCellData( BlackoilState::GASOILRATIO );
                    const auto& surface_vol = state_->getCellData( BlackoilState::SURFACEVOL );
                    for (int c = 0; c < numCells; ++c) {
                        // Uglyness 2: Here we explicitly use the layout of the saturation in the surface_vol field.
                        gor[c] = surface_vol[ c * numPhases + pu.phase_pos[Gas]] / surface_vol[ c * numPhases + pu.phase_pos[Oil]];
                    }
                }
            } else if (deck_->hasKeyword("EQUIL")) {
                // Which state class are we really using - what a f... mess?
                state_.reset( new ReservoirState( Opm::UgGridHelpers::numCells(grid),
                                                  Opm::UgGridHelpers::numFaces(grid),
                                                  props.numPhases()));

                initStateEquil(grid, props, *deck_, *Base::eclipse_state_, gravity_[2], *state_);
                //state_.faceflux().resize(Opm::UgGridHelpers::numFaces(grid), 0.0);
            } else {
                state_.reset( new ReservoirState( Opm::UgGridHelpers::numCells(grid),
                                                  Opm::UgGridHelpers::numFaces(grid),
                                                  props.numPhases()));
                initBlackoilStateFromDeck(Opm::UgGridHelpers::numCells(grid),
                                          Opm::UgGridHelpers::globalCell(grid),
                                          Opm::UgGridHelpers::numFaces(grid),
                                          Opm::UgGridHelpers::faceCells(grid),
                                          Opm::UgGridHelpers::beginFaceCentroids(grid),
                                          Opm::UgGridHelpers::beginCellCentroids(grid),
                                          Opm::UgGridHelpers::dimensions(grid),
                                          props, *deck_, gravity_[2], *state_);
            }

            // Threshold pressures.
            std::map<std::pair<int, int>, double> maxDp;
            computeMaxDp(maxDp, *deck_, *Base::eclipse_state_, grid_init_->grid(), *state_, props, gravity_[2]);
            threshold_pressures_ = thresholdPressures(*deck_, *Base::eclipse_state_, grid, maxDp);
            std::vector<double> threshold_pressures_nnc = thresholdPressuresNNC(*Base::eclipse_state_, Base::eclipse_state_->getInputNNC(), maxDp);
            threshold_pressures_.insert(threshold_pressures_.end(), threshold_pressures_nnc.begin(), threshold_pressures_nnc.end());

            // The capillary pressure is scaled in fluidprops_ to match the scaled capillary pressure in props.
            if (deck_->hasKeyword("SWATINIT")) {
                const int numCells = Opm::UgGridHelpers::numCells(grid);
                std::vector<int> cells(numCells);
                for (int c = 0; c < numCells; ++c) { cells[c] = c; }
                std::vector<double> pc = state_->saturation();
                props.capPress(numCells, state_->saturation().data(), cells.data(), pc.data(), nullptr);
                fluidprops_->setSwatInitScaling(state_->saturation(), pc);
            }
            initHydroCarbonState(*state_, pu, Opm::UgGridHelpers::numCells(grid), deck_->hasKeyword("DISGAS"), deck_->hasKeyword("VAPOIL"));
        }

        // Setup linear solver.
        // Writes to:
        //   fis_solver_
        void setupLinearSolver()
        {
            typedef typename BlackoilModelEbos :: ISTLSolverType ISTLSolverType;
            Base::fis_solver_.reset( new ISTLSolverType( Base::param_, Base::parallel_information_ ) );
        }

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            Base::simulator_.reset(new Simulator(*ebosSimulator_,
                                                 Base::param_,
                                                 *Base::geoprops_,
                                                 *Base::fluidprops_,
                                                 /*rockComp=*/nullptr,
                                                 *Base::fis_solver_,
                                                 Base::gravity_.data(),
                                                 Base::deck_->hasKeyword("DISGAS"),
                                                 Base::deck_->hasKeyword("VAPOIL"),
                                                 Base::eclipse_state_,
                                                 *Base::output_writer_,
                                                 Base::threshold_pressures_));
        }

    private:
        std::unique_ptr<EbosSimulator> ebosSimulator_;
    };
} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
