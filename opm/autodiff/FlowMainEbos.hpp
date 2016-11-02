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
            IOConfig& ioConfig = Base::eclipse_state_->getIOConfig();
            ioConfig.setOutputDir(Base::output_dir_);

            // Possibly override IOConfig setting (from deck) for how often RESTART files should get written to disk (every N report step)
            if (Base::param_.has("output_interval")) {
                const int output_interval = Base::param_.get<int>("output_interval");
                eclipse_state_->getRestartConfig().overrideRestartWriteInterval( size_t( output_interval ) );

            }

            // Possible to force initialization only behavior (NOSIM).
            if (Base::param_.has("nosim")) {
                const bool nosim = Base::param_.get<bool>("nosim");
                ioConfig.overrideNOSIM( nosim );
            }
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
                                                 Base::rock_comp_->isActive() ? Base::rock_comp_.get() : nullptr,
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
