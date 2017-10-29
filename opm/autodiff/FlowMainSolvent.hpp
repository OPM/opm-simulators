/*
  Copyright 2015 IRIS AS
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FLOWMAINSOLVENT_HEADER_INCLUDED
#define OPM_FLOWMAINSOLVENT_HEADER_INCLUDED



#include <opm/autodiff/FlowMain.hpp>
#include <opm/autodiff/SolventPropsAdFromDeck.hpp>



namespace Opm
{

    // The FlowMainSolvent class is for a black-oil simulator with solvent.
    template <class Grid, class Simulator>
    class FlowMainSolvent : public FlowMainBase<FlowMainSolvent<Grid, Simulator>, Grid, Simulator>
    {
    protected:
        using Base = FlowMainBase<FlowMainSolvent<Grid, Simulator>, Grid, Simulator>;
        friend Base;

        // Set in setupGridAndProps()
        std::unique_ptr<SolventPropsAdFromDeck> solvent_props_;

        // ------------   Methods   ------------


        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (Base::output_cout_) {
                const std::string version = moduleVersionName();
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*                   This is Flow-Solvent (version " << version << ")"
                          << std::string(18 - version.size(), ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*     Flow-Solvent is a simulator for fully implicit three-phase,    *\n";
                std::cout << "*    four-component (black-oil + solvent) flow, and is part of OPM.  *\n";
                std::cout << "*           For more information see http://opm-project.org          *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";
            }
        }





        // Set up grid and property objects, by calling base class
        // version and then creating solvent property object.
        void setupGridAndProps()
        {
            Base::setupGridAndProps();

            const Grid& grid = Base::grid_init_->grid();
            solvent_props_.reset(new SolventPropsAdFromDeck(*Base::deck_,
                                                            *Base::eclipse_state_,
                                                            UgGridHelpers::numCells(grid),
                                                            UgGridHelpers::globalCell(grid)));
        }





        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            Base::simulator_.reset(new Simulator(Base::param_,
                                                 Base::grid_init_->grid(),
                                                 *Base::geoprops_,
                                                 *Base::fluidprops_,
                                                 *solvent_props_,
                                                 Base::rock_comp_->isActive() ? Base::rock_comp_.get() : nullptr,
                                                 *Base::fis_solver_,
                                                 Base::gravity_.data(),
                                                 Base::deck_->hasKeyword("DISGAS"),
                                                 Base::deck_->hasKeyword("VAPOIL"),
                                                 Base::eclipse_state_,
                                                 Base::schedule_,
                                                 Base::summary_config_,
                                                 *Base::output_writer_,
                                                 Base::deck_,
                                                 Base::threshold_pressures_,
                                                 Base::deck_->hasKeyword("SOLVENT")));
        }

        void setupLinearSolver()
        {
            // require_full_sparsity_pattern as default for solvent runs
            if (Base::deck_->hasKeyword("SOLVENT") && !Base::param_.has("require_full_sparsity_pattern") )  {
                Base::param_.insertParameter("require_full_sparsity_pattern","true");
            }
            Base::setupLinearSolver();
        }

    };


} // namespace Opm


#endif // OPM_FLOWMAINSOLVENT_HEADER_INCLUDED
