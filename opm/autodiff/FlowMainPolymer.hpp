/*
  Copyright 2014, 2015 STATOIL ASA.
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

#ifndef OPM_FLOWMAINPOLYMER_HEADER_INCLUDED
#define OPM_FLOWMAINPOLYMER_HEADER_INCLUDED



#include <opm/autodiff/FlowMain.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>



namespace Opm
{

    // The FlowMainPolymer class is for a black-oil simulator with polymer.
    template <class Grid, class Simulator>
    class FlowMainPolymer : public FlowMainBase<FlowMainPolymer<Grid, Simulator>, Grid, Simulator>
    {
    protected:
        using Base = FlowMainBase<FlowMainPolymer<Grid, Simulator>, Grid, Simulator>;
        using Base::eclipse_state_;
        using Base::param_;
        using Base::fis_solver_;
        using Base::parallel_information_;
        friend Base;

        // Set in setupGridAndProps()
        std::unique_ptr<PolymerProperties> polymer_props_legacy_; // Held by reference in polymer_props_
        std::unique_ptr<PolymerPropsAd> polymer_props_;

        // ------------   Methods   ------------


        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (Base::output_cout_) {
                const std::string version = moduleVersionName();
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*                   This is Flow-Polymer (version " << version << ")"
                          << std::string(18 - version.size(), ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*     Flow-Polymer is a simulator for fully implicit three-phase,    *\n";
                std::cout << "*    four-component (black-oil + polymer) flow, and is part of OPM.  *\n";
                std::cout << "*           For more information see http://opm-project.org          *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";

            }
        }





        // Set up grid and property objects, by calling base class
        // version and then creating polymer property objects.
        // Writes to:
        //   polymer_props_legacy_
        //   polymer_props_
        void setupGridAndProps()
        {
            Base::setupGridAndProps();

            if (Base::deck_->hasKeyword("POLYMER")) {
                polymer_props_legacy_.reset(new PolymerProperties(*Base::deck_, *Base::eclipse_state_));
                polymer_props_.reset(new PolymerPropsAd(*polymer_props_legacy_));
            }
        }

        // Setup linear solver.
        // Writes to:
        //   fis_solver_
        // Currently, the CPR solver is not ready for polymer solver yet
        void setupLinearSolver()
        {
            const std::string cprSolver = "cpr";
            const std::string interleavedSolver = "interleaved";
            const std::string directSolver = "direct";
            const std::string flowDefaultSolver = interleavedSolver;

            const Opm::SimulationConfig& simCfg = eclipse_state_->getSimulationConfig();
            std::string solver_approach = flowDefaultSolver;

            if (param_.has("solver_approach")) {
                solver_approach = param_.template get<std::string>("solver_approach");
            }  else {
                if (simCfg.useCPR()) {
                    solver_approach = cprSolver;
                }
            }

            if (solver_approach == cprSolver) {
                OPM_THROW( std::runtime_error , "CPR solver is not ready for use with polymer solver yet.");
            } else if (solver_approach == interleavedSolver) {
                fis_solver_.reset(new NewtonIterationBlackoilInterleaved(param_, parallel_information_));
            } else if (solver_approach == directSolver) {
                fis_solver_.reset(new NewtonIterationBlackoilSimple(param_, parallel_information_));
            } else {
                OPM_THROW( std::runtime_error , "Internal error - solver approach " << solver_approach << " not recognized.");
            }
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
                                                 *polymer_props_,
                                                 Base::rock_comp_->isActive() ? Base::rock_comp_.get() : nullptr,
                                                 *Base::fis_solver_,
                                                 Base::gravity_.data(),
                                                 Base::deck_->hasKeyword("DISGAS"),
                                                 Base::deck_->hasKeyword("VAPOIL"),
                                                 Base::deck_->hasKeyword("POLYMER"),
                                                 Base::deck_->hasKeyword("PLYSHLOG"),
                                                 Base::deck_->hasKeyword("SHRATE"),
                                                 Base::eclipse_state_,
                                                 *Base::output_writer_,
                                                 Base::deck_,
                                                 Base::threshold_pressures_));
        }


    };


} // namespace Opm




#endif // OPM_FLOWMAINPOLYMER_HEADER_INCLUDED
