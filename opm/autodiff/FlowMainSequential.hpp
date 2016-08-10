/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FLOWMAINSEQUENTIAL_HEADER_INCLUDED
#define OPM_FLOWMAINSEQUENTIAL_HEADER_INCLUDED



#include <opm/autodiff/FlowMain.hpp>



namespace Opm
{

    // The FlowMainSequential class is for a black-oil simulator using the sequential models.
    template <class Grid, class Simulator>
    class FlowMainSequential : public FlowMainBase<FlowMainSequential<Grid, Simulator>, Grid, Simulator>
    {
    protected:
        using Base = FlowMainBase<FlowMainSequential<Grid, Simulator>, Grid, Simulator>;
        using Base::eclipse_state_;
        using Base::param_;
        using Base::fis_solver_;
        using Base::parallel_information_;
        friend Base;

        // ------------   Methods   ------------


        // Print startup message if on output rank.
        void printStartupMessage()
        {
            if (Base::output_cout_) {
                const std::string version = moduleVersionName();
                std::cout << "**********************************************************************\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*                 This is Flow-Sequential (version " << version << ")"
                          << std::string(17 - version.size(), ' ') << "*\n";
                std::cout << "*                                                                    *\n";
                std::cout << "*    Flow-Sequential is a simulator for fully implicit three-phase,  *\n";
                std::cout << "*                  black-oil flow, and is part of OPM.               *\n";
                std::cout << "*           For more information see http://opm-project.org          *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";

            }
        }





        // Setup linear solver.
        // Writes to:
        //   fis_solver_
        //   param_ (conditionally)
        // The CPR solver cannot be used with the sequential model.
        // Also, the interleaved solver requires the full sparsity pattern option.
        void setupLinearSolver()
        {
            const std::string cprSolver = "cpr";
            const std::string interleavedSolver = "interleaved";
            const std::string directSolver = "direct";
            std::string flowDefaultSolver = interleavedSolver;

            if (!param_.has("solver_approach")) {
                if (eclipse_state_->getSimulationConfig().useCPR()) {
                    flowDefaultSolver = cprSolver;
                }
            }

            const std::string solver_approach = param_.getDefault("solver_approach", flowDefaultSolver);

            if (solver_approach == cprSolver) {
                OPM_THROW( std::runtime_error , "CPR solver is not ready for use with sequential simulator.");
            } else if (solver_approach == interleavedSolver) {
                if (!param_.has("require_full_sparsity_pattern")) {
                    param_.insertParameter("require_full_sparsity_pattern", "true");
                }
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
            // We must override the min_iter argument unless it was already supplied, to avoid requiring iteration.
            if (!param_.has("min_iter")) {
                param_.insertParameter("min_iter", "0");
            }

            // Create the simulator instance.
            Base::simulator_.reset(new Simulator(Base::param_,
                                                 Base::grid_init_->grid(),
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


    };


} // namespace Opm


#endif // OPM_FLOWMAINSEQUENTIAL_HEADER_INCLUDED
