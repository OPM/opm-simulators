// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef BLACK_OIL_MODEL_TPSA_HPP
#define BLACK_OIL_MODEL_TPSA_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/BlackoilModel.hpp>

#include <stdexcept>
#include <string>

#include <fmt/format.h>


namespace Opm {

/*!
* \brief Black oil model for coupling Flow simulations with TPSA geomechanics
*/
template <class TypeTag>
class BlackoilModelTPSA : public BlackoilModel<TypeTag>
{
    using ParentType = BlackoilModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    using ModelParameters = typename ParentType::ModelParameters;

    /*!
    * \brief Constructor
    *
    * \param simulator Reference to simulator object
    * \param param Reference to parameters for model
    * \param well_model Refenerence to well model
    * \param terminal_output Bool for terminal output
    */
    explicit BlackoilModelTPSA(Simulator& simulator,
                               const ModelParameters& param,
                               BlackoilWellModel<TypeTag>& well_model,
                               const bool terminal_output)
        : ParentType(simulator, param, well_model, terminal_output)
    {}

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    * \returns Report for simulator performance
    *
    * \note Strategies of coupling Flow and TPSA currently implemented:
    * \li fixed-stress: fixed-stress algorithm, i.e. iteratively solving Flow and TPSA equations in sequence
    * \li lagged:       one-way coupling where Flow is solved with TPSA info from previous time step
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIteration(const int iteration,
                                             const SimulatorTimerInterface& timer,
                                             NonlinearSolverType& nonlinear_solver)
    {
        SimulatorReportSingle report {};
        const auto& problem = this->simulator_.problem();
        if (problem.fixedStressScheme()) {
            report = nonlinearIterationFixedStressTPSA(iteration, timer, nonlinear_solver);
        }
        else if (problem.laggedScheme()) {
            report = nonlinearIterationLaggedTPSA(iteration, timer, nonlinear_solver);
        }
        else {
            std::string msg = "Unknown Flow-TPSA coupling scheme!";
            OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        return report;
    }

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics in a fixed-stress, iterative loop
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    * \returns Report for simulator performance
    *
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationFixedStressTPSA(const int iteration,
                                                            const SimulatorTimerInterface& timer,
                                                            NonlinearSolverType& nonlinear_solver)
    {
        // Runtime parameters
        const auto& [minSeqIter, maxSeqIter] = this->simulator_.problem().fixedStressParameters();
        SimulatorReportSingle reportFlow;

        // Max. no. of fixed-stress iterations reached: warn and move on
        if (seqIter_ >= maxSeqIter) {
            // Warning
            std::string msg = fmt::format("TPSA: Fixed-stress scheme reached max iterations (={})!", maxSeqIter);
            OpmLog::warning(msg);

            // Return true Flow convergence to move to next time step and reset other variables
            reportFlow.converged = true;
            seqIter_ = 0;

            return reportFlow;
        }

        // Prepare before first iteration
        if (seqIter_ == 0) {
            this->simulator_.problem().geoMechModel().prepareTPSA();
        }

        // Run Flow nonlinear iteration
        reportFlow = ParentType::nonlinearIteration(iteration, timer, nonlinear_solver);

        // Solve TPSA equations if:
        // (i)  Flow has converged and run more than min number of Newton iterations
        // (ii) we have run at least min. number of fixed-stress iterations
        if (reportFlow.converged && iteration >= this->param_.newton_min_iter_) {
            // Solve TPSA equations
            bool tpsaConv = solveTpsaEquations();
            ++seqIter_;

            // Fixed-stress convergence check:
            // If the initial residual error, hence check for no. linearizations == 1, was small enough, we have
            // convergence in the fixed-stress iterations
            if (tpsaConv
                && this->simulator_.problem().geoMechModel().newtonMethod().numLinearizations() == 1
                && seqIter_ >= minSeqIter) {
                // Info
                std::string msg = fmt::format("TPSA: Fixed-stress scheme converged in {} iterations", seqIter_);
                OpmLog::info(msg);

                // Reset
                seqIter_ = 0;

                return reportFlow;
            }
            // Throw error if TPSA did not converge. Will force time step cuts in the outer Flow loop.
            else if (!tpsaConv) {
                // Reset
                seqIter_ = 0;

                throw std::runtime_error("TPSA: Fixed stress scheme update failed!");
            }

            // Return Flow convergence false to do another fixed-stress iteration
            reportFlow.converged = false;
        }

        return reportFlow;
    }

    /*!
    * \brief Perform a nonlinear iteration updating Flow and TPSA geomechanics in a lagged scheme
    *
    * \param iteration Flow nonlinear iteration
    * \param timer Simulation timer
    * \param nonlinear_solver Nonlinear solver type
    * \returns Report for simulator performance
    *
    */
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationLaggedTPSA(const int iteration,
                                                       const SimulatorTimerInterface& timer,
                                                       NonlinearSolverType& nonlinear_solver)
    {
        // Run Flow nonlinear iteration
        auto reportFlow = ParentType::nonlinearIteration(iteration, timer, nonlinear_solver);

        // Update TPSA geomechanics from successful Flow run
        if (reportFlow.converged && iteration >= this->param_.newton_min_iter_) {
            // Prepare before TPSA solve
            this->simulator_.problem().geoMechModel().prepareTPSA();

            // Solve TPSA equations
            bool tpsaConv = solveTpsaEquations();

            // Throw error if TPSA did not converge. Will force time step cuts in the outer Flow loop.
            if (!tpsaConv) {
                throw std::runtime_error("TPSA: Lagged scheme update failed!");
            }
        }

        return reportFlow;
    }

    /*!
    * \brief Solve TPSA geomechanics equations
    *
    * \returns Bool indicating TPSA convergence
    *
    * \note Calls Newton method for TPSA
    */
    bool solveTpsaEquations()
    {
        // Run Newthon method for TPSA equations
        return this->simulator_.problem().geoMechModel().newtonMethod().apply();
    }

private:
    int seqIter_{0};
};  // class BlackoilModelTPSA

}  // namespace Opm

#endif
