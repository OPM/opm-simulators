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
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>


namespace Opm {

/*!
* \brief Internalize runtime parameters
*
* \note Overloading FlowLinearSolverParameters::init() to read TPSA specific runtime/default parameters
*/
void TpsaLinearSolverParameters::init()
{
    // Runtime parameters
    linear_solver_reduction_ = Parameters::Get<Parameters::TpsaLinearSolverReduction>();
    relaxed_linear_solver_reduction_ = Parameters::Get<Parameters::TpsaRelaxedLinearSolverReduction>();
    linear_solver_maxiter_ = Parameters::Get<Parameters::TpsaLinearSolverMaxIter>();
    linear_solver_restart_ = Parameters::Get<Parameters::TpsaLinearSolverRestart>();
    linear_solver_verbosity_ = Parameters::Get<Parameters::TpsaLinearSolverVerbosity>();
    ilu_relaxation_ = Parameters::Get<Parameters::TpsaIluRelaxation>();
    ilu_fillin_level_ = Parameters::Get<Parameters::TpsaIluFillinLevel>();
    newton_use_gmres_ = Parameters::Get<Parameters::TpsaUseGmres>();
    ignoreConvergenceFailure_ = Parameters::Get<Parameters::TpsaLinearSolverIgnoreConvergenceFailure>();
    linsolver_ = Parameters::Get<Parameters::TpsaLinearSolver>();
    linear_solver_print_json_definition_ = Parameters::Get<Parameters::TpsaLinearSolverPrintJsonDefinition>();

    // Hardcode use of CPU linear solvers (?)
    linear_solver_accelerator_ = Parameters::LinearSolverAcceleratorType::CPU;
}

/*!
* \brief Register TPSA linear solver runtime parameters
*/
void TpsaLinearSolverParameters::registerParameters()
{
    Parameters::Register<Parameters::TpsaLinearSolverReduction>
        ("Minimum residual reduction in TPSA linear solver for convergenc");
    Parameters::Register<Parameters::TpsaRelaxedLinearSolverReduction>
        ("A relaxed version of --tpsa-linear-solver-reduction (use with care!)");
    Parameters::Register<Parameters::TpsaLinearSolverMaxIter>
        ("Maximum TPSA linear iterations");
    Parameters::Register<Parameters::TpsaLinearSolverRestart>
        ("Number of iterations before restarting GMRES if --tpsa-use-gmres=true");
    Parameters::Register<Parameters::TpsaLinearSolverVerbosity>
        ("Level of verbosity in TPSA linear solver: 0 = off, 2 = all");
    Parameters::Register<Parameters::TpsaIluRelaxation>
        ("Relaxation factor for TPSA linear solver ILU preconditioner");
    Parameters::Register<Parameters::TpsaIluFillinLevel>
        ("Fill-in level of TPSA linear solver ILU preconditioner");
    Parameters::Register<Parameters::TpsaUseGmres>
        ("Use GMRES linear solver. If false, BiCGStab is used.");
    Parameters::Register<Parameters::TpsaLinearSolverIgnoreConvergenceFailure>
        ("Continue simulation even if TPSA linear solver did not converge");
    Parameters::Register<Parameters::TpsaLinearSolver>
        ("Configuration for linear solver. Valid preset options are: ilu0, dilu, amg or umfpack. "
         "Alternatively, you can request a configuration to be read from a JSON file by giving the filename here, "
         "ending with '.json.'");
    Parameters::Register<Parameters::TpsaLinearSolverPrintJsonDefinition>
        ("Print JSON formatted configuration of the TPSA linear solver. Can be used to make configuration JSON file "
        "for --tpsa-linear-solver");
}

/*!
* \brief Reset TPSA linear solver parameters
*
* \warning Should be the same as the corresponding Parameters defaults!
*/
void TpsaLinearSolverParameters::reset()
{
    linear_solver_reduction_ = 1e-3;
    relaxed_linear_solver_reduction_ = 1e-3;
    linear_solver_maxiter_ = 200;
    linear_solver_restart_ = 40;
    linear_solver_verbosity_ = 0;
    ilu_relaxation_ = 9;
    ilu_fillin_level_ = 0;
    newton_use_gmres_ = false;
    ignoreConvergenceFailure_ = false;
    linsolver_ = "ilu0";
    linear_solver_print_json_definition_ = false;
}

}  // namespace Opm
