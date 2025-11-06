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
#ifndef TPSA_LINEAR_SOLVER_PARAMETERS_HPP
#define TPSA_LINEAR_SOLVER_PARAMETERS_HPP

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include <string>


// Default runtime parameters
namespace Opm::Parameters {

struct TpsaLinearSolverReduction { static constexpr double value = 1e-3; };
struct TpsaRelaxedLinearSolverReduction { static constexpr double value = 1e-3; };
struct TpsaLinearSolverMaxIter { static constexpr int value = 200; };
struct TpsaLinearSolverRestart { static constexpr int value = 40; };
struct TpsaLinearSolverVerbosity { static constexpr int value = 0; };
struct TpsaIluRelaxation { static constexpr double value = 0.9; };
struct TpsaIluFillinLevel { static constexpr int value = 0; };
struct TpsaUseGmres { static constexpr bool value = false; };
struct TpsaLinearSolverIgnoreConvergenceFailure { static constexpr bool value = false; };
struct TpsaLinearSolver { static constexpr auto value = "ilu0"; };
struct TpsaLinearSolverPrintJsonDefinition { static constexpr auto value = false; };

}  // namespace Opm::Parameters

namespace Opm {

struct TpsaLinearSolverParameters : public FlowLinearSolverParameters
{
    void init();
    static void registerParameters();
    void reset();
};

}  // namespace Opm

#endif