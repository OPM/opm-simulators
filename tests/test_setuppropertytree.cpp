/*
  Copyright 2025 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#define BOOST_TEST_MODULE TestSetupPropertyTree

#ifndef HAVE_MPI
#define HAVE_MPI 0
#endif

#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <stdexcept>
#include <string>

BOOST_AUTO_TEST_SUITE(SystemCPR)

// Unit test 1: JSON file with system_cpr but missing well_solver must throw.
//
// options_system_cpr_missing_well.json contains reservoir_smoother and
// reservoir_solver but no well_solver. validateSystemCPRTree must catch this
// at setup time so the error is clear rather than a cryptic PropertyTree
// exception from deep inside SystemPreconditioner::initSubSolvers.
BOOST_AUTO_TEST_CASE(JSONMissingWellSolver)
{
    Opm::PropertyTree prm("options_system_cpr_missing_well.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm), std::invalid_argument);
}

// Unit test 2: JSON file with system_cpr but missing reservoir_solver must throw.
BOOST_AUTO_TEST_CASE(JSONMissingReservoirSolver)
{
    Opm::PropertyTree prm("options_system_cpr_missing_ressolver.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm), std::invalid_argument);
}

// Unit test 3: JSON file with system_cpr but missing reservoir_smoother must throw.
BOOST_AUTO_TEST_CASE(JSONMissingReservoirSmoother)
{
    Opm::PropertyTree prm("options_system_cpr_missing_smoother.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm), std::invalid_argument);
}

// Unit test 3: JSON has reservoir_solver but no preconditioner.type (or no preconditioner at all)
BOOST_AUTO_TEST_CASE(JSONReservoirSolverPreconditionerTypeMissingOrWrong)
{
    // Test missing type
    Opm::PropertyTree prm("options_system_cpr_missing_precond_type.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm), std::invalid_argument);

    // Test wrong type (e.g., "ilu")
    Opm::PropertyTree prm2("options_system_cpr_res_precond_not_cpr.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm2), std::invalid_argument);
}

// Regression test: --matrix-add-well-contributions=true is incompatible
// with system_cpr and must throw in both serial and parallel (MPI) runs.
//
// system_cpr handles reservoir-well coupling in the outer three-stage
// SystemPreconditioner. Enabling matrix well contributions merges well rows
// into the reservoir block before that preconditioner runs, which would
// double-count the coupling and corrupt the linear system.
// ISTLSolverRuntimeOptionProxy::createSolver delegates to this function,
// making the invariant unit-testable without the full simulator framework.
BOOST_AUTO_TEST_CASE(MatrixAddWellContributionsIncompatible)
{
    BOOST_CHECK_THROW(Opm::checkSystemCPRMatrixAddWell(true),  std::invalid_argument);
    BOOST_CHECK_NO_THROW(Opm::checkSystemCPRMatrixAddWell(false));
}

// With add_wells the pressure stage is assembled and solved by the system
// preconditioner itself, taking its solver settings from the coarsesolver
// sub-tree. Without that sub-tree there is nothing to solve the CPRW pressure
// system with, so it must be rejected at setup time rather than falling over
// inside SystemCprwPressureStage::buildStructure.
BOOST_AUTO_TEST_CASE(JSONAddWellsRequiresCoarseSolver)
{
    Opm::PropertyTree prm("options_system_cprw_missing_coarsesolver.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(prm), std::invalid_argument);

    Opm::PropertyTree complete("options_system_cprw_complete.json");
    BOOST_CHECK_NO_THROW(Opm::validateSystemCPRTree(complete));
}

// An approximate (Krylov) well solver stops on a tolerance and therefore does
// a different number of inner iterations per right-hand side, so the system
// preconditioner is no longer a fixed operator. Only a flexible outer solver
// may be combined with it; bicgstab must be rejected.
BOOST_AUTO_TEST_CASE(ApproximateWellSolverRequiresFlexibleOuterSolver)
{
    Opm::PropertyTree ok("options_system_cprw_approx_wells.json");
    BOOST_CHECK_NO_THROW(Opm::validateSystemCPRTree(ok));

    Opm::PropertyTree bad("options_system_cprw_approx_wells_bad_outer.json");
    BOOST_CHECK_THROW(Opm::validateSystemCPRTree(bad), std::invalid_argument);

    // The stationary well solvers stay valid with the default outer solver.
    Opm::PropertyTree exact("options_system_cprw_complete.json");
    BOOST_CHECK_EQUAL(exact.get<std::string>("preconditioner.well_solver.solver"), "umfpack");
    BOOST_CHECK_NO_THROW(Opm::validateSystemCPRTree(exact));
}

// system_cprw must produce the same tree as system_cpr apart from add_wells.
BOOST_AUTO_TEST_CASE(SystemCPRWEnablesAddWells)
{
    const Opm::FlowLinearSolverParameters p;
    const auto cpr = Opm::setupSystemCPR("system_cpr", p);
    const auto cprw = Opm::setupSystemCPR("system_cprw", p);

    const std::string key = "preconditioner.reservoir_solver.preconditioner.add_wells";
    BOOST_CHECK_EQUAL(cpr.get<bool>(key), false);
    BOOST_CHECK_EQUAL(cprw.get<bool>(key), true);

    // Same coarse solver in both, since that is what solves the pressure
    // system in either case.
    const std::string coarse
        = "preconditioner.reservoir_solver.preconditioner.coarsesolver.preconditioner.type";
    BOOST_CHECK_EQUAL(cpr.get<std::string>(coarse), cprw.get<std::string>(coarse));
}

BOOST_AUTO_TEST_SUITE_END() // SystemCPR
