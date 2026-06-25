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

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <stdexcept>

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

BOOST_AUTO_TEST_SUITE_END() // SystemCPR
