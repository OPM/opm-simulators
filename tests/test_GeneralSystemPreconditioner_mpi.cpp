/*
  Copyright Equinor ASA 2026

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
#define BOOST_TEST_MODULE GeneralSystemPreconditionerMpi
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/system/GeneralSystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <dune/common/parallel/mpihelper.hh>

#include <string>

BOOST_AUTO_TEST_CASE(ParallelSweepMatchesFixedAndKeepsCopiesConsistent)
{
    const auto& world = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_EQUAL(world.size(), 2);
    using Comm = Opm::ParResComm;
    using LocalIndex = Comm::ParallelIndexSet::LocalIndex;
    using Attributes = Dune::OwnerOverlapCopyAttributeSet;
    using General = Opm::GeneralSystemPreconditioner<double, Opm::ParResOperator<double>, Comm>;
    using Fixed = Opm::SystemPreconditioner<double, Opm::ParResOperator<double>, Comm>;
    constexpr auto res = Dune::Indices::_0;
    constexpr auto well = Dune::Indices::_1;
    const int ownerCell = world.rank();
    const int copyCell = 1 - ownerCell;

    // Each rank owns one cell and has a copy of its neighbour. Its one well
    // is local and perforates only the owned cell.
    Comm comm(world);
    auto& indices = comm.indexSet();
    indices.beginResize();
    for (int cell = 0; cell < 2; ++cell)
        indices.add(cell, LocalIndex(cell, cell == ownerCell ? Attributes::owner : Attributes::copy, true));
    indices.endResize();
    comm.remoteIndices().template rebuild<false>();

    Opm::RRMatrix<double> A(2, 2, Opm::RRMatrix<double>::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(0);
        row.insert(1);
    }
    A = 0.0;
    for (int eq = 0; eq < Opm::numResDofs; ++eq) {
        A[ownerCell][ownerCell][eq][eq] = 4.0 + ownerCell + eq;
        A[ownerCell][copyCell][eq][eq] = -1.0;
        A[copyCell][copyCell][eq][eq] = 1.0;
    }
    Opm::RWMatrix<double> C(2, 1, Opm::RWMatrix<double>::row_wise);
    for (auto row = C.createbegin(); row != C.createend(); ++row)
        if (row.index() == static_cast<std::size_t>(ownerCell)) row.insert(0);
    C = 0.0;
    C[ownerCell][0][0][3] = 0.2;
    Opm::WRMatrix<double> B(1, 2, Opm::WRMatrix<double>::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) row.insert(ownerCell);
    B = 0.0;
    B[0][ownerCell][3][0] = 0.3;
    Opm::WWMatrix<double> D(1, 1, Opm::WWMatrix<double>::row_wise);
    for (auto row = D.createbegin(); row != D.createend(); ++row) row.insert(0);
    D = 0.0;
    for (int eq = 0; eq < Opm::numWellDofs; ++eq) D[0][0][eq][eq] = 2.0 + eq;

    Opm::WellDofLayout layout;
    layout.wellBlockOffsets = {0, 1};
    Opm::SystemMatrix<double> system;
    system.A = &A;
    system.B = &B;
    system.C = &C;
    system.D = &D;
    system.wellLayout = &layout;
    Opm::SystemVector<double> weights, rhs;
    weights[res].resize(2);
    weights[res] = 0.0;
    weights[res][0][0] = weights[res][1][0] = 1.0;
    weights[well].resize(1);
    weights[well] = 0.0;
    weights[well][0][3] = 1.0;
    const auto calculateWeights = [weights]() { return weights; };
    rhs[res].resize(2);
    rhs[well].resize(1);
    rhs[res][ownerCell] = 1.0 + ownerCell;
    rhs[res][copyCell] = -1000.0 - ownerCell;
    rhs[well] = 0.5 + ownerCell;

    const Opm::FlowLinearSolverParameters parameters;
    for (const bool coarse : {false, true}) {
        auto fixedPrm = Opm::setupSystemCPR(coarse ? "system_cprw" : "system_cpr", parameters)
                            .get_child("preconditioner");
        auto generalPrm = Opm::setupGeneralSystemCPR(
            coarse ? "general_system_cprw" : "general_system_cpr", parameters)
                              .get_child("preconditioner");
        // Use a single ILU application on the coarse system so this small
        // test exercises copy residuals without AMG hierarchy heuristics.
        if (coarse) {
            fixedPrm.put("reservoir_solver.preconditioner.coarsesolver.preconditioner.type",
                         std::string{"ilu0"});
            generalPrm.put("coarsesolver.preconditioner.type", std::string{"ilu0"});
        }
        Fixed fixed(system, calculateWeights, 0, fixedPrm, comm);
        General general(system, calculateWeights, 0, generalPrm, comm);
        auto actual = rhs, expected = rhs;
        actual = 0.0;
        expected = 0.0;
        fixed.apply(expected, rhs);
        general.apply(actual, rhs);
        auto consistent = actual[res];
        comm.copyOwnerToAll(consistent, consistent);
        BOOST_TEST_CONTEXT("coarse=" << coarse << " rank=" << world.rank()) {
            // The owner correction must agree with the fixed implementation,
            // and the returned copies must agree with their owners.
            for (int cell = 0; cell < 2; ++cell)
                for (int eq = 0; eq < Opm::numResDofs; ++eq) {
                    BOOST_CHECK_SMALL(actual[res][cell][eq] - expected[res][cell][eq], 1e-12);
                    BOOST_CHECK_EQUAL(actual[res][cell][eq], consistent[cell][eq]);
                }
            for (int eq = 0; eq < Opm::numWellDofs; ++eq)
                BOOST_CHECK_SMALL(actual[well][0][eq] - expected[well][0][eq], 1e-12);
            BOOST_CHECK_GT(actual.two_norm(), 0.0);
        }
    }
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    [[maybe_unused]] const auto& helper = Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
