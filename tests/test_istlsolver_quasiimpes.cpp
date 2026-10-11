/*
  Copyright 2026 Equinor ASA

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

#include "config.h"

#define BOOST_TEST_MODULE ISTLSolverQuasiImpes

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <array>
#include <cstddef>
#include <string>

namespace
{

using TypeTag = Opm::Properties::TTag::TestTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

class TestISTLSolver : public Opm::ISTLSolver<TypeTag>
{
    using Parent = Opm::ISTLSolver<TypeTag>;

public:
    using Matrix = Parent::Matrix;
    using Vector = Parent::Vector;
    using Preconditioner = Parent::AbstractPreconditionerType;

    static constexpr auto pressureIndex = Parent::pressureIndex;

    TestISTLSolver(const Simulator& simulator,
                   const std::string& preconditionerType,
                   const bool threadedWeights)
        : Parent(simulator, solverParameters(threadedWeights), /*forceSerial=*/true)
    {
        auto& prm = this->prm_.front();
        prm.put("preconditioner.type", preconditionerType);
        prm.put("preconditioner.weight_type", std::string {"quasiimpes"});

        // Apply only the CPR pressure stage, with an exact coarse solve, so that
        // the correction has a simple analytical value and a fine smoother
        // cannot mask stale weights.
        prm.put("preconditioner.pre_smooth", 0);
        prm.put("preconditioner.post_smooth", 0);
        prm.put("preconditioner.coarsesolver.solver", std::string {"loopsolver"});
        prm.put("preconditioner.coarsesolver.maxiter", 1);
        prm.put("preconditioner.coarsesolver.tol", 1.0e-12);
        prm.put("preconditioner.coarsesolver.preconditioner.type", std::string {"paroverilu0"});
        prm.put("preconditioner.coarsesolver.preconditioner.relaxation", 1.0);
    }

    Preconditioner* preconditioner()
    {
        return this->flexibleSolver_.front().pre_;
    }

    bool rebuildsOnNextPrepare() const
    {
        return this->shouldCreateSolver();
    }

    Vector applyPreconditioner(Vector rhs)
    {
        Vector correction(rhs.size());
        correction = 0.0;
        auto& pre = *preconditioner();
        pre.pre(correction, rhs);
        pre.apply(correction, rhs);
        pre.post(correction);
        return correction;
    }

private:
    static Opm::FlowLinearSolverParameters solverParameters(const bool threadedWeights)
    {
        Opm::FlowLinearSolverParameters parameters;
        parameters.init(false);
        parameters.linsolver_ = "cpr_quasiimpes";
        parameters.cpr_reuse_setup_ = 3;
        parameters.cpr_weights_thread_parallel_ = threadedWeights;
        parameters.linear_solver_print_json_definition_ = false;
        return parameters;
    }
};

TestISTLSolver::Matrix
diagonalMatrix(const std::size_t size)
{
    using Matrix = TestISTLSolver::Matrix;
    Matrix matrix(size, size, size, Matrix::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        row.insert(row.index());
    }
    matrix = 0.0;
    for (std::size_t row = 0; row < size; ++row) {
        auto& diagonal = matrix[row][row];
        for (int component = 0; component < diagonal.rows; ++component) {
            diagonal[component][component] = 1.0;
        }
    }
    return matrix;
}

void
checkCorrection(const TestISTLSolver::Vector& correction,
                const double pressureCorrection,
                const std::size_t otherIndex,
                const double otherCorrection)
{
    for (std::size_t row = 0; row < correction.size(); ++row) {
        for (std::size_t component = 0; component < correction[row].size(); ++component) {
            const auto expected = component == TestISTLSolver::pressureIndex ? pressureCorrection
                : component == otherIndex                                    ? otherCorrection
                                                                             : 0.0;
            BOOST_CHECK_SMALL(correction[row][component] - expected, 1.0e-12);
        }
    }
}

} // namespace

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(QuasiImpesWeightsFollowMatrixValuesWhenPreconditionerIsReused)
{
    const auto simulator = Opm::initSimulator<TypeTag>("drsdt_initial_rs.DATA",
                                                       "test_istlsolver_quasiimpes",
                                                       /*threads_per_process=*/2);
    const auto size = simulator->model().numGridDof();
    BOOST_REQUIRE_GT(size, 0u);

    // The synthetic matrix contains the complete system. There is no
    // separate well operator in this test, including during initPrepare().
    Opm::Parameters::SetDefault<Opm::Parameters::MatrixAddWellContributions>(true);

    using VectorBlock = TestISTLSolver::Vector::block_type;
    static_assert(VectorBlock::size() > 1);
    constexpr auto pressureIndex = TestISTLSolver::pressureIndex;
    constexpr auto otherIndex = (pressureIndex + 1) % VectorBlock::size();

    for (const auto& preconditionerType : std::array {"cpr", "cprt"}) {
        const bool transpose = std::string {preconditionerType} == "cprt";
        for (const bool threadedWeights : {false, true}) {
            BOOST_TEST_CONTEXT("preconditioner=" << preconditionerType
                                                 << ", threaded weights=" << threadedWeights)
            {
                auto matrix = diagonalMatrix(size);
                TestISTLSolver::Vector rhs(size);
                rhs = 0.0;
                for (auto& block : rhs) {
                    block[transpose ? pressureIndex : otherIndex] = 1.0;
                }

                TestISTLSolver solver(*simulator, preconditionerType, threadedWeights);
                solver.prepare(matrix, rhs);
                BOOST_REQUIRE(solver.preconditioner() != nullptr);
                checkCorrection(
                    solver.applyPreconditioner(rhs), transpose ? 1.0 : 0.0, otherIndex, 0.0);

                // Change values in the same matrix object. In (other, pressure)
                // order its diagonal sub-block is now [[2, 1], [1, 3]], whose
                // inverse is [[3, -1], [-1, 2]] / 5. CPR therefore gives -1/5
                // in pressure for e_other; CPRT gives (-1/5, 2/5) for e_pressure.
                for (std::size_t row = 0; row < size; ++row) {
                    auto& diagonal = matrix[row][row];
                    diagonal[otherIndex][otherIndex] = 2.0;
                    diagonal[otherIndex][pressureIndex] = 1.0;
                    diagonal[pressureIndex][otherIndex] = 1.0;
                    diagonal[pressureIndex][pressureIndex] = 3.0;
                }

                // A rebuild would recompute the weights from the new values in
                // any case, so the second prepare() has to update in place.
                BOOST_REQUIRE(!solver.rebuildsOnNextPrepare());
                solver.prepare(matrix, rhs);
                checkCorrection(solver.applyPreconditioner(rhs),
                                transpose ? 0.4 : -0.2,
                                otherIndex,
                                transpose ? -0.2 : 0.0);
            }
        }
    }
}
