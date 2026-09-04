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
#ifndef OPM_SYSTEMPRECONDITIONER_HEADER_INCLUDED
#define OPM_SYSTEMPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/system/MultiComm.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

namespace Opm {

// Preconditioner for the coupled reservoir-well system.
//
// Templated on scalar type, reservoir operator and communication types to
// unify sequential and parallel implementations. The 3-stage algorithm:
//   1. Reservoir CPR solve
//   2. Well solve + reservoir smoothing
//   3. Final well solve
//
// For parallel runs, copyOwnerToAll synchronises overlap DOFs before
// each reservoir sub-solve.
template <class Scalar, class ResOp, class ResComm = Dune::Amg::SequentialInformation>
class SystemPreconditioner : public Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOp>;
    using WellOperator = Dune::MatrixAdapter<WWMatrix<Scalar>, WellVector<Scalar>, WellVector<Scalar>>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    // Sequential constructor (enabled only for non-parallel specializations).
    SystemPreconditioner(const SystemMatrix<Scalar>& S,
                         const std::function<ResVector<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm) requires (!isParallel);

    // Parallel constructor (enabled only for parallel specializations).
    SystemPreconditioner(const SystemMatrix<Scalar>& S,
                         const std::function<ResVector<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm,
                         const ResComm& resComm) requires (isParallel);

    void pre(SystemVector<Scalar>&, SystemVector<Scalar>&) override
    {
    }

    void post(SystemVector<Scalar>&) override
    {
    }

    Dune::SolverCategory::Category category() const override
    {
        if constexpr (isParallel)
            return Dune::SolverCategory::overlapping;
        else
            return Dune::SolverCategory::sequential;
    }

    void update() override;

    void updateForChangedWellStructure();

    bool hasPerfectUpdate() const override
    {
        return true;
    }

//   System matrix block structure:
//
//       [ A  C ] [ x_res ]   [ resRes ]
//   S = [ B  D ] [ x_well ] = [ wRes  ]
//
//   A = reservoir-reservoir (top-left)
//   C = reservoir-well coupling (top-right)
//   B = well-reservoir coupling (bottom-left)
//   D = well-well (bottom-right)
    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d) override;

private:
    const SystemMatrix<Scalar>& S_;
    const ResComm* resComm_ = nullptr;
    int pressureIndex_ = 0;
    static constexpr int dummyWellPressureIndex = std::numeric_limits<int>::min();
    Opm::PropertyTree wellprm_;

    std::unique_ptr<ResOp> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<ResFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    WellVector<Scalar> wSol_;
    ResVector<Scalar> resSol_;
    ResVector<Scalar> dresSol_;
    WellVector<Scalar> dwSol_;
    ResVector<Scalar> tmp_resRes_;
    WellVector<Scalar> tmp_wRes_;
    ResVector<Scalar> resRes_;
    WellVector<Scalar> wRes_;

    void syncResVector(ResVector<Scalar>& v);

    void initWellSolver();

    void initSubSolvers(const Opm::PropertyTree& prm,
                        const std::function<ResVector<Scalar>()>& weightsCalculator);

    void initWorkVectors();

    void resizeReservoirWorkVectors();

    void resizeWellWorkVectors();
};

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONER_HEADER_INCLUDED
