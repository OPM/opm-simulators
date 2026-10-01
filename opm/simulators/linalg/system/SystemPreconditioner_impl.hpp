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
#ifndef OPM_SYSTEMPRECONDITIONER_IMPL_HEADER_INCLUDED
#define OPM_SYSTEMPRECONDITIONER_IMPL_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>

namespace Opm {

template <class Scalar, class ResOp, class ResComm>
SystemPreconditioner<Scalar,ResOp,ResComm>::
SystemPreconditioner(const SystemMatrix<Scalar>& S,
                     const std::function<ResVector<Scalar>()>& weightsCalculator,
                     int pressureIndex,
                     const Opm::PropertyTree& prm) requires (!isParallel)
    : S_(S)
    , pressureIndex_(pressureIndex)
{
    initSubSolvers(prm, weightsCalculator);
    initWorkVectors();
}

template <class Scalar, class ResOp, class ResComm>
SystemPreconditioner<Scalar,ResOp,ResComm>::
SystemPreconditioner(const SystemMatrix<Scalar>& S,
                     const std::function<ResVector<Scalar>()>& weightsCalculator,
                     int pressureIndex,
                     const Opm::PropertyTree& prm,
                     const ResComm& resComm) requires (isParallel)
    : S_(S)
    , resComm_(&resComm)
    , pressureIndex_(pressureIndex)
{
    initSubSolvers(prm, weightsCalculator);
    initWorkVectors();
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
update()
{
    resSolver_->preconditioner().update();
    resSmoother_->preconditioner().update();
    wellSolver_->preconditioner().update();
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
updateForChangedWellStructure()
{
    resSolver_->preconditioner().update();
    resSmoother_->preconditioner().update();
    initWellSolver();
    resizeWellWorkVectors();
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d)
{
    // Extract blocks using the agreed convention
    const auto& A = S_[_0][_0];
    const auto& C = S_[_0][_1];
    const auto& B = S_[_1][_0];
    const auto& D = S_[_1][_1];

    resRes_ = d[_0];
    wRes_ = d[_1];
    resSol_ = 0.0;
    wSol_ = 0.0;

    // Stage 1: Reservoir CPR solve
    {
        Dune::InverseOperatorResult res_result;
        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        syncResVector(tmp_resRes_);
        resSolver_->apply(dresSol_, tmp_resRes_, res_result);
        resSol_ += dresSol_;
        // resRes_ -= A * dresSol_
        A.mmv(dresSol_, resRes_);
        // wRes_ -= B * dresSol_
        B.mmv(dresSol_, wRes_);
    }

    // Stage 2: Well solve + reservoir system smoothing
    {
        Dune::InverseOperatorResult well_result;
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_result);
        wSol_ += dwSol_;
        // resRes_ -= C * dwSol_
        C.mmv(dwSol_, resRes_);
        // resRes_ -= D * dwSol_
        D.mmv(dwSol_, wRes_);

        Dune::InverseOperatorResult res_result;
        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        syncResVector(tmp_resRes_);
        resSmoother_->apply(dresSol_, tmp_resRes_, res_result);
        resSol_ += dresSol_;
        // wRes_ -= B * dresSol_
        B.mmv(dresSol_, wRes_);
    }

    // Stage 3: Final well solve
    {
        Dune::InverseOperatorResult well_result;
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_result);
        wSol_ += dwSol_;
    }

    syncResVector(resSol_);
    v[_0] = resSol_;
    v[_1] = wSol_;
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
syncResVector(ResVector<Scalar>& v)
{
    if constexpr (isParallel) {
        resComm_->copyOwnerToAll(v, v);
    }
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
initWellSolver()
{
    wop_ = std::make_unique<WellOperator>(S_[_1][_1]);
    std::function<WellVector<Scalar>()> weightsCalculatorWell;
    wellSolver_ = std::make_unique<WellFlexibleSolverType>(
        *wop_, wellprm_, weightsCalculatorWell, dummyWellPressureIndex);
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
initSubSolvers(const PropertyTree& prm,
               const std::function<ResVector<Scalar>()>& weightsCalculator)
{
    auto resprm = prm.get_child("reservoir_solver");
    auto resprmsmoother = prm.get_child("reservoir_smoother");
    wellprm_ = prm.get_child("well_solver");

    if constexpr (isParallel) {
        rop_ = std::make_unique<ResOp>(S_[_0][_0], *resComm_);
        resSolver_ = std::make_unique<ResFlexibleSolverType>(
            *rop_, *resComm_, resprm, weightsCalculator, pressureIndex_);
        resSmoother_ = std::make_unique<ResFlexibleSolverType>(
            *rop_, *resComm_, resprmsmoother, weightsCalculator, pressureIndex_);
    } else {
        rop_ = std::make_unique<ResOp>(S_[_0][_0]);
        resSolver_ = std::make_unique<ResFlexibleSolverType>(
            *rop_, resprm, weightsCalculator, pressureIndex_);
        resSmoother_ = std::make_unique<ResFlexibleSolverType>(
            *rop_, resprmsmoother, weightsCalculator, pressureIndex_);
    }

    initWellSolver();
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
initWorkVectors()
{
    resizeReservoirWorkVectors();
    resizeWellWorkVectors();
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
resizeReservoirWorkVectors()
{
    const auto numRes = S_[_0][_0].N();
    resSol_.resize(numRes);
    dresSol_.resize(numRes);
    tmp_resRes_.resize(numRes);
    resRes_.resize(numRes);
}

template <class Scalar, class ResOp, class ResComm>
void
SystemPreconditioner<Scalar,ResOp,ResComm>::
resizeWellWorkVectors()
{
    const auto numWell = S_[_1][_1].N();
    wSol_.resize(numWell);
    dwSol_.resize(numWell);
    tmp_wRes_.resize(numWell);
    wRes_.resize(numWell);
}

} // namespace Opm

#endif
