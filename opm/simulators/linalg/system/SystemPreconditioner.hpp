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

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/system/CpuSystemBackend.hpp>
#include <opm/simulators/linalg/system/MultiComm.hpp>
#include <opm/simulators/linalg/system/SystemCprwPressureStage.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>

namespace Opm {

// Reservoir operator/comm types used as template arguments.
template<typename Scalar>
using SeqResOperator = Dune::MatrixAdapter<RRMatrix<Scalar>, ResVector<Scalar>, ResVector<Scalar>>;

#if HAVE_MPI
using ParResComm = Dune::OwnerOverlapCopyCommunication<int, int>;
template<typename Scalar>
using ParResOperator = Dune::OverlappingSchwarzOperator<RRMatrix<Scalar>, ResVector<Scalar>, ResVector<Scalar>, ParResComm>;
#endif

// Preconditioner for the coupled reservoir-well system.
//
// The Backend traits struct supplies the matrix, vector and sub-solver types,
// so that sequential, parallel and GPU runs share one implementation.  It
// defaults to CpuSystemBackend<Scalar, ResOp, ResComm>; the GPU solver uses
// gpuistl::GpuSystemBackend (see gpuistl::GpuSystemPreconditioner).
// The 3-stage algorithm:
//   1. Reservoir CPR solve
//   2. Well solve + reservoir smoothing
//   3. Final well solve
//
// For parallel runs, copyOwnerToAll synchronises overlap DOFs before
// each reservoir sub-solve.
//
// On the GPU backend the reservoir sub-solvers run on the device, while the
// well sub-solver stays on the CPU: the well system is small, so the cost of
// moving its right-hand side and solution across is negligible. The CPRW
// pressure stage (add_wells) is CPU only.
template <class Scalar,
          class ResOp,
          class ResComm = Dune::Amg::SequentialInformation,
          class Backend = CpuSystemBackend<Scalar, ResOp, ResComm>>
class SystemPreconditioner : public Dune::PreconditionerWithUpdate<typename Backend::SysVector,
                                                                   typename Backend::SysVector>
{
    static_assert(std::is_same_v<typename Backend::Scalar, Scalar>);
    static_assert(std::is_same_v<typename Backend::ResComm, ResComm>);

public:
    using SysMatrix = typename Backend::SysMatrix;
    using SysVector = typename Backend::SysVector;
    using ResWork = typename Backend::ResWork;
    using WellWork = typename Backend::WellWork;

    static constexpr bool isGpu = Backend::isGpu;
    static constexpr bool isParallel = Backend::isParallel;

    using ResFlexibleSolverType = typename Backend::ResSolver;
    using WellOperator = Dune::MatrixAdapter<WWMatrix<Scalar>, WellVector<Scalar>, WellVector<Scalar>>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    // Sequential constructor (enabled only for non-parallel specializations).
    SystemPreconditioner(const SysMatrix& S,
                         const std::function<SysVector()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm)
        requires(!isParallel)
        : S_(S)
        , pressureIndex_(pressureIndex)
    {
        initSubSolvers(prm, weightsCalculator);
        initWorkVectors();
    }

    // Parallel constructor (enabled only for parallel specializations).
    SystemPreconditioner(const SysMatrix& S,
                         const std::function<SysVector()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm,
                         const ResComm& resComm)
        requires(isParallel)
        : S_(S)
        , resComm_(&resComm)
        , pressureIndex_(pressureIndex)
    {
        initSubSolvers(prm, weightsCalculator);
        initWorkVectors();
    }

    void pre(SysVector&, SysVector&) override
    {
    }

    void post(SysVector&) override
    {
    }

    Dune::SolverCategory::Category category() const override
    {
        if constexpr (isParallel)
            return Dune::SolverCategory::overlapping;
        else
            return Dune::SolverCategory::sequential;
    }

    void update() override
    {
        updatePressureStage();
        Backend::updateResSolver(*resSmoother_);
        wellSolver_->preconditioner().update();
    }

    void updateForChangedWellStructure()
    {
        if constexpr (!isGpu) {
            if (cprwStage_) {
                // The coarse system carries one unknown per well, so a changed
                // well structure changes its dimension and pattern.
                weights_ = weightsCalculator_();
                cprwStage_->buildStructure(weights_);
            } else {
                Backend::updateResSolver(*resSolver_);
            }
        } else {
            Backend::updateResSolver(*resSolver_);
        }
        Backend::updateResSolver(*resSmoother_);
        initWellSolver();
        resizeWellWorkVectors();
    }

    bool hasPerfectUpdate() const override
    {
        return Backend::hasPerfectUpdate;
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
    void apply(SysVector& v, const SysVector& d) override
    {
        // Extract blocks using the agreed convention
        const auto& A = *S_.A;
        const auto& C = *S_.C;
        const auto& B = *S_.B;
        const auto& D = *S_.D;

        resRes_ = Backend::res(d);
        wRes_ = Backend::well(d);
        Backend::setZero(resSol_);
        Backend::setZero(wSol_);

        // Stage 1: pressure solve.  Either reservoir-only CPR, or -- with
        // add_wells -- CPRW on the full system, in which case the stage also
        // produces a correction for the well unknowns.
        applyPressureStage(A, B, C, D);

        // Stage 2: Well solve + reservoir system smoothing
        {
            Backend::setZero(dwSol_);
            applyWellSolver(dwSol_, wRes_);
            Backend::add(wSol_, dwSol_);
            // resRes_ -= C * dwSol_
            Backend::mmv(C, dwSol_, resRes_);
            // resRes_ -= D * dwSol_
            Backend::mmv(D, dwSol_, wRes_);

            Dune::InverseOperatorResult res_result;
            Backend::setZero(dresSol_);
            tmp_resRes_ = resRes_;
            syncResVector(tmp_resRes_);
            resSmoother_->apply(dresSol_, tmp_resRes_, res_result);
            resSol_ += dresSol_;
            // wRes_ -= B * dresSol_
            Backend::mmv(B, dresSol_, wRes_);
        }

        // Stage 3: Final well solve
        {
            Backend::setZero(dwSol_);
            applyWellSolver(dwSol_, wRes_);
            Backend::add(wSol_, dwSol_);
        }

        syncResVector(resSol_);
        Backend::res(v) = resSol_;
        Backend::well(v) = wSol_;
    }

private:
    const SysMatrix& S_;
    const ResComm* resComm_ = nullptr;
    int pressureIndex_ = 0;
    static constexpr int dummyWellPressureIndex = std::numeric_limits<int>::min();
    Opm::PropertyTree wellprm_;

    // Only used on the CPU: the GPU reservoir sub-solvers build their own operator.
    std::unique_ptr<typename Backend::ResOperator> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<ResFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    // Non-null when the pressure stage includes the well unknowns (CPRW).
    // Then resSolver_ is not built and stage 1 goes through cprwStage_.
    // CPU only.
    using CprwStage = SystemCprwPressureStage<Scalar, ResComm>;
    std::unique_ptr<CprwStage> cprwStage_;
    std::function<SysVector()> weightsCalculator_;
    SysVector weights_;

    WellWork wSol_;
    ResWork resSol_;
    ResWork dresSol_;
    WellWork dwSol_;
    ResWork tmp_resRes_;
    WellWork tmp_wRes_;
    ResWork resRes_;
    WellWork wRes_;

    // GPU backend only: host-side copies of the well right-hand side and
    // correction, used to run the CPU well sub-solver.
    ConditionalStorage<isGpu, WellVector<Scalar>> hostWRes_;
    ConditionalStorage<isGpu, WellVector<Scalar>> hostDwSol_;

    template <class MatA, class MatB, class MatC, class MatD>
    void applyPressureStage(const MatA& A, const MatB& B, const MatC& C, const MatD& D)
    {
        if constexpr (!isGpu) {
            if (cprwStage_) {
                dresSol_ = 0.0;
                dwSol_ = 0.0;
                tmp_resRes_ = resRes_;
                syncResVector(tmp_resRes_);
                cprwStage_->apply(tmp_resRes_, wRes_, weights_, dresSol_, dwSol_);
                resSol_ += dresSol_;
                // resRes_ -= A * dresSol_
                A.mmv(dresSol_, resRes_);
                // wRes_ -= B * dresSol_
                B.mmv(dresSol_, wRes_);
                if (cprwStage_->prolongatesWellPressure()) {
                    wSol_ += dwSol_;
                    // resRes_ -= C * dwSol_ ;  wRes_ -= D * dwSol_
                    C.mmv(dwSol_, resRes_);
                    D.mmv(dwSol_, wRes_);
                }
                return;
            }
        }

        Dune::InverseOperatorResult res_result;
        Backend::setZero(dresSol_);
        tmp_resRes_ = resRes_;
        syncResVector(tmp_resRes_);
        resSolver_->apply(dresSol_, tmp_resRes_, res_result);
        resSol_ += dresSol_;
        // resRes_ -= A * dresSol_
        Backend::mmv(A, dresSol_, resRes_);
        // wRes_ -= B * dresSol_
        Backend::mmv(B, dresSol_, wRes_);
    }

    void updatePressureStage()
    {
        if constexpr (!isGpu) {
            if (cprwStage_) {
                weights_ = weightsCalculator_();
                cprwStage_->update(weights_);
                return;
            }
        }
        Backend::updateResSolver(*resSolver_);
    }

    // Solve with the well sub-solver.  The solver may modify its right-hand
    // side, so it always works on a copy.
    void applyWellSolver(WellWork& sol, const WellWork& rhs)
    {
        Dune::InverseOperatorResult well_result;
        if constexpr (isGpu) {
            if (hostWRes_->size() == 0) {
                return;
            }
            rhs.copyToHost(*hostWRes_);
            *hostDwSol_ = 0.0;
            wellSolver_->apply(*hostDwSol_, *hostWRes_, well_result);
            sol.copyFromHost(*hostDwSol_);
        } else {
            tmp_wRes_ = rhs;
            wellSolver_->apply(sol, tmp_wRes_, well_result);
        }
    }

    void syncResVector(ResWork& v)
    {
        if constexpr (isParallel) {
            resComm_->copyOwnerToAll(v, v);
        }
    }

    void initWellSolver()
    {
        wop_ = std::make_unique<WellOperator>(Backend::cpuWellMatrix(S_));
        std::function<WellVector<Scalar>()> weightsCalculatorWell;
        wellSolver_ = std::make_unique<WellFlexibleSolverType>(
            *wop_, wellprm_, weightsCalculatorWell, dummyWellPressureIndex);
    }

    void initSubSolvers(const Opm::PropertyTree& prm,
                        const std::function<SysVector()>& weightsCalculator)
    {
        auto resprm = prm.get_child("reservoir_solver");
        auto resprmsmoother = prm.get_child("reservoir_smoother");
        wellprm_ = prm.get_child("well_solver");

        // add_wells is the same switch the classic CPR/CPRW pair uses: it
        // promotes the pressure stage from reservoir-only CPR to CPRW over the
        // full (reservoir, well) system.
        const bool addWells = resprm.get("preconditioner.add_wells", false);

        // The weights arrive from the outer layer for the whole system; the
        // reservoir-only sub-solvers want just their own part of them.
        std::function<ResWork()> resWeightCalc;
        if (weightsCalculator) {
            resWeightCalc
                = [weightsCalculator]() { return ResWork(Backend::res(weightsCalculator())); };
        }

        if constexpr (isGpu) {
            if (addWells) {
                OPM_THROW(std::invalid_argument,
                          "The CPRW pressure stage (add_wells) of system_cpr is not "
                          "supported with --linear-solver-accelerator=gpu.");
            }
            // The GPU wrapper builds its own operator around the reservoir matrix.
            resSmoother_ = std::make_unique<ResFlexibleSolverType>(*S_.A,
                                                                   /*parallel=*/false,
                                                                   resprmsmoother,
                                                                   pressureIndex_,
                                                                   resWeightCalc,
                                                                   /*forceSerial=*/true,
                                                                   /*comm=*/nullptr);
            resSolver_ = std::make_unique<ResFlexibleSolverType>(*S_.A,
                                                                 /*parallel=*/false,
                                                                 resprm,
                                                                 pressureIndex_,
                                                                 resWeightCalc,
                                                                 /*forceSerial=*/true,
                                                                 /*comm=*/nullptr);
        } else {
            if constexpr (isParallel) {
                rop_ = std::make_unique<typename Backend::ResOperator>(*S_.A, *resComm_);
                resSmoother_ = std::make_unique<ResFlexibleSolverType>(
                    *rop_, *resComm_, resprmsmoother, resWeightCalc, pressureIndex_);
            } else {
                rop_ = std::make_unique<typename Backend::ResOperator>(*S_.A);
                resSmoother_ = std::make_unique<ResFlexibleSolverType>(
                    *rop_, resprmsmoother, resWeightCalc, pressureIndex_);
            }

            if (addWells) {
                if (!weightsCalculator) {
                    OPM_THROW(
                        std::invalid_argument,
                        "The CPRW pressure stage (add_wells) needs a weights calculator, but "
                        "none was configured. Set reservoir_solver.preconditioner.weight_type.");
                }
                weightsCalculator_ = weightsCalculator;
                weights_ = weightsCalculator_();

                auto coarseprm = resprm.get_child_optional("preconditioner.coarsesolver")
                    ? resprm.get_child("preconditioner.coarsesolver")
                    : PropertyTree();
                const auto wellTransfer = wellTransferFromString(
                    // Same default as setupPropertyTree ships, so a JSON that omits
                    // the key gets the same preconditioner as the built-in setup.
                    prm.get("well_transfer", std::string {"classic"}));
                const auto diagonal = wellCoarseDiagonalFromString(
                    prm.get("well_coarse_diagonal", std::string {"contract_d"}));
                cprwStage_ = std::make_unique<CprwStage>(S_,
                                                         coarseprm,
                                                         pressureIndex_,
                                                         wellTransfer,
                                                         resComm_,
                                                         diagonal,
                                                         prm.get("verbosity", 0));
                cprwStage_->buildStructure(weights_);
            } else {
                if constexpr (isParallel) {
                    resSolver_ = std::make_unique<ResFlexibleSolverType>(
                        *rop_, *resComm_, resprm, resWeightCalc, pressureIndex_);
                } else {
                    resSolver_ = std::make_unique<ResFlexibleSolverType>(
                        *rop_, resprm, resWeightCalc, pressureIndex_);
                }
            }
        }

        initWellSolver();
    }

    void initWorkVectors()
    {
        resizeReservoirWorkVectors();
        resizeWellWorkVectors();
    }

    void resizeReservoirWorkVectors()
    {
        const auto numRes = S_.A->N();
        Backend::resize(resSol_, numRes, numResDofs);
        Backend::resize(dresSol_, numRes, numResDofs);
        Backend::resize(tmp_resRes_, numRes, numResDofs);
        Backend::resize(resRes_, numRes, numResDofs);
    }

    void resizeWellWorkVectors()
    {
        const auto numWell = Backend::cpuWellMatrix(S_).N();
        Backend::resize(wSol_, numWell, numWellDofs);
        Backend::resize(dwSol_, numWell, numWellDofs);
        Backend::resize(tmp_wRes_, numWell, numWellDofs);
        Backend::resize(wRes_, numWell, numWellDofs);
        if constexpr (isGpu) {
            hostWRes_->resize(numWell);
            hostDwSol_->resize(numWell);
        }
    }
};

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONER_HEADER_INCLUDED
