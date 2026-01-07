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
#include <opm/simulators/linalg/system/SystemMatrixWellOperator.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <functional>
#include <limits>
#include <memory>
#include <type_traits>


namespace Opm
{

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
// Templated on scalar type, reservoir operator and communication types to
// unify sequential and parallel implementations. The 3-stage algorithm:
//   1. Reservoir CPR (or CPRW when WithWellCoupling=true) solve
//   2. Well solve + reservoir smoothing
//   3. Final well solve
//
// When WithWellCoupling=true (system_cprw), the internal reservoir operator
// is WellModelMatrixAdapter<RRMatrix,ResVector,ResVector> which exposes the
// B/C/D blocks via SystemMatrixWellOperator, allowing the standard CPRW
// preconditioner to build the extended (nCells+nWells) scalar pressure system.
//
// For parallel runs, copyOwnerToAll synchronises overlap DOFs before
// each reservoir sub-solve.
template <class Scalar, class ResOp,
          class ResComm = Dune::Amg::SequentialInformation,
          bool WithWellCoupling = false>
class SystemPreconditioner : public Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    // When WithWellCoupling=true, the internal resSolver_ uses the "cprw"
    // preconditioner type, which is registered for WellModelMatrixAdapter
    // (sequential) and WellModelGhostLastMatrixAdapter (parallel, ghost-last
    // ordered). Select the matching adapter by parallelism.
    using SeqWellResOp = WellModelMatrixAdapter<RRMatrix<Scalar>,
                                                ResVector<Scalar>,
                                                ResVector<Scalar>>;
#if HAVE_MPI
    using ParWellResOp = WellModelGhostLastMatrixAdapter<RRMatrix<Scalar>,
                                                         ResVector<Scalar>,
                                                         ResVector<Scalar>,
                                                         true>;
#else
    // Unreachable when MPI is off (isParallel is always false); kept so the
    // std::conditional_t below stays well-formed.
    using ParWellResOp = SeqWellResOp;
#endif
    using WellResOp = std::conditional_t<isParallel, ParWellResOp, SeqWellResOp>;
    using InternalResOp = std::conditional_t<WithWellCoupling, WellResOp, ResOp>;
    using ResFlexibleSolverType = Dune::FlexibleSolver<InternalResOp>;
    // Smoother always uses ResOp (plain matrix adapter): standard preconditioners such as
    // paroverilu0 are not registered in the factory for WellModelMatrixAdapter.
    using SmoothFlexibleSolverType = Dune::FlexibleSolver<ResOp>;
    using WellOperatorExtra = SystemMatrixWellOperator<Scalar>;

    using WellOperator = Dune::MatrixAdapter<WWMatrix<Scalar>, WellVector<Scalar>, WellVector<Scalar>>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    // Sequential constructor (enabled only for non-parallel specializations).
    SystemPreconditioner(const SystemMatrix<Scalar>& S,
                         const std::function<ResVector<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm)
        requires (!isParallel)
        : S_(S)
        , pressureIndex_(pressureIndex)
    {
        initSubSolvers(prm, weightsCalculator);
        initWorkVectors();
    }

    // Parallel constructor (enabled only for parallel specializations).
    SystemPreconditioner(const SystemMatrix<Scalar>& S,
                         const std::function<ResVector<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm,
                         const ResComm& resComm)
        requires (isParallel)
        : S_(S)
        , resComm_(&resComm)
        , pressureIndex_(pressureIndex)
    {
        initSubSolvers(prm, weightsCalculator);
        initWorkVectors();
    }

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

    void update() override
    {
        resSolver_->preconditioner().update();
        resSmoother_->preconditioner().update();
        wellSolver_->preconditioner().update();
    }

    void updateForChangedWellStructure()
    {
        if constexpr (WithWellCoupling) {
            // Well structure change may alter coarse matrix dimensions (nWells changed):
            // recreate the well operator and reservoir solver from scratch.
            initResOpAndSolver();
        } else {
            resSolver_->preconditioner().update();
        }
        resSmoother_->preconditioner().update();
        initWellSolver();
        resizeWellWorkVectors();
    }

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
     void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d) override
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

private:
    const SystemMatrix<Scalar>& S_;
    const ResComm* resComm_ = nullptr;
    int pressureIndex_ = 0;
    static constexpr int dummyWellPressureIndex = std::numeric_limits<int>::min();
    Opm::PropertyTree wellprm_;

    // Used only when WithWellCoupling=true: stored to allow rebuilding resSolver_
    // when well structure changes (nWells may change, altering coarse matrix dims).
    Opm::PropertyTree resSolverPrm_;
    std::function<ResVector<Scalar>()> resSolverWeightsCalc_;

    // When WithWellCoupling=true, systemWellOper_ must outlive rop_ (held by ref),
    // and rop_ must outlive resSolver_. Declaration order below respects this so
    // destruction (reverse order) tears them down safely.
    std::unique_ptr<WellOperatorExtra> systemWellOper_;
    std::unique_ptr<InternalResOp> rop_;
    // When WithWellCoupling=true, rop_ is WellModelMatrixAdapter; the smoother needs a
    // plain ResOp so that standard preconditioners (paroverilu0 etc.) can be looked up.
    std::unique_ptr<ResOp> srop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<SmoothFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    WellVector<Scalar> wSol_;
    ResVector<Scalar> resSol_;
    ResVector<Scalar> dresSol_;
    WellVector<Scalar> dwSol_;
    ResVector<Scalar> tmp_resRes_;
    WellVector<Scalar> tmp_wRes_;
    ResVector<Scalar> resRes_;
    WellVector<Scalar> wRes_;

    void syncResVector(ResVector<Scalar>& v)
    {
        if constexpr (isParallel) {
            resComm_->copyOwnerToAll(v, v);
        }
    }

    void initWellSolver()
    {
        wop_ = std::make_unique<WellOperator>(S_[_1][_1]);
        std::function<WellVector<Scalar>()> weightsCalculatorWell;
        wellSolver_ = std::make_unique<WellFlexibleSolverType>(
            *wop_, wellprm_, weightsCalculatorWell, dummyWellPressureIndex);
    }

    // Build the internal reservoir operator and solver for the CPRW
    // (WithWellCoupling=true) case. The reservoir solver runs on a well-model
    // adapter wrapping A_rr plus the well coupling exposed by
    // SystemMatrixWellOperator, so it can use the "cprw" preconditioner type:
    // WellModelMatrixAdapter serially, WellModelGhostLastMatrixAdapter in parallel.
    void initResOpAndSolver()
    {
        if constexpr (WithWellCoupling) {
            systemWellOper_ = std::make_unique<WellOperatorExtra>(S_, pressureIndex_);
            if constexpr (isParallel) {
                // WellModelGhostLastMatrixAdapter assumes interior rows precede
                // ghost rows; interiorSize is the owner-row count derived from the
                // reservoir communication's index set.
                const std::size_t interiorSize = computeInteriorSize(*resComm_);
                rop_ = std::make_unique<InternalResOp>(S_[_0][_0], *systemWellOper_, interiorSize);
                resSolver_ = std::make_unique<ResFlexibleSolverType>(
                    *rop_, *resComm_, resSolverPrm_, resSolverWeightsCalc_, pressureIndex_);
            } else {
                rop_ = std::make_unique<InternalResOp>(S_[_0][_0], *systemWellOper_);
                resSolver_ = std::make_unique<ResFlexibleSolverType>(
                    *rop_, resSolverPrm_, resSolverWeightsCalc_, pressureIndex_);
            }
        }
    }

    // Owner-row count for the ghost-last reservoir adapter, derived from the
    // communication index set (mirrors GhostLastMatrixAdapter::setInteriorSize in
    // WellOperators.hpp). No grid is available here, so the comm is the source.
    static std::size_t computeInteriorSize(const ResComm& comm)
    {
        if constexpr (isParallel) {
            const auto& indexSet = comm.indexSet();
            if (indexSet.size() == 0)
                return 0;
            std::size_t is = 0;
            for (auto idx = indexSet.begin(); idx != indexSet.end(); ++idx) {
                if (idx->local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                    is = std::max(is, static_cast<std::size_t>(idx->local().local()));
                }
            }
            return is + 1;
        } else {
            return 0;
        }
    }

    void initSubSolvers(const Opm::PropertyTree& prm,
                        const std::function<ResVector<Scalar>()>& weightsCalculator)
    {
        auto resprm = prm.get_child("reservoir_solver");
        auto resprmsmoother = prm.get_child("reservoir_smoother");
        wellprm_ = prm.get_child("well_solver");

        if constexpr (WithWellCoupling) {
            resSolverPrm_ = resprm;
            resSolverWeightsCalc_ = weightsCalculator;
            initResOpAndSolver();
        } else if constexpr (isParallel) {
            rop_ = std::make_unique<ResOp>(S_[_0][_0], *resComm_);
            resSolver_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, *resComm_, resprm, weightsCalculator, pressureIndex_);
        } else {
            rop_ = std::make_unique<ResOp>(S_[_0][_0]);
            resSolver_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, resprm, weightsCalculator, pressureIndex_);
        }

        // Smoother always operates on a plain ResOp (A_rr only).
        if constexpr (WithWellCoupling) {
            // rop_ is a well-model adapter; paroverilu0 and other standard preconditioners
            // are not registered for that type. Use a plain ResOp wrapper on A_rr instead.
            // In parallel ResOp is OverlappingSchwarzOperator (full matrix + comm
            // projection), so no ghost-last assumption is needed for the smoother.
            if constexpr (isParallel) {
                srop_ = std::make_unique<ResOp>(S_[_0][_0], *resComm_);
                resSmoother_ = std::make_unique<SmoothFlexibleSolverType>(
                    *srop_, *resComm_, resprmsmoother, weightsCalculator, pressureIndex_);
            } else {
                srop_ = std::make_unique<ResOp>(S_[_0][_0]);
                resSmoother_ = std::make_unique<SmoothFlexibleSolverType>(
                    *srop_, resprmsmoother, weightsCalculator, pressureIndex_);
            }
        } else if constexpr (isParallel) {
            resSmoother_ = std::make_unique<SmoothFlexibleSolverType>(
                *rop_, *resComm_, resprmsmoother, weightsCalculator, pressureIndex_);
        } else {
            resSmoother_ = std::make_unique<SmoothFlexibleSolverType>(
                *rop_, resprmsmoother, weightsCalculator, pressureIndex_);
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
        const auto numRes = S_[_0][_0].N();
        resSol_.resize(numRes);
        dresSol_.resize(numRes);
        tmp_resRes_.resize(numRes);
        resRes_.resize(numRes);
    }

    void resizeWellWorkVectors()
    {
        const auto numWell = S_[_1][_1].N();
        wSol_.resize(numWell);
        dwSol_.resize(numWell);
        tmp_wRes_.resize(numWell);
        wRes_.resize(numWell);
    }
};

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONER_HEADER_INCLUDED
