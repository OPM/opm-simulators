#pragma once

#include "MultiComm.hpp"
#include "SystemTypes.hpp"

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm
{

// Reservoir operator/comm types used as template arguments.
template<typename Scalar>
using SeqResOperatorT = Dune::MatrixAdapter<RRMatrixT<Scalar>, ResVectorT<Scalar>, ResVectorT<Scalar>>;

#if HAVE_MPI
using ParResComm = Dune::OwnerOverlapCopyCommunication<int, int>;
template<typename Scalar>
using ParResOperatorT = Dune::OverlappingSchwarzOperator<RRMatrixT<Scalar>, ResVectorT<Scalar>, ResVectorT<Scalar>, ParResComm>;
#endif

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
class SystemPreconditioner : public Dune::PreconditionerWithUpdate<SystemVectorT<Scalar>, SystemVectorT<Scalar>>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOp>;
    using WellOperator = Dune::MatrixAdapter<WWMatrixT<Scalar>, WellVectorT<Scalar>, WellVectorT<Scalar>>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    // Sequential constructor (enabled only for non-parallel specializations).
    template <bool P = isParallel, std::enable_if_t<!P, int> = 0>
    SystemPreconditioner(const SystemMatrixT<Scalar>& S,
                         const std::function<ResVectorT<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm)
        : S_(S)
    {
        initSubSolvers(prm, weightsCalculator, pressureIndex);
        initWorkVectors();
    }

    // Parallel constructor (enabled only for parallel specializations).
    template <bool P = isParallel, std::enable_if_t<P, int> = 0>
    SystemPreconditioner(const SystemMatrixT<Scalar>& S,
                         const std::function<ResVectorT<Scalar>()>& weightsCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm,
                         const ResComm& resComm)
        : S_(S)
        , resComm_(&resComm)
    {
        initSubSolvers(prm, weightsCalculator, pressureIndex);
        initWorkVectors();
    }

    void pre(SystemVectorT<Scalar>&, SystemVectorT<Scalar>&) override
    {
    }
    void post(SystemVectorT<Scalar>&) override
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

    bool hasPerfectUpdate() const override
    {
        return true;
    }

    void apply(SystemVectorT<Scalar>& v, const SystemVectorT<Scalar>& d) override
    {
        const auto& A = S_[_0][_0];
        const auto& C = S_[_1][_0];
        const auto& B = S_[_0][_1];
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
            A.mmv(dresSol_, resRes_);
            C.mmv(dresSol_, wRes_);
        }

        // Stage 2: Well solve + reservoir system smoothing
        {
            Dune::InverseOperatorResult well_result;
            dwSol_ = 0.0;
            tmp_wRes_ = wRes_;
            wellSolver_->apply(dwSol_, tmp_wRes_, well_result);
            wSol_ += dwSol_;
            B.mmv(dwSol_, resRes_);
            D.mmv(dwSol_, wRes_);

            Dune::InverseOperatorResult res_result;
            dresSol_ = 0.0;
            tmp_resRes_ = resRes_;
            syncResVector(tmp_resRes_);
            resSmoother_->apply(dresSol_, tmp_resRes_, res_result);
            resSol_ += dresSol_;
            C.mmv(dresSol_, wRes_);
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
    const SystemMatrixT<Scalar>& S_;
    const ResComm* resComm_ = nullptr;

    std::unique_ptr<ResOp> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<ResFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    WellVectorT<Scalar> wSol_;
    ResVectorT<Scalar> resSol_;
    ResVectorT<Scalar> dresSol_;
    WellVectorT<Scalar> dwSol_;
    ResVectorT<Scalar> tmp_resRes_;
    WellVectorT<Scalar> tmp_wRes_;
    ResVectorT<Scalar> resRes_;
    WellVectorT<Scalar> wRes_;

    void syncResVector(ResVectorT<Scalar>& v)
    {
        if constexpr (isParallel) {
            resComm_->copyOwnerToAll(v, v);
        }
    }

    void initSubSolvers(const Opm::PropertyTree& prm,
                        const std::function<ResVectorT<Scalar>()>& weightsCalculator,
                        int pressureIndex)
    {
        auto resprm = prm.get_child("reservoir_solver");
        auto resprmsmoother = prm.get_child("reservoir_smoother");
        auto wellprm = prm.get_child("well_solver");

        if constexpr (isParallel) {
            rop_ = std::make_unique<ResOp>(S_[_0][_0], *resComm_);
            resSolver_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, *resComm_, resprm, weightsCalculator, pressureIndex);
            resSmoother_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, *resComm_, resprmsmoother, weightsCalculator, pressureIndex);
        } else {
            rop_ = std::make_unique<ResOp>(S_[_0][_0]);
            resSolver_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, resprm, weightsCalculator, pressureIndex);
            resSmoother_ = std::make_unique<ResFlexibleSolverType>(
                *rop_, resprmsmoother, weightsCalculator, pressureIndex);
        }

        auto wop = std::make_unique<WellOperator>(S_[_1][_1]);
        std::function<WellVectorT<Scalar>()> weightsCalculatorWell;
        wellSolver_ = std::make_unique<WellFlexibleSolverType>(
            *wop, wellprm, weightsCalculatorWell, pressureIndex);
        wop_ = std::move(wop);
    }

    void initWorkVectors()
    {
        const auto numRes = S_[_0][_0].N();
        const auto numWell = S_[_1][_1].N();
        wSol_.resize(numWell);
        resSol_.resize(numRes);
        dresSol_.resize(numRes);
        dwSol_.resize(numWell);
        tmp_resRes_.resize(numRes);
        tmp_wRes_.resize(numWell);
        resRes_.resize(numRes);
        wRes_.resize(numWell);
    }
};

} // namespace Opm
