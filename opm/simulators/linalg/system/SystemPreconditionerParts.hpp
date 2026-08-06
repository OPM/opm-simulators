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
#ifndef OPM_SYSTEMPRECONDITIONERPARTS_HEADER_INCLUDED
#define OPM_SYSTEMPRECONDITIONERPARTS_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemCprwPressureStage.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace Opm
{

// --------------------------------------------------------------------------
// Building blocks for a preconditioner on the coupled system
//
//     S = [ A  C ]
//         [ B  D ]
//
// A sweep carries two things: the correction accumulated so far, and the
// residual that correction leaves. Each step reads the residual, produces a
// correction for one block, and subtracts that correction's effect from both
// halves of the residual.
//
// The residual is maintained incrementally rather than recomputed as d - S v.
// The two are equal in exact arithmetic but not bit for bit, and the point of
// this decomposition is that composing the steps reproduces the hand-written
// 3-stage SystemPreconditioner exactly -- see GeneralSystemPreconditioner.
//
// Composition is therefore multiplicative: every step sees what the step
// before it left behind. "Reservoir smoother followed by a well solve" is a
// block Gauss-Seidel sweep, not a block Jacobi one.
// --------------------------------------------------------------------------

template <class Scalar>
struct SystemSweepState
{
    // Correction accumulated so far.
    SystemVector<Scalar>* v = nullptr;
    // Residual left by that correction, maintained incrementally.
    SystemVector<Scalar>* res = nullptr;
};

template <class Scalar>
class SystemSweepStep
{
public:
    virtual ~SystemSweepStep() = default;

    virtual void apply(SystemSweepState<Scalar>& state) = 0;

    // Refresh against changed matrix values, pattern unchanged.
    virtual void update() = 0;

    // Rebuild against a changed well structure, which changes D's dimension.
    // A changed well structure comes with changed matrix values, so a step
    // that owns nothing sized by the wells still has to refresh itself; only
    // steps that must be rebuilt outright override this.
    virtual void rebuildForChangedWellStructure()
    {
        update();
    }
};

// --------------------------------------------------------------------------
// Solve the well block:  corr = D^-1 res[_1],  v[_1] += corr,
// res[_0] -= C corr,  res[_1] -= D corr.
//
// Appending this to anything -- a coarse pressure correction, a reservoir
// smoother -- cleans up the well equations that step disturbed.
// --------------------------------------------------------------------------
template <class Scalar>
class SystemWellStep : public SystemSweepStep<Scalar>
{
public:
    using WellOperator = Dune::MatrixAdapter<WWMatrix<Scalar>, WellVector<Scalar>, WellVector<Scalar>>;
    using WellSolver = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemWellStep(const SystemMatrix<Scalar>& S, const PropertyTree& prm)
        : S_(S)
        , prm_(prm)
    {
        rebuildForChangedWellStructure();
    }

    void rebuildForChangedWellStructure() override
    {
        wop_ = std::make_unique<WellOperator>(*S_.D);
        std::function<WellVector<Scalar>()> noWeights;
        solver_ = std::make_unique<WellSolver>(*wop_, prm_, noWeights, dummyPressureIndex);
        rhs_.resize(S_.D->N());
        corr_.resize(S_.D->N());
    }

    void update() override
    {
        solver_->preconditioner().update();
    }

    void apply(SystemSweepState<Scalar>& state) override
    {
        auto& v = *state.v;
        auto& res = *state.res;

        rhs_ = res[_1];
        corr_ = 0.0;
        Dune::InverseOperatorResult result;
        solver_->apply(corr_, rhs_, result);

        v[_1] += corr_;
        S_.C->mmv(corr_, res[_0]);
        S_.D->mmv(corr_, res[_1]);
    }

private:
    // The well solver never uses a pressure index; pass something that would
    // fail loudly if it ever did.
    static constexpr std::size_t dummyPressureIndex = static_cast<std::size_t>(-1);

    const SystemMatrix<Scalar>& S_;
    PropertyTree prm_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<WellSolver> solver_;
    WellVector<Scalar> rhs_;
    WellVector<Scalar> corr_;
};

// --------------------------------------------------------------------------
// Apply a reservoir-only solver to the reservoir block:
//     corr = R^-1 res[_0],  v[_0] += corr,
//     res[_0] -= A corr,  res[_1] -= B corr.
//
// Lifts anything acting on a ResVector (ILU0, a CPR solve, ...) into a step on
// the coupled system. It leaves the well unknowns untouched; pair it with
// SystemWellStep for a full sweep.
//
// In parallel the residual handed to the solver is made consistent first; the
// correction is not, and the caller synchronises the accumulated reservoir
// correction once at the end instead.
// --------------------------------------------------------------------------
template <class Scalar, class ResOp, class ResComm = Dune::Amg::SequentialInformation>
class SystemReservoirStep : public SystemSweepStep<Scalar>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    using ResSolver = Dune::FlexibleSolver<ResOp>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemReservoirStep(const SystemMatrix<Scalar>& S,
                        const PropertyTree& prm,
                        const std::function<ResVector<Scalar>()>& weightsCalculator,
                        const int pressureIndex,
                        const ResComm* comm = nullptr)
        : S_(S)
    {
        if constexpr (isParallel) {
            comm_ = comm;
            rop_ = std::make_unique<ResOp>(*S_.A, *comm_);
            solver_ = std::make_unique<ResSolver>(*rop_, *comm_, prm,
                                                  weightsCalculator, pressureIndex);
        } else {
            rop_ = std::make_unique<ResOp>(*S_.A);
            solver_ = std::make_unique<ResSolver>(*rop_, prm,
                                                  weightsCalculator, pressureIndex);
        }
        rhs_.resize(S_.A->N());
        corr_.resize(S_.A->N());
    }

    void update() override
    {
        solver_->preconditioner().update();
    }

    void apply(SystemSweepState<Scalar>& state) override
    {
        auto& v = *state.v;
        auto& res = *state.res;

        rhs_ = res[_0];
        if constexpr (isParallel) {
            comm_->copyOwnerToAll(rhs_, rhs_);
        }

        corr_ = 0.0;
        Dune::InverseOperatorResult result;
        solver_->apply(corr_, rhs_, result);

        v[_0] += corr_;
        S_.A->mmv(corr_, res[_0]);
        S_.B->mmv(corr_, res[_1]);
    }

private:
    const SystemMatrix<Scalar>& S_;
    const ResComm* comm_ = nullptr;
    std::unique_ptr<ResOp> rop_;
    std::unique_ptr<ResSolver> solver_;
    ResVector<Scalar> rhs_;
    ResVector<Scalar> corr_;
};

// --------------------------------------------------------------------------
// The CPRW pressure stage as a sweep step: restrict the coupled residual to
// the scalar pressure system, solve it, prolong back.
//
// Whether the coarse correction reaches the well unknowns at all is the
// stage's own business (well_transfer); when it does not, there is nothing to
// subtract from the residual through C and D.
// --------------------------------------------------------------------------
template <class Scalar, class ResComm = Dune::Amg::SequentialInformation>
class SystemCprwStep : public SystemSweepStep<Scalar>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    using Stage = SystemCprwPressureStage<Scalar, ResComm>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemCprwStep(const SystemMatrix<Scalar>& S,
                   const PropertyTree& coarseSolverPrm,
                   const std::function<SystemVector<Scalar>()>& weightsCalculator,
                   const int pressureIndex,
                   const WellTransfer wellTransfer,
                   const WellCoarseDiagonal diagonal,
                   const int verbosity,
                   const ResComm* comm = nullptr)
        : S_(S)
        , weightsCalculator_(weightsCalculator)
    {
        if (!weightsCalculator_) {
            OPM_THROW(std::invalid_argument,
                      "The CPRW pressure stage (add_wells) needs a weights calculator, but "
                      "none was configured. Set the reservoir solver's weight_type.");
        }
        if constexpr (isParallel) {
            comm_ = comm;
        }
        stage_ = std::make_unique<Stage>(S_, coarseSolverPrm, pressureIndex,
                                         wellTransfer, comm_, diagonal, verbosity);
        weights_ = weightsCalculator_();
        stage_->buildStructure(weights_);
        rhs_.resize(S_.A->N());
        corrRes_.resize(S_.A->N());
        corrWell_.resize(S_.D->N());
    }

    void update() override
    {
        weights_ = weightsCalculator_();
        stage_->update(weights_);
    }

    void rebuildForChangedWellStructure() override
    {
        // The coarse system carries one unknown per well, so a changed well
        // structure changes its dimension and pattern.
        weights_ = weightsCalculator_();
        stage_->buildStructure(weights_);
        corrWell_.resize(S_.D->N());
    }

    void apply(SystemSweepState<Scalar>& state) override
    {
        auto& v = *state.v;
        auto& res = *state.res;

        rhs_ = res[_0];
        if constexpr (isParallel) {
            comm_->copyOwnerToAll(rhs_, rhs_);
        }

        corrRes_ = 0.0;
        corrWell_ = 0.0;
        stage_->apply(rhs_, res[_1], weights_, corrRes_, corrWell_);

        v[_0] += corrRes_;
        S_.A->mmv(corrRes_, res[_0]);
        S_.B->mmv(corrRes_, res[_1]);

        if (stage_->prolongatesWellPressure()) {
            v[_1] += corrWell_;
            S_.C->mmv(corrWell_, res[_0]);
            S_.D->mmv(corrWell_, res[_1]);
        }
    }

private:
    const SystemMatrix<Scalar>& S_;
    const ResComm* comm_ = nullptr;
    std::function<SystemVector<Scalar>()> weightsCalculator_;
    SystemVector<Scalar> weights_;
    std::unique_ptr<Stage> stage_;
    ResVector<Scalar> rhs_;
    ResVector<Scalar> corrRes_;
    WellVector<Scalar> corrWell_;
};


// --------------------------------------------------------------------------
// An ordered list of steps applied as one multiplicative sweep, presented as a
// Dune preconditioner so it can serve as the fine smoother of a two-level
// method.
//
// Dune hands a smoother a single (correction, defect) pair, so a sweep of more
// than one block has to carry its own defect internally. That is what this
// does: it starts from the defect it is given, and each step reduces it.
// --------------------------------------------------------------------------
template <class Scalar>
class SystemSweepPreconditioner
    : public Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>
{
public:
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    explicit SystemSweepPreconditioner(std::vector<std::unique_ptr<SystemSweepStep<Scalar>>> steps,
                                       const bool parallel = false)
        : steps_(std::move(steps))
        , parallel_(parallel)
    {
    }

    void pre(SystemVector<Scalar>&, SystemVector<Scalar>&) override {}
    void post(SystemVector<Scalar>&) override {}

    Dune::SolverCategory::Category category() const override
    {
        return parallel_ ? Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

    void update() override
    {
        for (const auto& step : steps_) {
            step->update();
        }
    }

    void rebuildForChangedWellStructure()
    {
        for (const auto& step : steps_) {
            step->rebuildForChangedWellStructure();
        }
    }

    bool empty() const
    {
        return steps_.empty();
    }

    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d) override
    {
        v[_0].resize(d[_0].size());
        v[_1].resize(d[_1].size());
        v[_0] = 0.0;
        v[_1] = 0.0;
        res_[_0] = d[_0];
        res_[_1] = d[_1];

        SystemSweepState<Scalar> state{&v, &res_};
        for (const auto& step : steps_) {
            step->apply(state);
        }
    }

private:
    std::vector<std::unique_ptr<SystemSweepStep<Scalar>>> steps_;
    bool parallel_ = false;
    SystemVector<Scalar> res_;
};

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONERPARTS_HEADER_INCLUDED
