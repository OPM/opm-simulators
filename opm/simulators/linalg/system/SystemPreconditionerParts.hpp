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

#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <functional>
#include <memory>
#include <type_traits>

namespace Opm
{

// --------------------------------------------------------------------------
// Building blocks for preconditioners on the coupled system
//
//     S = [ A  C ]
//         [ B  D ]
//
// Each part below is a PreconditionerWithUpdate over SystemVector, so they can
// be composed and handed to anything that takes a preconditioner -- the fine
// smoother of a two-level method, or the whole preconditioner.
//
// The composition is multiplicative: every part recomputes the defect left by
// the part before it, rather than working from the original right-hand side.
// That is what makes "reservoir smoother followed by a well solve" a
// block Gauss-Seidel sweep rather than a block Jacobi one.
// --------------------------------------------------------------------------

// Solve the well block only:  v[_1] += D^-1 (d[_1] - B v[_0] - D v[_1]).
//
// The defect is recomputed from the correction accumulated so far, so this can
// be appended to anything -- a reservoir smoother, or a coarse pressure
// correction -- and it will clean up the well equations that step disturbed.
template <class Scalar>
class SystemWellSolve
{
public:
    using WellOperator = Dune::MatrixAdapter<WWMatrix<Scalar>, WellVector<Scalar>, WellVector<Scalar>>;
    using WellSolver = Dune::FlexibleSolver<WellOperator>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemWellSolve(const SystemMatrix<Scalar>& S, const PropertyTree& prm)
        : S_(S)
        , prm_(prm)
    {
        build();
    }

    // Rebuild against the current D. Needed when the well structure changes,
    // since D changes dimension.
    void build()
    {
        wop_ = std::make_unique<WellOperator>(*S_.D);
        std::function<WellVector<Scalar>()> noWeights;
        wellSolver_ = std::make_unique<WellSolver>(*wop_, prm_, noWeights, dummyPressureIndex);
        work_.resize(S_.D->N());
        corr_.resize(S_.D->N());
    }

    void update()
    {
        wellSolver_->preconditioner().update();
    }

    // v is the correction so far, d the original right-hand side.
    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d)
    {
        // work = d[_1] - B v[_0] - D v[_1]
        work_ = d[_1];
        S_.B->mmv(v[_0], work_);
        S_.D->mmv(v[_1], work_);

        corr_ = 0.0;
        Dune::InverseOperatorResult result;
        wellSolver_->apply(corr_, work_, result);
        v[_1] += corr_;
    }

private:
    // The well solver never uses a pressure index; pass something that will
    // fail loudly if it ever is used.
    static constexpr std::size_t dummyPressureIndex = static_cast<std::size_t>(-1);

    const SystemMatrix<Scalar>& S_;
    PropertyTree prm_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<WellSolver> wellSolver_;
    WellVector<Scalar> work_;
    WellVector<Scalar> corr_;
};

// Apply a reservoir-only solver to the reservoir block:
//     v[_0] += R^-1 (d[_0] - A v[_0] - C v[_1]).
//
// Lifts anything that acts on ResVector (ILU0, a CPR solve, ...) into a step on
// the coupled system. On its own it leaves the well unknowns untouched; pair it
// with SystemWellSolve to get a full sweep.
template <class Scalar, class ResOp, class ResComm = Dune::Amg::SequentialInformation>
class SystemReservoirSolve
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    using ResSolver = Dune::FlexibleSolver<ResOp>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemReservoirSolve(const SystemMatrix<Scalar>& S,
                         const PropertyTree& prm,
                         const std::function<ResVector<Scalar>()>& weightsCalculator,
                         const int pressureIndex,
                         const ResComm* comm = nullptr)
        : S_(S)
        , prm_(prm)
        , weightsCalculator_(weightsCalculator)
        , pressureIndex_(pressureIndex)
        , comm_(comm)
    {
        build();
    }

    void build()
    {
        if constexpr (isParallel) {
            rop_ = std::make_unique<ResOp>(*S_.A, *comm_);
            resSolver_ = std::make_unique<ResSolver>(*rop_, *comm_, prm_,
                                                     weightsCalculator_, pressureIndex_);
        } else {
            rop_ = std::make_unique<ResOp>(*S_.A);
            resSolver_ = std::make_unique<ResSolver>(*rop_, prm_,
                                                     weightsCalculator_, pressureIndex_);
        }
        work_.resize(S_.A->N());
        corr_.resize(S_.A->N());
    }

    void update()
    {
        resSolver_->preconditioner().update();
    }

    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d)
    {
        // work = d[_0] - A v[_0] - C v[_1]
        work_ = d[_0];
        S_.A->mmv(v[_0], work_);
        S_.C->mmv(v[_1], work_);
        if constexpr (isParallel) {
            comm_->copyOwnerToAll(work_, work_);
        }

        corr_ = 0.0;
        Dune::InverseOperatorResult result;
        resSolver_->apply(corr_, work_, result);
        if constexpr (isParallel) {
            comm_->copyOwnerToAll(corr_, corr_);
        }
        v[_0] += corr_;
    }

private:
    const SystemMatrix<Scalar>& S_;
    PropertyTree prm_;
    std::function<ResVector<Scalar>()> weightsCalculator_;
    int pressureIndex_ = 0;
    const ResComm* comm_ = nullptr;

    std::unique_ptr<ResOp> rop_;
    std::unique_ptr<ResSolver> resSolver_;
    ResVector<Scalar> work_;
    ResVector<Scalar> corr_;
};

// A reservoir solve followed by a well solve -- the "reservoir smoother plus a
// well solve" pairing. Both halves recompute their own defect, so this is one
// multiplicative sweep over the two blocks.
template <class Scalar, class ResOp, class ResComm = Dune::Amg::SequentialInformation>
class SystemReservoirWellSweep
{
public:
    using Reservoir = SystemReservoirSolve<Scalar, ResOp, ResComm>;

    SystemReservoirWellSweep(const SystemMatrix<Scalar>& S,
                             const PropertyTree& prm,
                             const std::function<ResVector<Scalar>()>& weightsCalculator,
                             const int pressureIndex,
                             const ResComm* comm = nullptr)
        : reservoir_(S, prm.get_child("reservoir_smoother"), weightsCalculator, pressureIndex, comm)
        , well_(S, prm.get_child("well_solver"))
    {
    }

    void buildForChangedWellStructure()
    {
        well_.build();
    }

    void update()
    {
        reservoir_.update();
        well_.update();
    }

    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d)
    {
        reservoir_.apply(v, d);
        well_.apply(v, d);
    }

private:
    Reservoir reservoir_;
    SystemWellSolve<Scalar> well_;
};

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONERPARTS_HEADER_INCLUDED
