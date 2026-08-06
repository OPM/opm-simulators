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
#ifndef OPM_GENERALSYSTEMPRECONDITIONER_HEADER_INCLUDED
#define OPM_GENERALSYSTEMPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditionerParts.hpp>
#include <opm/simulators/linalg/system/SystemPressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <opm/common/ErrorMacros.hpp>

#include <dune/istl/paamg/pinfo.hh>

#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace Opm
{

// --------------------------------------------------------------------------
// A preconditioner for the coupled (reservoir, well) system assembled from
// named parts:
//
//   coarse_solver { reservoir_solver, well_solver }
//   smoother      { reservoir_smoother, well_solver }
//
// When the coarse solver is the CPRW pressure stage (add_wells), this is a
// genuine two-level method on the system and is run by Dune's
// TwoLevelMethodCpr: SystemPressureBhpTransferPolicy is the transfer, the
// scalar pressure system of dimension nCells + nWells is the coarse level, and
// everything after the coarse correction is the fine smoother, with
// preSteps = 0 and postSteps = 1.
//
// Dune hands a smoother one (correction, defect) pair, so the parts that
// follow the coarse correction -- the coarse solver's own well solve, then the
// smoother's parts -- are composed into a single SystemSweepPreconditioner
// that carries its defect internally. The sequence applied is therefore
//
//     coarse -> well -> reservoir smoother -> well
//
// which is what SystemPreconditioner does, and the two agree bit for bit.
// That works because TwoLevelMethodCpr's bookkeeping is the same one the
// sweep uses: postsmooth() reduces the defect by applyscaleadd(-1, lhs, rhs),
// which on SystemMatrix accumulates A,C into the reservoir half and B,D into
// the well half -- the same operations in the same order as the hand-written
// version's four mmv calls.
//
// Without add_wells there is no system-wide coarse space: the reservoir solve
// is an ordinary block solve, so all the parts go into one sweep and no
// two-level method is built.
// --------------------------------------------------------------------------
template <class Scalar, class ResOp, class ResComm = Dune::Amg::SequentialInformation>
class GeneralSystemPreconditioner
    : public Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>
{
public:
    static constexpr bool isParallel = !std::is_same_v<ResComm, Dune::Amg::SequentialInformation>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    using ReservoirStep = SystemReservoirStep<Scalar, ResOp, ResComm>;
    using WellStep = SystemWellStep<Scalar>;
    using Sweep = SystemSweepPreconditioner<Scalar>;

    // The fine level of the two-level method is the whole coupled system.
    using FineOperator = Dune::MatrixAdapter<SystemMatrix<Scalar>,
                                             SystemVector<Scalar>,
                                             SystemVector<Scalar>>;
    using TransferPolicy = SystemPressureBhpTransferPolicy<FineOperator, ResComm, Scalar>;
    using CoarseOperator = typename TransferPolicy::CoarseOperator;
    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperator,
                                                               Dune::FlexibleSolver<CoarseOperator>,
                                                               TransferPolicy>;
    using TwoLevel = Dune::Amg::TwoLevelMethodCpr<
        FineOperator,
        CoarseSolverPolicy,
        Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>>;

    GeneralSystemPreconditioner(const SystemMatrix<Scalar>& S,
                                const std::function<SystemVector<Scalar>()>& weightsCalculator,
                                int pressureIndex,
                                const PropertyTree& prm)
        requires(!isParallel)
        : S_(S)
        , pressureIndex_(pressureIndex)
        , prm_(prm)
        , weightsCalculator_(weightsCalculator)
    {
        build();
    }

    GeneralSystemPreconditioner(const SystemMatrix<Scalar>& S,
                                const std::function<SystemVector<Scalar>()>& weightsCalculator,
                                int pressureIndex,
                                const PropertyTree& prm,
                                const ResComm& resComm)
        requires(isParallel)
        : S_(S)
        , resComm_(&resComm)
        , pressureIndex_(pressureIndex)
        , prm_(prm)
        , weightsCalculator_(weightsCalculator)
    {
        build();
    }

    void pre(SystemVector<Scalar>&, SystemVector<Scalar>&) override
    {
        // Not forwarded to the two-level method: its pre() resizes u_ and rhs_,
        // and MultiTypeBlockVector has no resize. apply() assigns to both, which
        // sizes them, so there is nothing for pre() to do here.
    }

    void post(SystemVector<Scalar>&) override
    {
    }

    Dune::SolverCategory::Category category() const override
    {
        if constexpr (isParallel) {
            return Dune::SolverCategory::overlapping;
        } else {
            return Dune::SolverCategory::sequential;
        }
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

    void update() override
    {
        if (twoLevel_) {
            weights_ = weightsCalculator_();
            twoLevel_->updatePreconditioner(sweep_, *coarsePolicy_);
        } else {
            sweep_->update();
        }
    }

    void updateForChangedWellStructure()
    {
        // A changed well structure changes D's dimension and the dimension of
        // the coarse system, so everything is built again from the property
        // tree. That is deliberately blunt: refreshing the reservoir solves in
        // place instead would keep a CPR hierarchy rather than re-aggregate it,
        // and the two are different preconditioners -- on Norne that choice
        // alone moves general_system_cpr between 4303 and 4253 linear
        // iterations. Rebuilding is the safe default until it is measured.
        build();
    }

    void apply(SystemVector<Scalar>& v, const SystemVector<Scalar>& d) override
    {
        v[_0].resize(d[_0].size());
        v[_1].resize(d[_1].size());
        v[_0] = 0.0;
        v[_1] = 0.0;

        if (twoLevel_) {
            // TwoLevelMethodCpr accumulates into v, so it has to start at zero.
            twoLevel_->apply(v, d);
        } else {
            sweep_->apply(v, d);
        }

        // The steps leave the reservoir correction with stale overlap entries;
        // make it consistent once, as the fixed three-stage version does.
        if constexpr (isParallel) {
            resComm_->copyOwnerToAll(v[_0], v[_0]);
        }
    }

private:
    const SystemMatrix<Scalar>& S_;
    const ResComm* resComm_ = nullptr;
    int pressureIndex_ = 0;
    PropertyTree prm_;
    std::function<SystemVector<Scalar>()> weightsCalculator_;
    SystemVector<Scalar> weights_;

    std::shared_ptr<Sweep> sweep_;
    std::unique_ptr<FineOperator> fineOp_;
    std::unique_ptr<TransferPolicy> transfer_;
    std::unique_ptr<CoarseSolverPolicy> coarsePolicy_;
    std::unique_ptr<TwoLevel> twoLevel_;
    // Kept so the coarse level can be built again when the well count changes.
    std::optional<PropertyTree> coarseResPrm_;

    void build()
    {
        twoLevel_.reset();
        coarsePolicy_.reset();
        transfer_.reset();
        fineOp_.reset();
        sweep_.reset();
        coarseResPrm_.reset();

        // The weights arrive from the outer layer for the whole system; the
        // reservoir-only sub-solvers want just their own part of them.
        std::function<ResVector<Scalar>()> resWeightCalc;
        if (weightsCalculator_) {
            const auto calc = weightsCalculator_;
            resWeightCalc = [calc]() { return calc()[_0]; };
        }

        const auto coarse = prm_.get_child_optional("coarse_solver");
        const auto smoother = prm_.get_child_optional("smoother");
        if (!coarse && !smoother) {
            OPM_THROW(std::invalid_argument,
                      "general_system_cpr needs at least one of the sub-trees "
                      "'coarse_solver' and 'smoother'.");
        }

        std::vector<std::unique_ptr<SystemSweepStep<Scalar>>> steps;
        const auto addReservoir = [&](const PropertyTree& p) {
            steps.push_back(std::make_unique<ReservoirStep>(
                S_, p, resWeightCalc, pressureIndex_, resComm_));
        };
        const auto addWell = [&](const PropertyTree& p) {
            steps.push_back(std::make_unique<WellStep>(S_, p));
        };

        bool twoLevelCoarse = false;
        if (coarse) {
            if (const auto res = coarse->get_child_optional("reservoir_solver")) {
                // add_wells is the same switch the classic CPR/CPRW pair uses:
                // it promotes the pressure stage from reservoir-only CPR to
                // CPRW over the full (reservoir, well) system, and only then
                // is there a coarse space for the whole system.
                if (res->get("preconditioner.add_wells", false)) {
                    coarseResPrm_ = *res;
                    buildCoarse();
                    twoLevelCoarse = true;
                } else {
                    addReservoir(*res);
                }
            }
            if (const auto well = coarse->get_child_optional("well_solver")) {
                addWell(*well);
            }
        }
        if (smoother) {
            if (const auto res = smoother->get_child_optional("reservoir_smoother")) {
                addReservoir(*res);
            }
            if (const auto well = smoother->get_child_optional("well_solver")) {
                addWell(*well);
            }
        }

        if (steps.empty() && !twoLevelCoarse) {
            OPM_THROW(std::invalid_argument,
                      "general_system_cpr was configured with no parts at all.");
        }

        sweep_ = std::make_shared<Sweep>(std::move(steps), isParallel);

        if (twoLevelCoarse) {
            twoLevel_ = std::make_unique<TwoLevel>(*fineOp_, sweep_, *transfer_,
                                                   *coarsePolicy_,
                                                   /*preSteps=*/0, /*postSteps=*/1);
        }
    }

    void buildCoarse()
    {
        const auto& resPrm = *coarseResPrm_;
        const auto coarsePrm = resPrm.get_child_optional("preconditioner.coarsesolver")
            ? resPrm.get_child("preconditioner.coarsesolver")
            : PropertyTree();
        const auto wellTransfer = wellTransferFromString(
            prm_.get("well_transfer", std::string{"full"}));
        const auto diagonal = wellCoarseDiagonalFromString(
            prm_.get("well_coarse_diagonal", std::string{"contract_d"}));

        if (!weightsCalculator_) {
            OPM_THROW(std::invalid_argument,
                      "The CPRW pressure stage (add_wells) needs a weights calculator, but "
                      "none was configured. Set the reservoir solver's weight_type.");
        }
        weights_ = weightsCalculator_();

        fineOp_ = std::make_unique<FineOperator>(S_);
        transfer_ = std::make_unique<TransferPolicy>(S_, weights_, coarsePrm, pressureIndex_,
                                                     wellTransfer, diagonal,
                                                     prm_.get("verbosity", 0), resComm_);
        coarsePolicy_ = std::make_unique<CoarseSolverPolicy>(coarsePrm);
    }
};

} // namespace Opm

#endif // OPM_GENERALSYSTEMPRECONDITIONER_HEADER_INCLUDED
