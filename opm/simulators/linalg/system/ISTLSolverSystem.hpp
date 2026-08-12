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
#ifndef OPM_ISTLSOLVERSYSTEM_HEADER_INCLUDED
#define OPM_ISTLSOLVERSYSTEM_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/system/GeneralSystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>
#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm
{

template <class TypeTag>
class ISTLSolverSystem : public ISTLSolver<TypeTag>
{
protected:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    // Compile-time validation: SystemPreconditionerFactory and related types
    // are hardcoded for standard 3-phase blackoil (3 reservoir equations, 4 well equations).
    // See SystemTypes.hpp for details.
    static_assert(Indices::numEq == 3,
                  "ISTLSolverSystem (with system_cpr preconditioner) only supports "
                  "3-equation blackoil models. This model has different equation count.");

    constexpr static std::size_t pressureIndex
        = Indices::pressureSwitchIdx;

    enum { enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>() };
    constexpr static bool isIncompatibleWithCprw = enablePolymerMolarWeight;

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif
    using Parent = ISTLSolver<TypeTag>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

public:
    ISTLSolverSystem(const Simulator& simulator,
                     const FlowLinearSolverParameters& parameters,
                     bool forceSerial = false)
        : Parent(simulator, parameters, forceSerial)
    {
    }

    explicit ISTLSolverSystem(const Simulator& simulator)
        : Parent(simulator)
    {
    }

    ~ISTLSolverSystem()
    {
        if (thpEliminationTried_ > 0) {
            OpmLog::debug(fmt::format("cellavg_vfp: control-row elimination "
                                      "applied {} of {} times",
                                      thpEliminationTried_ - thpEliminationRejected_,
                                      thpEliminationTried_));
        }
    }

    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        OPM_TIMEBLOCK(istlSolverPrepare);
        this->initPrepare(M.istlMatrix(), b);
        prepareSystemSolver();
    }

    void prepare(const Matrix& M, Vector& b) override
    {
        OPM_TIMEBLOCK(istlSolverPrepare);
        this->initPrepare(M, b);
        prepareSystemSolver();
    }

    bool solve(Vector& x) override
    {
        OPM_TIMEBLOCK(istlSolverSolve);
        ++this->solveCount_;

        // Same fine-system dump as ISTLSolver::solve(), which this overrides,
        // so the reservoir matrix and rhs can be diffed against the classic path.
        if (this->prm_[this->activeSolverNum_].get("verbosity", 0) > 10) {
            Helper::writeSystem(this->simulator_,
                                this->getMatrix(),
                                *Parent::rhs_,
                                this->comm_.get());
        }

        const std::size_t numRes = Parent::matrix_->N();
        const std::size_t numWell = cachedWellStructure_.totalWellBlocks;

        sysX_[_0].resize(numRes);
        sysX_[_0] = 0.0;
        sysX_[_1].resize(numWell);
        sysX_[_1] = 0.0;

        sysRhs_[_0].resize(numRes);
        sysRhs_[_0] = *Parent::rhs_;
        sysRhs_[_1].resize(numWell);
        sysRhs_[_1] = 0.0;

        Dune::InverseOperatorResult result;
        sysSolver_->apply(sysX_, sysRhs_, result);
        this->iterations_ = result.iterations;

        x = sysX_[_0];

        return this->checkConvergence(result);
    }

private:
    bool sysInitialized_ = false;
    WellMatrixStructure cachedWellStructure_;

    // Aggregation of merged well block rows into wells, and the well weights
    // used by the CPRW pressure stage.  Both are produced here, in the outer
    // layer, from data already extracted from the well model; the
    // preconditioner consumes them as plain numbers.
    WellDofLayout wellLayout_;
    std::string wellWeightType_ = "quasiimpes";
    std::string wellThpTreatment_ = "bhp";
    std::vector<char> thpControlled_;  // per well, order as wellDMatrices_
    // Diagnostics for cellavg_vfp: how often the control-row elimination was
    // attempted and how often the diagonal guard refused it.
    mutable long thpEliminationTried_ = 0;
    mutable long thpEliminationRejected_ = 0;

    // Current per-well B/C/D blocks for the explicit 2x2 system matrix.
    std::vector<WRMatrix<Scalar>> wellBMatrices_;
    std::vector<RWMatrix<Scalar>> wellCMatrices_;
    std::vector<WWMatrix<Scalar>> wellDMatrices_;
    Opm::SparseTable<int> wellCells_;

    // Owned storage for merged well matrices; SystemMatrix points into these.
    WRMatrix<Scalar> mergedB_;
    RWMatrix<Scalar> mergedC_;
    WWMatrix<Scalar> mergedD_;

    SystemMatrix<Scalar> sysMatrix_;
    SystemVector<Scalar> sysX_;
    SystemVector<Scalar> sysRhs_;

    // Serial solver components
    std::unique_ptr<SystemSeqOp<Scalar>> sysOp_;
    std::unique_ptr<Dune::FlexibleSolver<SystemSeqOp<Scalar>>> sysFlexSolverSeq_;

    // Parallel solver components
#if HAVE_MPI
    using WellComm = Dune::JacComm;
    std::unique_ptr<WellComm> wellComm_;
    std::unique_ptr<SystemComm> systemComm_;
    std::unique_ptr<SystemParOp<Scalar>> sysOpPar_;
    std::unique_ptr<Dune::FlexibleSolver<SystemParOp<Scalar>>> sysFlexSolverPar_;
#endif

    using SysSolverType = Dune::InverseOperator<SystemVector<Scalar>, SystemVector<Scalar>>;
    using SysPrecondType = Dune::PreconditionerWithUpdate<SystemVector<Scalar>, SystemVector<Scalar>>;
    using SeqSysPrecondType = SystemPreconditioner<Scalar, SeqResOperator<Scalar>>;
    using SeqGeneralSysPrecondType = GeneralSystemPreconditioner<Scalar, SeqResOperator<Scalar>>;
#if HAVE_MPI
    using ParSysPrecondType = SystemPreconditioner<Scalar, ParResOperator<Scalar>, ParResComm>;
    using ParGeneralSysPrecondType = GeneralSystemPreconditioner<Scalar, ParResOperator<Scalar>, ParResComm>;
#endif
    SysSolverType* sysSolver_ = nullptr;
    SysPrecondType* sysPrecond_ = nullptr;

    void prepareSystemSolver()
    {
        OPM_TIMEBLOCK(flexibleSolverPrepare);

        wellBMatrices_.clear();
        wellCMatrices_.clear();
        wellDMatrices_.clear();
        wellCells_.clear();

        this->simulator_.problem().wellModel().addBCDMatrix(
            wellBMatrices_, wellCMatrices_, wellDMatrices_, wellCells_);

        {
            // Read the coarse-space options before building the layout, so
            // the very first assembly already honours them.
            const auto wellOpts = coarseSpaceTree(this->prm_[this->activeSolverNum_]);
            wellWeightType_ = wellOpts.get("well_weight_type", std::string{"cellavg"});
            // Give a pressure-controlled well a trivial coarse equation, as the
            // classic CPRW does.  Off keeps the contracted equation for every well.
            wellLayout_.identityOnPressureControl
                = wellOpts.get("well_identity_on_pressure_control", false);
            // How a THP-controlled well enters the coarse system.  The
            // trivial row assumes constant bhp, the contracted row constant
            // rates; under THP neither holds, and "vfp" eliminates the rate
            // through the control equation instead.
            wellThpTreatment_
                = wellOpts.get("well_thp_treatment", std::string{"bhp"});
            if (wellThpTreatment_ != "bhp" && wellThpTreatment_ != "rate"
                && wellThpTreatment_ != "vfp") {
                OPM_THROW(std::invalid_argument,
                          "well_thp_treatment must be bhp, rate or vfp, got "
                          + wellThpTreatment_);
            }
        }

        buildWellDofLayout();

        const Opm::WellMatrixMerger<Scalar> merger(
            Parent::matrix_->N(), wellBMatrices_, wellCMatrices_, wellDMatrices_, wellCells_);

        const bool localStructureChanged = !sysInitialized_
            || !merger.hasSameStructure(cachedWellStructure_);

        // All ranks must agree on whether to take the structure-change path,
        // because the distributed solver create and update paths use different
        // MPI-collective sequences.
#if HAVE_MPI
        const bool globalStructureChanged = this->comm_->communicator().max(
            static_cast<int>(localStructureChanged)) > 0;
#else
        const bool globalStructureChanged = localStructureChanged;
#endif
        const bool needStructureRefresh = !sysInitialized_ || globalStructureChanged;

        const auto& prm = this->prm_[this->activeSolverNum_];

        if (needStructureRefresh) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            merger.buildMatrices(mergedB_, mergedC_, mergedD_);
            sysMatrix_.A = Parent::matrix_;
            sysMatrix_.B = &mergedB_;
            sysMatrix_.C = &mergedC_;
            sysMatrix_.D = &mergedD_;
            sysMatrix_.wellLayout = &wellLayout_;

            const auto newStructure = merger.buildStructure();
            // A connection opening or closing inside an existing well keeps
            // every dimension; a well or segment appearing or vanishing does
            // not, and only the latter introduces unknowns the initial build
            // never saw.
            const auto change = (sysInitialized_ && newStructure.hasSameDimensions(cachedWellStructure_))
                ? WellStructureChange::Pattern
                : WellStructureChange::Dimension;
            cachedWellStructure_ = newStructure;

            refreshSystemSolverForChangedWellStructure(prm, change);
            sysInitialized_ = true;
        } else {
            OPM_TIMEBLOCK(flexibleSolverUpdate);
            // Pattern unchanged: write fresh values into the existing merged
            // matrices without any (de)allocation.
            merger.updateValues(mergedB_, mergedC_, mergedD_);

            // Refresh A pointer in case the reservoir matrix was reallocated.
            sysMatrix_.A = Parent::matrix_;
            sysMatrix_.B = &mergedB_;
            sysMatrix_.C = &mergedC_;
            sysMatrix_.D = &mergedD_;
            sysMatrix_.wellLayout = &wellLayout_;
            sysPrecond_->update();
        }
    }

    // Which merged well block rows belong to which well.  The merged D matrix
    // is the per-well D blocks concatenated, so this is a plain prefix sum
    // over their dimensions: one block for a standard well, one per segment
    // for a multisegment well.
    void buildWellDofLayout()
    {
        auto& offsets = wellLayout_.wellBlockOffsets;
        offsets.clear();
        offsets.reserve(wellDMatrices_.size() + 1);
        offsets.push_back(0);
        std::size_t total = 0;
        for (const auto& d : wellDMatrices_) {
            total += d.N();
            offsets.push_back(total);
        }

        // Which wells are on pressure control.  Asking this is the outer
        // layer's job; below here it is just a flag per well.  The order
        // matches addBCDMatrix, which walks the same well container.
        wellLayout_.pressureControlled.clear();
        thpControlled_.clear();
        {
            const auto& wellModel = this->simulator_.problem().wellModel();
            const auto& wellState = wellModel.wellState();
            thpControlled_.reserve(wellDMatrices_.size());
            wellLayout_.pressureControlled.reserve(wellDMatrices_.size());
            std::size_t j = 0;
            for (const auto& well : wellModel) {
                const auto& ws = wellState.well(well->indexOfWell());
                const bool thp = well->isInjector()
                    ? ws.injection_cmode == Well::InjectorCMode::THP
                    : ws.production_cmode == Well::ProducerCMode::THP;
                thpControlled_.push_back(thp ? 1 : 0);

                bool pc = well->isPressureControlled(wellState);
                if (pc && thp) {
                    if (wellThpTreatment_ == "rate") {
                        // A THP well keeps its contracted row; the trivial
                        // row is reserved for BHP control, where it is exact.
                        pc = false;
                    } else if (wellThpTreatment_ == "vfp") {
                        // The control-row elimination replaces the trivial row
                        // for a single-block well with a usable VFP derivative;
                        // a flat VFP behaves like BHP control, so keep the
                        // trivial row there (and for multisegment wells).
                        const bool singleBlock = wellDMatrices_[j].N() == 1;
                        pc = !(singleBlock
                               && wellDMatrices_[j].exists(0, 0)
                               && usableThpDerivative(wellDMatrices_[j][0][0]));
                    }
                    // "bhp": THP counts as pressure control, as classic cprw.
                }
                if (wellLayout_.identityOnPressureControl) {
                    wellLayout_.pressureControlled.push_back(pc ? 1 : 0);
                }
                ++j;
            }
        }
    }

    // The control-row elimination divides by the control equation's
    // total-rate derivative; refuse it when that entry is negligible
    // against the rest of the row (flat VFP, or BHP control).
    template <class Block>
    static bool usableThpDerivative(const Block& d)
    {
        const int k = numWellDofs - 1;
        Scalar rowMax = 0.0;
        for (int c = 0; c < numWellDofs; ++c) {
            rowMax = std::max(rowMax, std::abs(d[k][c]));
        }
        return rowMax > 0.0 && std::abs(d[k][0]) > 1e-8 * rowMax;
    }

    // Weights used to contract each well's equations down to the single scalar
    // the CPRW pressure system carries for that well.  Computed here rather
    // than inside the preconditioner so that the linear-solver core never sees
    // anything well-specific, and so that this can later be replaced by a
    // value obtained from the well model without touching the core.
    // Merged well block row -> well index.
    std::size_t wellOfBlock(const std::size_t blockRow) const
    {
        const auto& off = wellLayout_.wellBlockOffsets;
        const auto it = std::upper_bound(off.begin(), off.end(), blockRow);
        assert(it != off.begin() && it != off.end());
        return static_cast<std::size_t>(std::distance(off.begin(), it) - 1);
    }

    WellVector<Scalar> computeWellWeights(const ResVector<Scalar>& resWeights) const
    {
        const std::size_t numBlocks = mergedD_.N();
        const int q = wellLayout_.pressureDofIndex;

        WellVector<Scalar> weights(numBlocks);
        for (std::size_t wb = 0; wb < numBlocks; ++wb) {
            auto& lambda = weights[wb];
            lambda = 0.0;

            if (wellWeightType_ == "unit") {
                // Pick the pressure row of the well equations as-is.
                lambda[q] = 1.0;
                continue;
            }

            if (wellWeightType_ == "cellavg" || wellWeightType_ == "cellblockavg") {
                // The classic CPRW weighting (use_well_weights = false):
                // average the reservoir weights over perforated cells and use
                // them on the conservation equations only, weight zero on the
                // control equation.
                //
                // "cellavg" averages over every perforation of the whole well
                // and gives every block of that well the same weights, which is
                // what MultisegmentWellEquations::extractCPRPressureMatrix
                // does. "cellblockavg" averages per block row instead, which
                // is a finer but non-classic variant.
                const bool perWell = (wellWeightType_ == "cellavg");
                const std::size_t first = perWell ? wellLayout_.firstBlock(wellOfBlock(wb)) : wb;
                const std::size_t last = perWell ? wellLayout_.endBlock(wellOfBlock(wb)) : wb + 1;
                int nperf = 0;
                for (std::size_t b = first; b < last; ++b) {
                    for (auto col = mergedB_[b].begin(), end = mergedB_[b].end(); col != end; ++col) {
                        const auto& cw = resWeights[col.index()];
                        for (int i = 0; i < numResDofs; ++i) {
                            lambda[i] += cw[i];
                        }
                        ++nperf;
                    }
                }
                if (nperf > 0) {
                    for (int i = 0; i < numResDofs; ++i) {
                        lambda[i] /= nperf;
                    }
                } else {
                    // No perforations of this well on this rank; regularise
                    // rather than leaving an empty row.
                    for (int i = 0; i < numResDofs; ++i) {
                        lambda[i] = 1.0;
                    }
                }
                lambda[q] = 0.0;

                if (wellThpTreatment_ == "vfp") {
                    // Eliminate the total-rate correction through the control
                    // equation: weight mu on the control row cancels the
                    // WQTotal column of the contracted equation, which turns
                    // the coarse diagonal into the Schur value
                    //   w'D_cb - (w'D_cq) D_kb / D_kq.
                    // Rate control has D_kb = 0, so only THP wells change.
                    // Acts through the D contraction only, so it is a no-op
                    // under well_coarse_diagonal = row_sum.
                    const auto j = wellOfBlock(wb);
                    const bool singleBlock
                        = wellLayout_.endBlock(j) - wellLayout_.firstBlock(j) == 1;
                    if (singleBlock && j < thpControlled_.size() && thpControlled_[j]
                        && mergedD_.exists(wb, wb)) {
                        const auto& d = mergedD_[wb][wb];
                        if (usableThpDerivative(d)) {
                            const int k = numWellDofs - 1;
                            Scalar num = 0.0;
                            Scalar diag0 = 0.0;
                            for (int i = 0; i < k; ++i) {
                                num += lambda[i] * d[i][0];
                                diag0 += lambda[i] * d[i][k];
                            }
                            const Scalar mu = -num / d[k][0];
                            const Scalar diagc = diag0 + mu * d[k][k];
                            // On the unstable VFP branch the Schur value can
                            // shrink or flip the diagonal; keep the plain
                            // contracted row there.
                            if (diagc * diag0 > 0.01 * diag0 * diag0) {
                                lambda[k] = mu;
                            } else {
                                ++thpEliminationRejected_;
                            }
                            ++thpEliminationTried_;
                        }
                    }
                }
                continue;
            }

            // Quasi-IMPES well weights: lambda = D_ii^-T e_q, scaled to unit
            // max norm.  This is the analogue of the use_well_weights=true
            // branch of StandardWellEquations::extractCPRPressureMatrix, and
            // it needs no knowledge of the well's control mode.
            Dune::FieldVector<Scalar, numWellDofs> rhs(0.0);
            rhs[q] = 1.0;
            bool ok = false;
            if (mergedD_.exists(wb, wb)) {
                try {
                    const auto dt = mergedD_[wb][wb].transposed();
                    dt.solve(lambda, rhs);
                    Scalar absMax = 0.0;
                    for (int i = 0; i < numWellDofs; ++i) {
                        absMax = std::max(absMax, std::abs(lambda[i]));
                    }
                    if (absMax > 0.0 && std::isfinite(absMax)) {
                        lambda /= absMax;
                        ok = true;
                    }
                } catch (const Dune::FMatrixError&) {
                    ok = false;
                }
            }
            if (!ok) {
                // Singular or degenerate well block: fall back to the plain
                // pressure row rather than poisoning the coarse system.
                lambda = 0.0;
                lambda[q] = 1.0;
            }
        }
        return weights;
    }

    void refreshSystemSolverForChangedWellStructure(const Opm::PropertyTree& prm,
                                                   const WellStructureChange change)
    {
        if (!sysInitialized_ || !sysPrecond_) {
            createSystemSolver(prm);
            return;
        }

#if HAVE_MPI
        if (this->comm_->communicator().size() > 1) {
            if (auto* precond = dynamic_cast<ParSysPrecondType*>(sysPrecond_)) {
                precond->updateForChangedWellStructure();
            } else if (auto* general = dynamic_cast<ParGeneralSysPrecondType*>(sysPrecond_)) {
                general->updateForChangedWellStructure(change);
            } else
            { // Rebuild the parallel solver if the parallel preconditioner cannot be updated in-place.
                createSystemSolver(prm);
            }
            return;
        }
#endif

        if (auto* precond = dynamic_cast<SeqSysPrecondType*>(sysPrecond_)) {
            precond->updateForChangedWellStructure();
        } else if (auto* general = dynamic_cast<SeqGeneralSysPrecondType*>(sysPrecond_)) {
            general->updateForChangedWellStructure(change);
        } else
        { // Rebuild the solver if the sequential preconditioner cannot be updated in-place
            createSystemSolver(prm);
        }
    }

    // The keys describing the coarse space and its weighting.  general_system_cpr
    // keeps them beside the coarse solver; system_cpr keeps them at the top of
    // the preconditioner and its weighting inside the reservoir solver.
    Opm::PropertyTree coarseSpaceTree(const Opm::PropertyTree& prm) const
    {
        if (auto general = prm.get_child_optional("preconditioner.coarsesolver")) {
            return *general;
        }
        return prm.get_child("preconditioner");
    }

    void createSystemSolver(const Opm::PropertyTree& prm)
    {
        std::function<ResVector<Scalar>()> resWeightCalc;
        if (prm.get("preconditioner.type", std::string{}) == "general_system_cpr") {
            // The general layout names the weighting directly, at the top of
            // the preconditioner as cpr does, rather than through a nested CPR
            // preconditioner sub-tree.
            resWeightCalc = this->makeWeightsCalculator(
                prm.get("preconditioner.weight_type", std::string{"trueimpes"}),
                this->getMatrix(), pressureIndex);
        } else {
            // Derive weights from the reservoir sub-block config (which uses CPR internally)
            resWeightCalc = this->getWeightsCalculator(
                prm.get_child("preconditioner.reservoir_solver"), this->getMatrix(), pressureIndex);
        }

        // The well part of the weights is filled here too: the CPRW pressure
        // stage restricts the well rows with it, and re-reads it on every
        // update, so it has to track the current merged D.
        std::function<SystemVector<Scalar>()> sysWeightCalc;
        if (resWeightCalc) {
            sysWeightCalc = [this, resWeightCalc]() {
                SystemVector<Scalar> w;
                w[_0] = resWeightCalc();
                w[_1] = this->computeWellWeights(w[_0]);
                return w;
            };
        }

#if HAVE_MPI
        const bool is_parallel = this->comm_->communicator().size() > 1;
        if (is_parallel) {
            wellComm_ = std::make_unique<WellComm>();
            systemComm_ = std::make_unique<SystemComm>(*(this->comm_), *wellComm_);

            sysOpPar_ = std::make_unique<SystemParOp<Scalar>>(sysMatrix_, *systemComm_);

            sysFlexSolverPar_ = std::make_unique<Dune::FlexibleSolver<SystemParOp<Scalar>>>(
                *sysOpPar_, *systemComm_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverPar_.get();
            sysPrecond_ = &sysFlexSolverPar_->preconditioner();
        }
        else
#endif
        {
            sysOp_ = std::make_unique<SystemSeqOp<Scalar>>(sysMatrix_);

            sysFlexSolverSeq_ = std::make_unique<Dune::FlexibleSolver<SystemSeqOp<Scalar>>>(
                *sysOp_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverSeq_.get();
            sysPrecond_ = &sysFlexSolverSeq_->preconditioner();
        }
    }
};

} // namespace Opm

#endif // OPM_ISTLSOLVERSYSTEM_HEADER_INCLUDED
