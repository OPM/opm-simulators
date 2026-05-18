#pragma once

#include "SystemPreconditionerFactory.hpp"
#include "WellMatrixMerger.hpp"

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>

namespace Opm
{

template <class TypeTag>
class ISTLSolverSystem : public ISTLSolver<TypeTag>
{
protected:
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using WellModel = GetPropType<TypeTag, Properties::WellModel>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;
    using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using ElementChunksType = ElementChunks<GridView, Dune::Partitions::All>;

    constexpr static std::size_t pressureIndex
        = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

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

        const std::size_t numRes = Parent::matrix_->N();
        const std::size_t numWell = cachedWellStructure_.totalWellDofs;

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

    // Current per-well B/C/D blocks for the explicit 2x2 system matrix.
    std::vector<WRMatrixT<Scalar>> wellBMatrices_;
    std::vector<RWMatrixT<Scalar>> wellCMatrices_;
    std::vector<WWMatrixT<Scalar>> wellDMatrices_;
    std::vector<std::vector<int>> wellCells_;

    // Owned storage for merged well matrices; SystemMatrix points into these.
    WRMatrixT<Scalar> mergedB_;
    RWMatrixT<Scalar> mergedC_;
    WWMatrixT<Scalar> mergedD_;

    SystemMatrixT<Scalar> sysMatrix_;
    SystemVectorT<Scalar> sysX_;
    SystemVectorT<Scalar> sysRhs_;

    // Serial solver components
    std::unique_ptr<SystemSeqOpT<Scalar>> sysOp_;
    std::unique_ptr<Dune::FlexibleSolver<SystemSeqOpT<Scalar>>> sysFlexSolverSeq_;

    // Parallel solver components
#if HAVE_MPI
    using WellComm = Dune::JacComm;
    std::unique_ptr<WellComm> wellComm_;
    std::unique_ptr<SystemComm> systemComm_;
    std::unique_ptr<SystemParOpT<Scalar>> sysOpPar_;
    std::unique_ptr<Dune::FlexibleSolver<SystemParOpT<Scalar>>> sysFlexSolverPar_;
#endif

    using SysSolverType = Dune::InverseOperator<SystemVectorT<Scalar>, SystemVectorT<Scalar>>;
    using SysPrecondType = Dune::PreconditionerWithUpdate<SystemVectorT<Scalar>, SystemVectorT<Scalar>>;
    using SeqSysPrecondType = SystemPreconditioner<Scalar, SeqResOperatorT<Scalar>>;
#if HAVE_MPI
    using ParSysPrecondType = SystemPreconditioner<Scalar, ParResOperatorT<Scalar>, ParResComm>;
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

        const Opm::WellMatrixMerger<Scalar> merger(
            Parent::matrix_->N(), wellBMatrices_, wellCMatrices_, wellDMatrices_, wellCells_);

        const bool localStructureChanged = !sysInitialized_
            || !merger.hasSameStructure(cachedWellStructure_);

        // All ranks must agree on whether to take the structure-change path,
        // because the distributed solver create and update paths use different
        // MPI-collective sequences.
        const bool needStructureRefresh = !sysInitialized_
            || this->comm_->communicator().max(static_cast<int>(localStructureChanged)) > 0;

        const auto& prm = this->prm_[this->activeSolverNum_];

        if (needStructureRefresh) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            merger.buildMatrices(mergedB_, mergedC_, mergedD_);
            sysMatrix_.A = Parent::matrix_;
            sysMatrix_.B = &mergedB_;
            sysMatrix_.C = &mergedC_;
            sysMatrix_.D = &mergedD_;
            cachedWellStructure_ = merger.buildStructure();

            refreshSystemSolverForChangedWellStructure(prm);
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
            sysPrecond_->update();
        }
    }

    void refreshSystemSolverForChangedWellStructure(const Opm::PropertyTree& prm)
    {
        if (!sysInitialized_ || !sysPrecond_) {
            createSystemSolver(prm);
            return;
        }

#if HAVE_MPI
        if (this->comm_->communicator().size() > 1) {
            if (auto* precond = dynamic_cast<ParSysPrecondType*>(sysPrecond_)) {
                precond->updateForChangedWellStructure();
                return;
            }
            createSystemSolver(prm);
            return;
        }
#endif

        if (auto* precond = dynamic_cast<SeqSysPrecondType*>(sysPrecond_)) {
            precond->updateForChangedWellStructure();
            return;
        }

        createSystemSolver(prm);
    }

    void createSystemSolver(const Opm::PropertyTree& prm)
    {
        // Derive weights from the reservoir sub-block config (which uses CPR internally)
        auto resSolverPrm = prm.get_child("preconditioner.reservoir_solver");
        std::function<ResVectorT<Scalar>()> resWeightCalc
            = this->getWeightsCalculator(resSolverPrm, this->getMatrix(), pressureIndex);

        std::function<SystemVectorT<Scalar>()> sysWeightCalc;
        if (resWeightCalc) {
            sysWeightCalc = [resWeightCalc]() {
                SystemVectorT<Scalar> w;
                w[_0] = resWeightCalc();
                return w;
            };
        }

        const bool is_parallel = this->comm_->communicator().size() > 1;

        if (is_parallel) {
#if HAVE_MPI
            wellComm_ = std::make_unique<WellComm>();
            systemComm_ = std::make_unique<SystemComm>(*(this->comm_), *wellComm_);

            sysOpPar_ = std::make_unique<SystemParOpT<Scalar>>(sysMatrix_, *systemComm_);

            sysFlexSolverPar_ = std::make_unique<Dune::FlexibleSolver<SystemParOpT<Scalar>>>(
                *sysOpPar_, *systemComm_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverPar_.get();
            sysPrecond_ = &sysFlexSolverPar_->preconditioner();
#endif
        } else {
            sysOp_ = std::make_unique<SystemSeqOpT<Scalar>>(sysMatrix_);

            sysFlexSolverSeq_ = std::make_unique<Dune::FlexibleSolver<SystemSeqOpT<Scalar>>>(
                *sysOp_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverSeq_.get();
            sysPrecond_ = &sysFlexSolverSeq_->preconditioner();
        }
    }
};

} // namespace Opm
