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

        const size_t numRes = Parent::matrix_->N();
        sysX_[_0].resize(numRes);
        sysX_[_0] = 0.0;
        sysX_[_1].resize(cachedWellDofs_);
        sysX_[_1] = 0.0;

        sysRhs_[_0] = *Parent::rhs_;
        sysRhs_[_1] = mergedWellResidual_;

        Dune::InverseOperatorResult result;
        sysSolver_->apply(sysX_, sysRhs_, result);
        this->iterations_ = result.iterations;

        x = sysX_[_0];

        return this->checkConvergence(result);
    }

    std::optional<typename Parent::WellSolutionView> getWellSolution() const override
    {
        return typename Parent::WellSolutionView{sysX_[_1], wellDofOffsets_};
    }

private:
    bool sysInitialized_ = false;
    size_t cachedWellDofs_ = 0;
    WellVector mergedWellResidual_;
    std::vector<int> wellDofOffsets_;

    // Owned storage for merged well matrices; SystemMatrix points into these.
    WRMatrix mergedB_;
    RWMatrix mergedC_;
    WWMatrix mergedD_;

    SystemMatrix sysMatrix_;
    SystemVector sysX_;
    SystemVector sysRhs_;

    // Serial solver components
    std::unique_ptr<SystemSeqOp> sysOp_;
    std::unique_ptr<Dune::FlexibleSolver<SystemSeqOp>> sysFlexSolverSeq_;

    // Parallel solver components
#if HAVE_MPI
    using WellComm = Dune::JacComm;
    std::unique_ptr<WellComm> wellComm_;
    std::unique_ptr<SystemComm> systemComm_;
    std::unique_ptr<SystemParOp> sysOpPar_;
    std::unique_ptr<Dune::FlexibleSolver<SystemParOp>> sysFlexSolverPar_;
#endif

    using SysSolverType = Dune::InverseOperator<SystemVector, SystemVector>;
    using SysPrecondType = Dune::PreconditionerWithUpdate<SystemVector, SystemVector>;
    SysSolverType* sysSolver_ = nullptr;
    SysPrecondType* sysPrecond_ = nullptr;

    void prepareSystemSolver()
    {
        OPM_TIMEBLOCK(flexibleSolverPrepare);

        std::vector<WRMatrix> b_matrices;
        std::vector<RWMatrix> c_matrices;
        std::vector<WWMatrix> d_matrices;
        std::vector<std::vector<int>> wcells;
        std::vector<WellVector> well_residuals;

        this->simulator_.problem().wellModel().addBCDMatrix(
            b_matrices, c_matrices, d_matrices, wcells, well_residuals);

        Opm::WellMatrixMerger merger(Parent::matrix_->N());
        for (size_t i = 0; i < b_matrices.size(); ++i) {
            merger.addWell(b_matrices[i],
                           c_matrices[i],
                           d_matrices[i],
                           wcells[i],
                           static_cast<int>(i),
                           "Well" + std::to_string(i + 1),
                           well_residuals[i]);
        }
        merger.finalize();

        mergedWellResidual_ = std::move(merger.getMergedWellResidual());
        wellDofOffsets_ = std::move(merger.getWellDofOffsets());

        const size_t newWellDofs = merger.getMergedD().N();
        const bool localNeedRebuild = !sysInitialized_ || (newWellDofs != cachedWellDofs_);

        // All ranks must agree: rebuild if ANY rank needs it, since
        // create and update take different MPI-collective code paths
        // (AMG hierarchy construction vs. update).
        const bool needRebuild
            = this->comm_->communicator().max(static_cast<int>(localNeedRebuild)) > 0;

        mergedB_ = std::move(merger.getMergedB());
        mergedC_ = std::move(merger.getMergedC());
        mergedD_ = std::move(merger.getMergedD());
        sysMatrix_.A = Parent::matrix_;
        sysMatrix_.B = &mergedB_;
        sysMatrix_.C = &mergedC_;
        sysMatrix_.D = &mergedD_;
        cachedWellDofs_ = newWellDofs;

        const auto& prm = this->prm_[this->activeSolverNum_];

        if (needRebuild) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            createSystemSolver(prm);
            sysInitialized_ = true;
        } else {
            OPM_TIMEBLOCK(flexibleSolverUpdate);
            sysPrecond_->update();
        }
    }

    void createSystemSolver(const Opm::PropertyTree& prm)
    {
        // Derive weights from the reservoir sub-block config (which uses CPR internally)
        auto resSolverPrm = prm.get_child("preconditioner.reservoir_solver");
        std::function<ResVector()> resWeightCalc
            = this->getWeightsCalculator(resSolverPrm, this->getMatrix(), pressureIndex);

        std::function<SystemVector()> sysWeightCalc;
        if (resWeightCalc) {
            sysWeightCalc = [resWeightCalc]() {
                SystemVector w;
                w[_0] = resWeightCalc();
                return w;
            };
        }

        const bool is_parallel = this->comm_->communicator().size() > 1;

        if (is_parallel) {
#if HAVE_MPI
            wellComm_ = std::make_unique<WellComm>();
            systemComm_ = std::make_unique<SystemComm>(*(this->comm_), *wellComm_);

            sysOpPar_ = std::make_unique<SystemParOp>(sysMatrix_, *systemComm_);

            sysFlexSolverPar_ = std::make_unique<Dune::FlexibleSolver<SystemParOp>>(
                *sysOpPar_, *systemComm_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverPar_.get();
            sysPrecond_ = &sysFlexSolverPar_->preconditioner();
#endif
        } else {
            sysOp_ = std::make_unique<SystemSeqOp>(sysMatrix_);

            sysFlexSolverSeq_ = std::make_unique<Dune::FlexibleSolver<SystemSeqOp>>(
                *sysOp_, prm, sysWeightCalc, pressureIndex);

            sysSolver_ = sysFlexSolverSeq_.get();
            sysPrecond_ = &sysFlexSolverSeq_->preconditioner();
        }
    }
};

} // namespace Opm
