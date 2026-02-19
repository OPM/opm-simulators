#pragma once

#include "SystemPreconditioner.hpp"
#include "WellMatrixMerger.hpp"

#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/solvers.hh>

#include <opm/simulators/linalg/ISTLSolver.hpp>

namespace Opm
{

template <class TypeTag>
class ISTLSolverExperiment : public ISTLSolver<TypeTag>
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
    ISTLSolverExperiment(const Simulator& simulator,
                         const FlowLinearSolverParameters& parameters,
                         bool forceSerial = false)
        : Parent(simulator, parameters, forceSerial)
    {
    }

    explicit ISTLSolverExperiment(const Simulator& simulator)
        : Parent(simulator)
    {
    }

    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        const auto& prm = this->prm_[this->activeSolverNum_];
        if (prm.get("use_system_solver", false)) {
            this->initPrepare(M.istlMatrix(), b);
        } else {
            Parent::prepare(M, b);
        }
    }

    void prepare(const Matrix& M, Vector& b) override
    {
        const auto& prm = this->prm_[this->activeSolverNum_];
        if (prm.get("use_system_solver", false)) {
            this->initPrepare(M, b);
        } else {
            Parent::prepare(M, b);
        }
    }

    bool solve(Vector& x) override
    {
        OPM_TIMEBLOCK(ISTLSolverExperiment_solve);

        const auto& prm = this->prm_[this->activeSolverNum_];
        const bool solve_system = prm.get("use_system_solver", false);

        if (!solve_system) {
            return Parent::solve(x);
        }

        ++this->solveCount_;

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
                           "Well" + std::to_string(i + 1));
        }
        merger.finalize();

        const size_t newWellDofs = merger.getMergedD().N();
        const bool needRebuild = !sysInitialized_ || (newWellDofs != cachedWellDofs_);

        sysMatrix_[_0][_0] = *Parent::matrix_;
        sysMatrix_[_0][_1] = merger.getMergedC();
        sysMatrix_[_1][_0] = merger.getMergedB();
        sysMatrix_[_1][_1] = merger.getMergedD();

        const auto& prm_system = prm.get_child("system_solver");

        if (needRebuild) {
            createSystemSolver(prm_system);
            sysInitialized_ = true;
            cachedWellDofs_ = newWellDofs;
        } else {
            updateSystemSolver();
        }

        const size_t numRes = Parent::matrix_->N();
        sysX_[_0].resize(numRes);
        sysX_[_0] = 0.0;
        sysX_[_1].resize(newWellDofs);
        sysX_[_1] = 0.0;

        sysRhs_[_0] = *Parent::rhs_;
        sysRhs_[_1].resize(newWellDofs);
        sysRhs_[_1] = 0.0;

        Dune::InverseOperatorResult result;
        sysSolver_->apply(sysX_, sysRhs_, result);
        this->iterations_ = result.iterations;

        x = sysX_[_0];

        return true;
    }

private:
    using SysSolverType = Dune::InverseOperator<SystemVector, SystemVector>;

    bool sysInitialized_ = false;
    size_t cachedWellDofs_ = 0;

    SystemMatrix sysMatrix_;
    SystemVector sysX_;
    SystemVector sysRhs_;

    // Serial solver components
    std::unique_ptr<Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector>> sysOp_;
    std::unique_ptr<SystemPreconditioner> sysPrecondSeq_;
    std::unique_ptr<SysSolverType> sysSolver_;

    // Parallel solver components
    using WellComm = Dune::JacComm;
    using SystemComm = Dune::MultiCommunicator<const CommunicationType&, const WellComm&>;
    std::unique_ptr<WellComm> wellComm_;
    std::unique_ptr<SystemComm> systemComm_;
    std::unique_ptr<
        Dune::OverlappingSchwarzOperator<SystemMatrix, SystemVector, SystemVector, SystemComm>>
        sysOpPar_;
    std::unique_ptr<SystemPreconditionerParallel> sysPrecondPar_;
    std::unique_ptr<Dune::BlockPreconditioner<SystemVector,
                                              SystemVector,
                                              SystemComm,
                                              SystemPreconditionerParallel>>
        sysBlockPrecond_;
    std::shared_ptr<Dune::ScalarProduct<SystemVector>> sysScalarProduct_;

    void createSystemSolver(const Opm::PropertyTree& prm_system)
    {
        std::function<ResVector()> weightCalculator
            = this->getWeightsCalculator(prm_system.get_child("preconditioner.reservoir_solver"),
                                         this->getMatrix(),
                                         pressureIndex);

        const Opm::PropertyTree precond_prm = prm_system.get_child("preconditioner");

        const bool is_parallel = this->comm_->communicator().size() > 1;
        const bool is_iorank = this->comm_->communicator().rank() == 0;
        const int verbosity = is_iorank ? prm_system.get<int>("verbosity", 0) : 0;

        if (is_parallel) {
            wellComm_ = std::make_unique<WellComm>();
            systemComm_ = std::make_unique<SystemComm>(*(this->comm_), *wellComm_);

            sysOpPar_ = std::make_unique<Dune::OverlappingSchwarzOperator<SystemMatrix,
                                                                          SystemVector,
                                                                          SystemVector,
                                                                          SystemComm>>(
                sysMatrix_, *systemComm_);

            sysScalarProduct_ = Dune::createScalarProduct<SystemVector, SystemComm>(
                *systemComm_, sysOpPar_->category());

            sysPrecondPar_ = std::make_unique<SystemPreconditionerParallel>(
                sysMatrix_, weightCalculator, pressureIndex, precond_prm, *systemComm_);

            sysBlockPrecond_
                = std::make_unique<Dune::BlockPreconditioner<SystemVector,
                                                             SystemVector,
                                                             SystemComm,
                                                             SystemPreconditionerParallel>>(
                    *sysPrecondPar_, *systemComm_);

            sysSolver_ = createKrylovSolver(
                *sysOpPar_, *sysBlockPrecond_, sysScalarProduct_.get(), prm_system, verbosity);
        } else {
            sysOp_
                = std::make_unique<Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector>>(
                    sysMatrix_);

            sysPrecondSeq_ = std::make_unique<SystemPreconditioner>(
                sysMatrix_, weightCalculator, pressureIndex, precond_prm);

            sysSolver_
                = createKrylovSolver(*sysOp_, *sysPrecondSeq_, nullptr, prm_system, verbosity);
        }
    }

    void updateSystemSolver()
    {
        const bool is_parallel = this->comm_->communicator().size() > 1;
        if (is_parallel) {
            sysPrecondPar_->update();
        } else {
            sysPrecondSeq_->update();
        }
    }

    template <class Operator, class Preconditioner>
    std::unique_ptr<SysSolverType> createKrylovSolver(Operator& op,
                                                      Preconditioner& precond,
                                                      Dune::ScalarProduct<SystemVector>* sp,
                                                      const Opm::PropertyTree& prm,
                                                      int verbosity)
    {
        const double tol = prm.get<double>("tol");
        const int maxiter = prm.get<int>("maxiter");
        const std::string solver_type = prm.get<std::string>("solver");

        if (solver_type == "bicgstab") {
            if (sp) {
                return std::make_unique<Dune::BiCGSTABSolver<SystemVector>>(
                    op, *sp, precond, tol, maxiter, verbosity);
            }
            return std::make_unique<Dune::BiCGSTABSolver<SystemVector>>(
                op, precond, tol, maxiter, verbosity);
        } else if (solver_type == "fgmres") {
            const int restart = prm.get<int>("restart", 15);
            if (sp) {
                return std::make_unique<Dune::RestartedGMResSolver<SystemVector>>(
                    op, *sp, precond, tol, restart, maxiter, verbosity);
            }
            return std::make_unique<Dune::RestartedGMResSolver<SystemVector>>(
                op, precond, tol, restart, maxiter, verbosity);
        }
        OPM_THROW(std::invalid_argument, "Properties: Solver " + solver_type + " not known.");
    }
};

} // namespace Opm
