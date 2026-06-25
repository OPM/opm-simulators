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
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>
#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>

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
#if HAVE_MPI
    using ParSysPrecondType = SystemPreconditioner<Scalar, ParResOperator<Scalar>, ParResComm>;
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
            } else
            { // Rebuild the parallel solver if the parallel preconditioner cannot be updated in-place.
                createSystemSolver(prm);
            }
            return;
        }
#endif

        if (auto* precond = dynamic_cast<SeqSysPrecondType*>(sysPrecond_)) {
            precond->updateForChangedWellStructure();
        } else
        { // Rebuild the solver if the sequential preconditioner cannot be updated in-place
            createSystemSolver(prm);
        }
    }

    void createSystemSolver(const Opm::PropertyTree& prm)
    {
        // Derive weights from the reservoir sub-block config (which uses CPR internally)
        auto resSolverPrm = prm.get_child("preconditioner.reservoir_solver");
        std::function<ResVector<Scalar>()> resWeightCalc
            = this->getWeightsCalculator(resSolverPrm, this->getMatrix(), pressureIndex);

        std::function<SystemVector<Scalar>()> sysWeightCalc;
        if (resWeightCalc) {
            sysWeightCalc = [resWeightCalc]() {
                SystemVector<Scalar> w;
                w[_0] = resWeightCalc();
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