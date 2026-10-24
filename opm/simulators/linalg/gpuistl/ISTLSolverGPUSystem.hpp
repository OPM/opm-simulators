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
#ifndef OPM_ISTLSOLVERGPUSYSTEM_HEADER_INCLUDED
#define OPM_ISTLSOLVERGPUSYSTEM_HEADER_INCLUDED

#include <opm/simulators/linalg/gpuistl/GpuSystemBackend.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSystemPreconditionerFactory.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSystemTypes.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <opm/simulators/linalg/AbstractISTLSolver.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/printlinearsolverparameter.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/common/CriticalError.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <dune/common/parallel/communication.hh>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::gpuistl
{

// GPU counterpart of ISTLSolverSystem: the full system solver with System
// CPR preconditioner, with the reservoir sub-solvers on the GPU.  Selected with --linear-solver=system_cpr (or a system_cpr JSON
// file) and --linear-solver-accelerator=gpu.
//
// Sequential runs only.  The CPRW pressure stage (system_cprw, add_wells) is
// not supported.
template <class TypeTag>
class ISTLSolverGPUSystem : public AbstractISTLSolver<GetPropType<TypeTag, Properties::SparseMatrixAdapter>,
                                                      GetPropType<TypeTag, Properties::GlobalEqVector>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ElementChunksType = Opm::ElementChunks<GridView, Dune::Partitions::All>;
    using Parent = AbstractISTLSolver<SparseMatrixAdapter, Vector>;

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    using GpuVectorType = GpuVector<Scalar>;
    using SysOp = GpuSystemSeqOp<Scalar>;
    using SysSolver = Dune::FlexibleSolver<SysOp>;
    using SysPrecond = GpuSystemPreconditioner<Scalar>;

    static_assert(Indices::numEq == numResDofs,
                  "ISTLSolverGPUSystem (with system_cpr preconditioner) only supports "
                  "3-equation blackoil models. This model has different equation count.");

    constexpr static std::size_t pressureIndex = Indices::pressureSwitchIdx;

public:
    ISTLSolverGPUSystem(const Simulator& simulator,
                        const FlowLinearSolverParameters& parameters,
                        bool forceSerial = false)
        : simulator_(simulator)
        , parameters_(parameters)
        , elementChunks_(simulator.vanguard().gridView(), Dune::Partitions::all, ThreadManager::maxThreads())
    {
#if HAVE_MPI
        comm_ = std::make_shared<CommunicationType>(simulator.vanguard().grid().comm());
#else
        comm_ = std::make_shared<CommunicationType>(simulator.gridView().comm());
#endif
        if (!forceSerial && simulator.gridView().comm().size() > 1) {
            OPM_THROW(std::invalid_argument,
                      "The system_cpr preconditioner with --linear-solver-accelerator=gpu "
                      "does not support parallel runs yet.");
        }

        parameters_.init(simulator.vanguard().eclState().getSimulationConfig().useCPR());
        prm_ = setupPropertyTree(parameters_,
                                 Parameters::IsSet<Parameters::LinearSolverMaxIter>(),
                                 Parameters::IsSet<Parameters::LinearSolverReduction>());

        Opm::detail::printLinearSolverParameters(parameters_, prm_, simulator.gridView().comm());
    }

    explicit ISTLSolverGPUSystem(const Simulator& simulator)
        : ISTLSolverGPUSystem(simulator, FlowLinearSolverParameters(), false)
    {
    }

    void eraseMatrix() override
    {
    }

    void setActiveSolver(int num) override
    {
        if (num != 0) {
            OPM_THROW(std::logic_error, "Only one solver available for the GPU system solver.");
        }
    }

    int numAvailableSolvers() const override
    {
        return 1;
    }

    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        prepare(M.istlMatrix(), b);
    }

    void prepare(const Matrix& M, Vector& b) override
    {
        OPM_TIMEBLOCK(istlSolverGpuSystemPrepare);
        try {
            prepareSystemSolver(M, b);
        }
        OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR("This is likely due to a faulty linear solver JSON specification. "
                                                "Check for errors related to missing nodes.");
    }

    void setResidual(Vector&) override
    {
        // Handled in prepare().
    }

    void getResidual(Vector& b) const override
    {
        if (!sysRhs_) {
            OPM_THROW(std::runtime_error, "The right-hand side is not set; prepare(matrix, rhs) needs to be called");
        }
        sysRhs_->res.copyToHost(b);
    }

    void setMatrix(const SparseMatrixAdapter&) override
    {
        // Handled in prepare().
    }

    bool solve(Vector& x) override
    {
        OPM_TIMEBLOCK(istlSolverGpuSystemSolve);
        ++solveCount_;

        const std::size_t numResScalars = gpuA_->N() * numResDofs;
        const std::size_t numWellScalars = mergedD_.N() * numWellDofs;
        if (!sysX_ || sysX_->res.dim() != numResScalars || sysX_->well.dim() != numWellScalars) {
            sysX_.emplace(numResScalars, numWellScalars);
        }
        *sysX_ = 0.0;

        Dune::InverseOperatorResult result;
        sysSolver_->apply(*sysX_, *sysRhs_, result);
        iterations_ = result.iterations;

        sysX_->res.copyToHost(x);

        return Parent::checkConvergence(result, parameters_);
    }

    int iterations() const override
    {
        return iterations_;
    }

    const CommunicationType* comm() const override
    {
        return comm_.get();
    }

    int getSolveCount() const override
    {
        return solveCount_;
    }

private:
    const Simulator& simulator_;
    FlowLinearSolverParameters parameters_;
    PropertyTree prm_;
    std::shared_ptr<CommunicationType> comm_;
    ElementChunksType elementChunks_;
    int iterations_ = 0;
    int solveCount_ = 0;

    bool sysInitialized_ = false;
    WellMatrixStructure cachedWellStructure_;

    // Current per-well B/C/D blocks and their merged host versions.  The
    // merged D also serves the CPU well sub-solver.
    std::vector<WRMatrix<Scalar>> wellBMatrices_;
    std::vector<RWMatrix<Scalar>> wellCMatrices_;
    std::vector<WWMatrix<Scalar>> wellDMatrices_;
    Opm::SparseTable<int> wellCells_;
    WRMatrix<Scalar> mergedB_;
    RWMatrix<Scalar> mergedC_;
    WWMatrix<Scalar> mergedD_;

    // Device copies of the four blocks; GpuSystemMatrix points into these.
    std::optional<GpuRRMatrix<Scalar>> gpuA_;
    std::optional<GpuRWMatrix<Scalar>> gpuC_;
    std::optional<GpuWRMatrix<Scalar>> gpuB_;
    std::optional<GpuWWMatrix<Scalar>> gpuD_;
    std::vector<Scalar> wellBlockValues_;

    GpuSystemMatrix<Scalar> sysMatrix_;
    std::optional<GpuSystemVector<Scalar>> sysX_;
    std::optional<GpuSystemVector<Scalar>> sysRhs_;

    std::unique_ptr<SysOp> sysOp_;
    std::unique_ptr<SysSolver> sysSolver_;
    SysPrecond* sysPrecond_ = nullptr;

    // CPR weights of the reservoir block.
    Vector cpuWeights_;
    std::optional<GpuVectorType> gpuWeights_;
    std::optional<GpuVector<int>> gpuDiagIdx_;

    void prepareSystemSolver(const Matrix& M, const Vector& b)
    {
        OPM_TIMEBLOCK(flexibleSolverPrepare);

        wellBMatrices_.clear();
        wellCMatrices_.clear();
        wellDMatrices_.clear();
        wellCells_.clear();

        simulator_.problem().wellModel().addBCDMatrix(
            wellBMatrices_, wellCMatrices_, wellDMatrices_, wellCells_);

        const Opm::WellMatrixMerger<Scalar> merger(
            M.N(), wellBMatrices_, wellCMatrices_, wellDMatrices_, wellCells_);

        const bool structureChanged = !sysInitialized_
            || !merger.hasSameStructure(cachedWellStructure_);

        if (structureChanged) {
            merger.buildMatrices(mergedB_, mergedC_, mergedD_);
            cachedWellStructure_ = merger.buildStructure();
        } else {
            // Pattern unchanged: write fresh values into the existing merged
            // matrices without any (de)allocation.
            merger.updateValues(mergedB_, mergedC_, mergedD_);
        }

        uploadMatrices(M, structureChanged);
        uploadRhs(b);

        if (!sysInitialized_) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            createSystemSolver();
            sysInitialized_ = true;
        } else if (structureChanged) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            sysPrecond_->updateForChangedWellStructure();
        } else {
            OPM_TIMEBLOCK(flexibleSolverUpdate);
            sysPrecond_->update();
        }
    }

    // Copy the current blocks to the device.  The reservoir block keeps its
    // pattern for the whole run; the well blocks are reallocated only when
    // the well structure changes.
    void uploadMatrices(const Matrix& M, bool structureChanged)
    {
        if (!gpuA_) {
            gpuA_.emplace(GpuRRMatrix<Scalar>::fromMatrix(M));
        } else {
            gpuA_->updateNonzeroValues(M);
        }

        if (structureChanged) {
            gpuB_.emplace(makeScalarGpuMatrix(mergedB_));
            gpuC_.emplace(makeScalarGpuMatrix(mergedC_));
            gpuD_.emplace(makeScalarGpuMatrix(mergedD_));
        } else {
            updateScalarGpuMatrix(mergedB_, *gpuB_, wellBlockValues_);
            updateScalarGpuMatrix(mergedC_, *gpuC_, wellBlockValues_);
            updateScalarGpuMatrix(mergedD_, *gpuD_, wellBlockValues_);
        }

        sysMatrix_.A = &*gpuA_;
        sysMatrix_.B = &*gpuB_;
        sysMatrix_.C = &*gpuC_;
        sysMatrix_.D = &*gpuD_;
        *sysMatrix_.hostD = &mergedD_;
    }

    // As in ISTLSolverSystem, the well part of the right-hand side is zero.
    void uploadRhs(const Vector& b)
    {
        const std::size_t numResScalars = b.size() * numResDofs;
        const std::size_t numWellScalars = mergedD_.N() * numWellDofs;
        if (!sysRhs_ || sysRhs_->res.dim() != numResScalars || sysRhs_->well.dim() != numWellScalars) {
            sysRhs_.emplace(numResScalars, numWellScalars);
        }
        sysRhs_->res.copyFromHost(b);
        if (numWellScalars > 0) {
            sysRhs_->well = 0.0;
        }
    }

    void createSystemSolver()
    {
        const auto resSolverPrm = prm_.get_child("preconditioner.reservoir_solver");
        auto sysWeightCalc = makeWeightsCalculator(resSolverPrm);

        sysOp_ = std::make_unique<SysOp>(sysMatrix_);
        sysSolver_ = std::make_unique<SysSolver>(*sysOp_, prm_, sysWeightCalc, pressureIndex);
        sysPrecond_ = dynamic_cast<SysPrecond*>(&sysSolver_->preconditioner());
        if (!sysPrecond_) {
            OPM_THROW(std::invalid_argument,
                      "The GPU system solver requires preconditioner.type = system_cpr.");
        }
    }

    // CPR weights for the reservoir block, returned in the reservoir part of
    // a system vector; the well part is left empty.
    std::function<GpuSystemVector<Scalar>()>
    makeWeightsCalculator(const PropertyTree& resSolverPrm)
    {
        using namespace std::string_literals;
        const auto weightsType = resSolverPrm.get("preconditioner.weight_type"s, "quasiimpes"s);

        if (weightsType == "quasiimpes") {
            gpuWeights_.emplace(gpuA_->N() * numResDofs);
            gpuDiagIdx_.emplace(Amg::precomputeDiagonalIndices(*gpuA_));
            return [this]() {
                Amg::getQuasiImpesWeights<Scalar, /*transpose=*/false>(
                    *gpuA_, pressureIndex, *gpuWeights_, *gpuDiagIdx_);
                return weightsAsSystemVector();
            };
        }

        if (weightsType == "trueimpes" || weightsType == "trueimpesanalytic") {
            cpuWeights_.resize(gpuA_->N());
            gpuWeights_.emplace(cpuWeights_);
            const bool analytic = weightsType == "trueimpesanalytic";
            const bool threadParallel = parameters_.cpr_weights_thread_parallel_;
            return [this, analytic, threadParallel]() {
                ElementContext elemCtx(simulator_);
                if (analytic) {
                    Amg::getTrueImpesWeightsAnalytic(pressureIndex, cpuWeights_, elemCtx,
                                                     simulator_.model(), elementChunks_,
                                                     threadParallel);
                } else {
                    Amg::getTrueImpesWeights(pressureIndex, cpuWeights_, elemCtx,
                                             simulator_.model(), elementChunks_,
                                             threadParallel);
                }
                gpuWeights_->copyFromHost(cpuWeights_);
                return weightsAsSystemVector();
            };
        }

        OPM_THROW(std::invalid_argument,
                  "Weights type " + weightsType + " not implemented for the GPU system solver.");
    }

    GpuSystemVector<Scalar> weightsAsSystemVector() const
    {
        GpuSystemVector<Scalar> w(gpuWeights_->dim(), 0);
        w.res = *gpuWeights_;
        return w;
    }
};

} // namespace Opm::gpuistl

#endif // OPM_ISTLSOLVERGPUSYSTEM_HEADER_INCLUDED
