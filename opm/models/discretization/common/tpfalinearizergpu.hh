/*
  Copyright 2026 Equinor AS

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
#ifndef TPFA_LINEARIZER_GPU_HH
#define TPFA_LINEARIZER_GPU_HH

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/tpfalinearizerstructs.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/linalg/linalgproperties.hh>

#if HAVE_CUDA
#include <opm/models/discretization/common/tpfalinearizergpukernels.hh>
#include <opm/models/discretization/common/tpfalinearizergpuparams.hh>
#include <opm/simulators/flow/SimpleFIBlackOilModel.hpp>
#include <opm/simulators/flow/ThermalGasWaterFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl_assembly_if_available.hpp>
#if OPM_IS_COMPILING_WITH_GPU_COMPILER
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/detail/gpusparse_matrix_utilities.hpp>
#else
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_utilities.hpp>
#endif
#endif
#include <opm/simulators/wells/WellMatrixSink.hpp>
#endif

#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <dune/common/fvector.hh>

#include <cassert>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace Opm::Parameters
{
struct EnableGpuAssemblySolverBridge {
    static constexpr bool value = true;
};
} // namespace Opm::Parameters

namespace Opm
{
template <class TypeTag>
class TpfaLinearizer;

#if HAVE_CUDA
namespace detail {

template<class Block>
struct GpuWellSystemEntry
{
    std::size_t index;
    Block value;
};

#if OPM_IS_COMPILING_WITH_GPU_COMPILER
template<int numEq, class Scalar, class MatrixBlock, class VectorBlock>
__global__ void applyGpuWellContributions(
    const GpuWellSystemEntry<MatrixBlock>* matrixDeltas,
    std::size_t matrixCount,
    const GpuWellSystemEntry<VectorBlock>* residualEntries,
    std::size_t residualCount,
    Scalar* matrix,
    VectorBlock* residual)
{
    const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < matrixCount) {
        const auto& [blockIndex, block] = matrixDeltas[i];
        auto* values = matrix + blockIndex * numEq * numEq;
        for (int row = 0; row < numEq; ++row) {
            for (int col = 0; col < numEq; ++col) {
                values[row * numEq + col] += block[row][col];
            }
        }
    }
    if (i < residualCount) {
        const auto& [cell, value] = residualEntries[i];
        residual[cell] = value;
    }
}
#endif

} // namespace detail
#endif

// Owns the GPU assembly lifecycle. The layout is identical in host and device
// translation units; only GPU compilers instantiate kernel-dependent methods.
template <class TypeTag>
class TpfaLinearizerGpu
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    static constexpr int numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool runAssemblyOnGpu = getPropValue<TypeTag, Properties::RunAssemblyOnGpu>();
    using MatrixBlockCPU = typename SparseMatrixAdapter::MatrixBlock;
    using VectorBlockCPU = Dune::FieldVector<Scalar, numEq>;
    using ADVectorBlockCPU = GetPropType<TypeTag, Properties::RateVector>;

public:
#if HAVE_CUDA
    using MatrixBlockGPU = gpuistl::MiniMatrix<Scalar, numEq>;
    using VectorBlockGPU = gpuistl::MiniVector<Scalar, numEq>;
    using ADVectorBlockGPU = gpuistl::MiniVector<Evaluation, numEq>;
    using MatrixBlock = std::conditional_t<runAssemblyOnGpu, MatrixBlockGPU, MatrixBlockCPU>;
    using VectorBlock = std::conditional_t<runAssemblyOnGpu, VectorBlockGPU, VectorBlockCPU>;
    using ADVectorBlock = std::conditional_t<runAssemblyOnGpu, ADVectorBlockGPU, ADVectorBlockCPU>;
#else
    using MatrixBlock = MatrixBlockCPU;
    using VectorBlock = VectorBlockCPU;
    using ADVectorBlock = ADVectorBlockCPU;
#endif

    struct WellContributionCounts {
        std::size_t matrixBlocks = 0;
        std::size_t residualRows = 0;
    };

    void eraseMatrix()
    {
#if HAVE_CUDA
        gpuParams_.reset();
        gpuBufferDiagMatAddress_.reset();
        gpuJacobian_.reset();
#endif
        gpuSystemActive_ = false;
    }

    void createMatrix([[maybe_unused]] const SparseMatrixAdapter& jacobian)
    {
        if constexpr (runAssemblyOnGpu) {
            requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
            const auto [columns, rows] = gpuistl::detail::extractSparsityPattern(jacobian.istlMatrix());
            gpuJacobian_ = std::make_unique<gpuistl::GpuSparseMatrixWrapper<Scalar>>(
                gpuistl::GpuVector<int>(rows), gpuistl::GpuVector<int>(columns), numEq);
            gpuBufferDiagMatAddress_ = std::make_unique<gpuistl::GpuBuffer<MatrixBlockGPU*>>(
                gpuistl::detail::getDiagPtrsTyped<MatrixBlockGPU>(*gpuJacobian_));
#endif
        }
    }

    void resetSystem()
    {
        if constexpr (runAssemblyOnGpu) {
            requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
            gpuJacobian_->setToZero();
#endif
        }
    }

    [[nodiscard]] bool gpuSystemActive() const
    {
        return gpuSystemActive_;
    }

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    bool usesGpuPropertyAssemblyBridge() const
    {
        return gpuParams_ && gpuParams_->usesGpuPropertyAssemblyBridge();
    }
#endif

    decltype(auto) gpuJacobian()
    {
        requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
        return gpuParams_->gpuJacobian();
#else
        OPM_THROW(std::logic_error, "GPU matrix requires GPU support");
#endif
    }

    decltype(auto) gpuResidual()
    {
        requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
        return gpuParams_->residualBuffer();
#else
        OPM_THROW(std::logic_error, "GPU residual requires GPU support");
#endif
    }

    decltype(auto) flattenedGpuResidual()
    {
        requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
        return gpuParams_->flattenedResidual();
#else
        OPM_THROW(std::logic_error, "GPU residual requires GPU support");
#endif
    }

    void copyJacobianToHost([[maybe_unused]] SparseMatrixAdapter& jacobian)
    {
        if constexpr (runAssemblyOnGpu) {
            requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
            if (gpuSystemActive_) {
                gpuParams_->copyJacobianToHost(jacobian, *gpuJacobian_);
                OpmLog::debug("[GPU TPFA] full Jacobian transfer for system export");
            }
#endif
        }
    }

    [[nodiscard]] bool updateSolverBridge([[maybe_unused]] TpfaLinearizer<TypeTag>& linearizer)
    {
        gpuSystemActive_ = false;
        if constexpr (runAssemblyOnGpu) {
            requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
            // Multiple ranks still need a host matrix for makeOverlapRowsInvalid.
            const auto& wellModel = linearizer.problem_().wellModel();
            bool supportsWellSystem = wellModel.numLocalWellsEnd() == 0;
            if constexpr (requires {
                linearizer.model_().nohystGpuProperties();
                wellModel.empty();
                wellModel.addMatrixContributions();
            }) {
                supportsWellSystem = wellModel.empty() || wellModel.addMatrixContributions();
            }
            gpuSystemActive_ = Parameters::Get<Parameters::EnableGpuAssemblySolverBridge>()
                && !getPropValue<TypeTag, Properties::EnableConstraints>()
                && linearizer.simulator_().gridView().comm().size() == 1
                && supportsWellSystem
                && linearizer.simulator_().model().newtonMethod().linearSolver().hasGpuSolver();
#endif
        }
        return gpuSystemActive_;
    }

    WellContributionCounts linearizeWellContributions(
        [[maybe_unused]] TpfaLinearizer<TypeTag>& linearizer)
    {
        requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
        auto& wellModel = linearizer.problem_().wellModel();
        using Sink = WellMatrixSink<MatrixBlockCPU>;
        if constexpr (requires(Sink& sink, std::set<unsigned>& rows) {
            wellModel.linearizeWithMatrixSink(sink, linearizer.residual_, rows);
            wellModel.empty();
        }) {
            if (wellModel.empty()) {
                return {};
            }

            // GPU BSR has the same block order as the contiguous Dune BCRS
            // storage used to create it and to map the TPFA face blocks.
            const auto& matrix = linearizer.jacobian_->istlMatrix();
            const auto* firstBlock = &*matrix[0].begin();
            std::map<std::size_t, MatrixBlockCPU> matrixDeltas;
            Sink sink([&](std::size_t row, std::size_t col, const MatrixBlockCPU& value) {
                const auto entry = matrix[row].find(col);
                if (entry == matrix[row].end()) {
                    OPM_THROW(std::logic_error, "Missing TPFA matrix block for a well contribution");
                }
                const std::size_t index = &*entry - firstBlock;
                auto& block = matrixDeltas.try_emplace(index, Scalar{0}).first->second;
                block += value;
            });

            std::set<unsigned> residualRows;
            wellModel.linearizeWithMatrixSink(sink, linearizer.residual_, residualRows);
            if (matrixDeltas.empty() && residualRows.empty()) {
                return {};
            }

            wellMatrixStaging_.clear();
            wellMatrixStaging_.reserve(matrixDeltas.size());
            for (const auto& [index, block] : matrixDeltas) {
                wellMatrixStaging_.push_back({index, MatrixBlockGPU(block)});
            }
            wellResidualStaging_.clear();
            wellResidualStaging_.reserve(residualRows.size());
            for (const auto cell : residualRows) {
                wellResidualStaging_.push_back({cell, VectorBlockGPU(linearizer.residual_[cell])});
            }

            uploadWellEntries(wellMatrixStaging_, wellMatrixEntries_);
            uploadWellEntries(wellResidualStaging_, wellResidualEntries_);

            constexpr unsigned blockSize = 256;
            const auto workItems = std::max(matrixDeltas.size(), residualRows.size());
            detail::applyGpuWellContributions<numEq>
                <<<(workItems + blockSize - 1) / blockSize, blockSize>>>(
                    wellMatrixEntries_.data(), matrixDeltas.size(),
                    wellResidualEntries_.data(), residualRows.size(),
                    gpuJacobian_->getNonZeroValues().data(), gpuParams_->residualBuffer().data());
            OPM_CUDA_CHECK_LAST_ERROR;
            return {matrixDeltas.size(), residualRows.size()};
        }
#endif
        return {};
    }

    template <class SubDomainType>
    void linearize([[maybe_unused]] TpfaLinearizer<TypeTag>& linearizer,
                   [[maybe_unused]] const SubDomainType& domain,
                   [[maybe_unused]] unsigned int numCells,
                   [[maybe_unused]] Scalar dt,
                   [[maybe_unused]] bool dispersionActive)
    {
        requireGpuSupport();
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
        const bool useAssemblyLinearsolveBridge = gpuSystemActive_;
        // GPU assembly currently supports only the full domain.
        if constexpr (std::is_same_v<SubDomainType, FullDomain<>>) {

            if constexpr (requires { linearizer.model_().nohystGpuProperties(); }) {
                if (dispersionActive || !linearizer.boundaryInfo_.empty()) {
                    OPM_THROW(std::logic_error,
                              "Three-phase GPU reservoir assembly requires no dispersion or nontrivial boundary conditions.");
                }
                if (!gpuParams_) {
                    auto params = std::make_unique<GpuParams>(
                        domain, linearizer.neighborInfo_, *gpuJacobian_, *gpuBufferDiagMatAddress_,
                        linearizer.jacobian_->istlMatrix(), linearizer.residual_);
                    gpuParams_ = GpuParamsOwner(params.release(), [](GpuParams* ptr) {
                        std::default_delete<GpuParams> {}(ptr);
                    });
                }
                gpuParams_->resetResidual(numCells);
                linearizer.model_().nohystGpuProperties().assemble(*gpuParams_, linearizer.problem_(), dt);
            } else {
                int constexpr blockSize = 256; // Experimentally this is a good value for multiple
                                               // GPUs. Autotune this later.

                if (!gpuParams_) {
                    auto params = std::make_unique<GpuParams>(domain,
                                                              linearizer.neighborInfo_,
                                                              *gpuJacobian_,
                                                              *gpuBufferDiagMatAddress_,
                                                              linearizer.jacobian_->istlMatrix(),
                                                              linearizer.residual_,
                                                              linearizer.boundaryInfo_,
                                                              linearizer.model_(),
                                                              linearizer.problem_(),
                                                              numCells);
                    gpuParams_ = GpuParamsOwner(params.release(), [](GpuParams* ptr) {
                        std::default_delete<GpuParams> {}(ptr);
                    });
                } else {
                    gpuParams_->refreshForLinearization(
                        linearizer.model_(), linearizer.boundaryInfo_, numCells);
                }

                linearizer.prepare_source_intensive_quantities_gpu(numCells, domain);

                launchLinearization<typename GpuParams::LocalResidualGPU>(
                    numCells,
                    gpuParams_->domainView(),
                    gpuParams_->neighborInfoView(),
                    gpuParams_->diagMatAddressView(),
                    gpuParams_->residualView(),
                    gpuParams_->modelView(),
                    dt,
                    dispersionActive,
                    gpuParams_->flowProblemView());

                if (gpuParams_->boundaryInfoSize() > 0) {
                    auto boundaryInfoView = gpuParams_->boundaryInfoView();
                    linearize_bc_threadsafe<TpfaLinearizerGpu<TypeTag>,
                                            typename GpuParams::GPUBOIQ,
                                            decltype(gpuParams_->modelView()),
                                            typename GpuParams::LocalResidualGPU>
                        <<<((gpuParams_->boundaryInfoSize() + blockSize - 1) / blockSize), blockSize>>>(
                            gpuParams_->diagMatAddressView(),
                            gpuParams_->residualView(),
                            boundaryInfoView,
                            gpuParams_->modelView(),
                            gpuParams_->flowProblemView());
                }
            }

            if (useAssemblyLinearsolveBridge) {
                linearizeSourceTerms(linearizer, numCells, domain);

                // The convergence report consumes the residual on the host,
                // but the numerical Jacobian remains device authoritative.
                gpuParams_->copyResidualToHost(linearizer.residual_, numCells);
            } else {
                gpuParams_->copyResidualToHost(linearizer.residual_, numCells);
                gpuParams_->copyJacobianToHost(*linearizer.jacobian_, *gpuJacobian_);
                OpmLog::debug("[GPU TPFA] full Jacobian transfer for solver fallback");
            }
        } else {
            OPM_THROW(std::logic_error, "Only FullDomain is supported on GPU");
        }
#else
        OPM_THROW(std::logic_error,
                  "Trying to run GPU assembly without compiling with GPU support");
#endif
    }
#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    // Thread-safe on the GPU, but floating-point atomic additions make the numerical
    // result depend on the order in which boundary faces are calculated.
    template <class LocalIntensiveQuantities,
              class LocalResidualT,
              class ModelT,
              class DiagPtrType,
              class ResidualT,
              class BoundaryInfoT,
              class ProblemT>
    OPM_HOST_DEVICE static void
    linearize_bc_threadsafe_single_cell(DiagPtrType diagonalAddresses,
                                        ResidualT residual,
                                        const BoundaryInfoT boundaryInfoElement,
                                        ModelT model,
                                        ProblemT problem)
    {
        if (boundaryInfoElement.bcdata.type != BCType::NONE) {
            VectorBlock residualContribution(0.0);
            MatrixBlock jacobianContribution(0.0);
            ADVectorBlock flux(0.0);
            const unsigned cell = boundaryInfoElement.cell;
            const LocalIntensiveQuantities& insideIntensiveQuantities
                = model.intensiveQuantities(cell, /*timeIdx*/ 0);
            // LocalResidualT carries the device fluid system. The linearizer's
            // TypeTag still names the host fluid system, so it cannot gate this call.
            LocalResidualT::computeBoundaryFlux(
                flux, problem, boundaryInfoElement.bcdata, insideIntensiveQuantities, cell);
            flux *= boundaryInfoElement.bcdata.faceArea;
            TpfaLinearizer<TypeTag>::setResAndJacobi(
                residualContribution, jacobianContribution, flux);
            auto* residualBlock = &(residual.data()[cell]);
            for (int i = 0; i < numEq; ++i) {
                atomicAdd(&((*residualBlock)[i]), residualContribution[i]);
            }
            auto* matrixBlock = diagonalAddresses[cell];
            for (int row = 0; row < jacobianContribution.size(); ++row) {
                for (int col = 0; col < jacobianContribution.size(); ++col) {
                    Scalar* element = &((*matrixBlock)[row][col]);
                    atomicAdd(element, jacobianContribution[row][col]);
                }
            }
        }
    }
#endif

private:
    // CPU assembly also calls the lifecycle methods, which return before this
    // check. Actual GPU operations must never silently succeed without support.
    static void requireGpuSupport()
    {
        if constexpr (!runAssemblyOnGpu) {
            OPM_THROW(std::logic_error, "GPU assembly is disabled for this type tag");
        }
#if !HAVE_CUDA
        OPM_THROW(std::logic_error, "GPU assembly requires CUDA or HIP support");
#elif !OPM_IS_COMPILING_WITH_GPU_COMPILER
        OPM_THROW(std::logic_error, "GPU assembly requires compilation with a GPU compiler");
#endif
    }

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    template <class LocalResidualT,
              class ModelClass,
              class DiagPtrType,
              class DomainType,
              class NeighborSparseTable,
              class ResidualType,
              class ProblemT>
    void launchLinearization(const unsigned int numCells,
                             const DomainType& domain,
                             const NeighborSparseTable& neighborInfo,
                             DiagPtrType& diagMatAddress,
                             ResidualType& residual,
                             const ModelClass& model,
                             Scalar dt,
                             [[maybe_unused]] bool dispersionActive,
                             const ProblemT& problem)
    {
        static_assert(!getPropValue<TypeTag, Properties::EnableBioeffects>(),
                      "Bioeffects not yet supported on GPU");
        assert(!dispersionActive && "Dispersion not yet supported on GPU");
        int constexpr blockSize = 256;
        kernel_linearize<TpfaLinearizer<TypeTag>,
                         ModelClass,
                         LocalResidualT,
                         DiagPtrType,
                         Scalar,
                         DomainType,
                         NeighborSparseTable,
                         ResidualType,
                         ProblemT><<<((numCells + blockSize - 1) / blockSize), blockSize>>>(
            numCells, domain, neighborInfo, diagMatAddress, residual, model, dt, problem);
    }
    template <class SubDomainType>
    void linearizeSourceTerms(TpfaLinearizer<TypeTag>& linearizer,
                              unsigned int numCells,
                              const SubDomainType& domain)
    {
        // Evaluate the same CPU source calculation as the host path, uploading
        // only nonzero residual and diagonal-block contributions to gpuParams_.
        std::vector<unsigned> indices;
        std::vector<VectorBlockCPU> residualContributions;
        std::vector<MatrixBlockCPU> jacobianContributions;
        // Small initial capacity; larger source sets grow normally.
        constexpr std::size_t reservedSize = 20;
        indices.reserve(reservedSize);
        residualContributions.reserve(reservedSize);
        jacobianContributions.reserve(reservedSize);
        const bool needsDenseSources = [&] {
            if constexpr (!getPropValue<TypeTag, Properties::EnableBioeffects>()
                          && requires { linearizer.problem_().aquiferModel().hasAquifers(); }) {
                if (linearizer.separateSparseSourceTerms_) {
                    const auto& simulator = linearizer.simulator_();
                    const auto& schedule = simulator.vanguard().schedule();
                    return linearizer.problem_().aquiferModel().hasAquifers()
                        || schedule[simulator.episodeIndex()].source().size() != 0;
                }
            }
            return true;
        }();
        const auto denseCellCount = needsDenseSources ? numCells : 0u;
        for (unsigned ii = 0; ii < denseCellCount; ++ii) {
            const unsigned cell = domain.cells[ii];
            ADVectorBlockCPU source(0.0);
            const double volume = linearizer.model_().dofTotalVolume(cell);
            if constexpr (requires { LocalResidual::computeSource(source, linearizer.problem_(), cell, 0); }) {
                if (linearizer.separateSparseSourceTerms_) {
                    LocalResidual::computeSourceDense(source, linearizer.problem_(), cell, 0);
                } else {
                    LocalResidual::computeSource(source, linearizer.problem_(), cell, 0);
                }
            } else {
                const auto& iq = linearizer.model_().intensiveQuantities(cell, 0);
                if (linearizer.separateSparseSourceTerms_) {
                    LocalResidual::computeSourceDense(source, linearizer.problem_(), iq, cell, 0);
                } else {
                    LocalResidual::computeSource(source, linearizer.problem_(), iq, cell, 0);
                }
            }
            source *= -volume;
            VectorBlockCPU residualContribution(0.0);
            MatrixBlockCPU jacobianContribution(0.0);
            TpfaLinearizer<TypeTag>::setResAndJacobi(
                residualContribution, jacobianContribution, source);
            bool nonzero = false;
            for (unsigned eq = 0; eq < numEq; ++eq) {
                nonzero = nonzero || residualContribution[eq] != 0.0;
                for (unsigned pv = 0; pv < numEq; ++pv) {
                    nonzero = nonzero || jacobianContribution[eq][pv] != 0.0;
                }
            }
            if (nonzero) {
                indices.push_back(cell);
                residualContributions.push_back(residualContribution);
                jacobianContributions.push_back(jacobianContribution);
            }
        }
        gpuParams_->applySourceContributions(indices, residualContributions, jacobianContributions);

        // Dense source evaluation excludes wells in the sparse mode. Accumulate
        // shared completion cells before uploading, so each device row has one writer.
        const auto& wellModel = linearizer.problem_().wellModel();
        const auto addSparseSource = [&](unsigned cell, const VectorBlockCPU& residual,
                                         const MatrixBlockCPU& jacobian) {
            auto& blocks = wellSourceEntries_.try_emplace(cell, Scalar{0}, Scalar{0}).first->second;
            blocks.first += residual;
            blocks.second += jacobian;
        };
        if constexpr (requires { wellModel.addReservoirSourceTerms(addSparseSource); }) {
            if (linearizer.separateSparseSourceTerms_) {
                wellSourceEntries_.clear();
                wellModel.addReservoirSourceTerms(addSparseSource);
                indices.clear();
                residualContributions.clear();
                jacobianContributions.clear();
                for (const auto& [cell, blocks] : wellSourceEntries_) {
                    indices.push_back(cell);
                    residualContributions.push_back(blocks.first);
                    jacobianContributions.push_back(blocks.second);
                }
                gpuParams_->applySourceContributions(indices, residualContributions, jacobianContributions);
            }
        }
    }

    template<class Entry>
    static void uploadWellEntries(const std::vector<Entry>& staging, gpuistl::GpuBuffer<Entry>& buffer)
    {
        if (!staging.empty()) {
            if (buffer.size() < staging.size()) {
                buffer = gpuistl::GpuBuffer<Entry>(staging.size());
            }
            buffer.copyFromHost(staging.data(), staging.size());
        }
    }
#endif

#if HAVE_CUDA
    using WellMatrixEntry = detail::GpuWellSystemEntry<MatrixBlockGPU>;
    using WellResidualEntry = detail::GpuWellSystemEntry<VectorBlockGPU>;
    std::vector<WellMatrixEntry> wellMatrixStaging_;
    std::vector<WellResidualEntry> wellResidualStaging_;
    gpuistl::GpuBuffer<WellMatrixEntry> wellMatrixEntries_;
    gpuistl::GpuBuffer<WellResidualEntry> wellResidualEntries_;
    std::map<unsigned, std::pair<VectorBlockCPU, MatrixBlockCPU>> wellSourceEntries_;

    std::unique_ptr<gpuistl::GpuSparseMatrixWrapper<Scalar>> gpuJacobian_;
    std::unique_ptr<gpuistl::GpuBuffer<MatrixBlockGPU*>> gpuBufferDiagMatAddress_;
    // Type-erased destruction permits an incomplete state in host translation units.
    using GpuParams = TpfaLinearizerGpuParams<TypeTag>;
    using GpuParamsOwner = std::unique_ptr<GpuParams, void (*)(GpuParams*)>;
    GpuParamsOwner gpuParams_ {nullptr, nullptr};
#endif
    bool gpuSystemActive_ = false;
};
} // namespace Opm
#endif
