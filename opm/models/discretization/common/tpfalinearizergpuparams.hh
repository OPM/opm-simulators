// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*

  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * \brief GPU parameter setup class for TpfaLinearizer.
 *
 * \note This file is designed to be included from tpfalinearizer.hh inside a
 *       \c #if HAVE_CUDA block.  It is self-contained: all required GPU headers,
 *       model headers, and forward declarations are included below.
 */
#ifndef TPFA_LINEARIZER_GPU_PARAMS_HH
#define TPFA_LINEARIZER_GPU_PARAMS_HH

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER

#include <memory>
#include <string>
#include <vector>

#include <opm/common/utility/gpuistl_if_available.hpp>

#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl_hip/MiniMatrix.hpp>
#include <opm/simulators/linalg/gpuistl_hip/MiniVector.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/MiniMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#endif

#include <opm/grid/utility/SparseTable.hpp>
#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/simulators/flow/SimpleFIBlackOilModel.hpp>
#include <opm/simulators/flow/ThermalGasWaterFlowProblem.hpp>

#include <opm/models/discretization/common/tpfalinearizerstructs.hh>

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Manages GPU buffer allocation and view setup for one linearization pass.
 *
 * Constructed once per call to TpfaLinearizer::linearize_() when the GPU assembly
 * path is active.  The constructor performs all CPU-to-GPU copies and stores the
 * resulting GPU buffers as members.  Callers retrieve lightweight GpuView objects
 * via the accessor methods and pass them directly to the GPU kernels.  After the
 * kernels complete, copyResidualToHost() and copyJacobianToHost() transfer the
 * computed results back to the CPU.
 */
template <class TypeTag>
class TpfaLinearizerGpuParams
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };

    using MatrixBlockCPU = typename SparseMatrixAdapter::MatrixBlock;
    using VectorBlockCPU = Dune::FieldVector<Scalar, numEq>;
    using MatrixBlockGPU = gpuistl::MiniMatrix<Scalar, numEq>;
    using VectorBlockGPU = gpuistl::MiniVector<Scalar, numEq>;

    using ResidualNBInfo = typename LocalResidual::ResidualNBInfo;
    using NeighborInfoCPU = NeighborInfoStruct<ResidualNBInfo, MatrixBlockCPU>;
    using NeighborInfoGPU = NeighborInfoStruct<ResidualNBInfo, MatrixBlockGPU>;

    using ScalarFluidState = typename IntensiveQuantities::ScalarFluidState;
    using BoundaryConditionDataCPU = BoundaryConditionData<VectorBlockCPU, ScalarFluidState>;
    using BoundaryInfoCPU = BoundaryInfo<BoundaryConditionDataCPU>;

    // Fluid system buffer/view/ptr types derived without spelling out the fluid-system
    // template parameters explicitly.
    using NonStaticFluidSystem = std::decay_t<decltype(FluidSystem::getNonStaticInstance())>;
    using GpuFluidSystemBuffer = std::decay_t<decltype(::Opm::gpuistl::copy_to_gpu(
        std::declval<NonStaticFluidSystem&>()))>;
    using GpuFluidSystemView
        = std::decay_t<decltype(::Opm::gpuistl::make_view(std::declval<GpuFluidSystemBuffer&>()))>;
    using GpuFluidSystemPtr = std::shared_ptr<GpuFluidSystemView>;

public:
    // Type aliases required by callers to instantiate the GPU kernel templates.
    using CorrectTypeTagView =
        typename ::Opm::Properties::TTag::to_gpu_type_t<TypeTag, gpuistl::GpuView>;

    // GPUBOIQ and LocalResidualGPU are not taken from typetag as we need the
    // outer class to template correctly here
    using GPUBOIQ = BlackOilIntensiveQuantities<CorrectTypeTagView>;
    using LocalResidualGPU = BlackOilLocalResidualTPFA<CorrectTypeTagView>;

private:
    using GpuScalarFluidState = typename GPUBOIQ::ScalarFluidState;
    using BoundaryConditionDataGPU = BoundaryConditionData<VectorBlockGPU, GpuScalarFluidState>;
    using BoundaryInfoGPU = BoundaryInfo<BoundaryConditionDataGPU>;

    using GpuModelBufferType = SimpleFIBlackOilModel<CorrectTypeTagView, gpuistl::GpuBuffer>;
    using GpuFlowProblemBufferType = ThermalGasWaterFlowProblem<Scalar, gpuistl::GpuBuffer>;

    using GpuModel = GetPropType<TypeTag, Properties::GpuFIBlackOilModel>;

    FullDomain<gpuistl::GpuBuffer<int>> domainBuffer_;
    SparseTable<NeighborInfoGPU, gpuistl::GpuBuffer> neighborInfoBuffer_;
    // diagMatAddressView_ is non-owning: the underlying buffer lives in TpfaLinearizer.
    gpuistl::GpuView<MatrixBlockGPU*> diagMatAddressView_;
    gpuistl::GpuBuffer<VectorBlockGPU> residualBuffer_;
    gpuistl::GpuView<VectorBlockGPU>
        residualView_; // stored as lvalue because mutable ref must be provided
    // dynamicGpuFluidSystemBuffer_ must be declared before dynamicGpuFluidSystemPtr_
    // because the ptr holds a GpuView into the buffer's GPU memory.
    GpuFluidSystemBuffer dynamicGpuFluidSystemBuffer_;
    // The fluid-system ptr is kept alive because gpuModelBuffer_ and
    // boundaryInfoBuffer_ store raw GPU pointers into it.
    GpuFluidSystemPtr dynamicGpuFluidSystemPtr_;
    GpuModelBufferType gpuModelBuffer_;
    GpuFlowProblemBufferType gpuFlowProblemBuffer_;
    gpuistl::GpuBuffer<BoundaryInfoGPU> boundaryInfoBuffer_;

public:
    /*!
     * \brief Construct and populate all GPU buffers needed for one linearization pass.
        To avoid optionals we use initializerlists, although some are a bit ugly
        with the lambdas
     *
     * \param domain                  The full domain (CPU).
     * \param neighborInfo            CPU neighbor-info sparse table.
     * \param gpuJacobian             Persistent GPU sparse-matrix wrapper (owned by
     * TpfaLinearizer).
     * \param gpuBufferDiagMatAddress Persistent GPU diagonal-pointer buffer (owned by
     * TpfaLinearizer).
     * \param cpuJacobian             CPU ISTL matrix (used to compute GPU pointer offsets).
     * \param residual                CPU residual vector (read to initialise the GPU residual
     * buffer).
     * \param boundaryInfo            CPU boundary-condition info.
     * \param model                   CPU model (provides intensive quantities and volumes).
     * \param problem                 CPU problem (provides thermal transmissibilities and module
     * params).
     * \param numCells                Number of cells in the domain.
     */
    TpfaLinearizerGpuParams(const FullDomain<>& domain,
                            const SparseTable<NeighborInfoCPU>& neighborInfo,
                            gpuistl::GpuSparseMatrixWrapper<Scalar>& gpuJacobian,
                            gpuistl::GpuBuffer<MatrixBlockGPU*>& gpuBufferDiagMatAddress,
                            typename SparseMatrixAdapter::IstlMatrix& cpuJacobian,
                            GlobalEqVector& residual,
                            const std::vector<BoundaryInfoCPU>& boundaryInfo,
                            Model& model,
                            const Problem& problem,
                            unsigned int numCells)
        : domainBuffer_(copy_to_gpu(domain))
        , neighborInfoBuffer_(
              gpuistl::copy_to_gpu<MatrixBlockGPU>(neighborInfo, gpuJacobian, cpuJacobian))
        , diagMatAddressView_(gpuistl::make_view(gpuBufferDiagMatAddress))
        , residualBuffer_(gpuistl::copy_to_gpu_residual<GlobalEqVector, VectorBlockGPU>(residual))
        , residualView_(gpuistl::make_view(residualBuffer_))
        , dynamicGpuFluidSystemBuffer_(
              ::Opm::gpuistl::copy_to_gpu(FluidSystem::getNonStaticInstance()))
        , dynamicGpuFluidSystemPtr_(
              gpuistl::make_gpu_shared_ptr(::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer_)))
        , gpuModelBuffer_([&]() -> GpuModelBufferType {
            std::vector<Scalar> volumes(numCells);
            for (unsigned i = 0; i < numCells; ++i) {
                volumes[domain.cells[i]] = model.dofTotalVolume(domain.cells[i]);
            }
            GpuModel gpuModel(
                model.intensiveQuantityCache()[0], model.intensiveQuantityCache()[1], volumes);
            return gpuistl::copy_to_gpu(gpuModel, *dynamicGpuFluidSystemPtr_.get());
        }())
        , gpuFlowProblemBuffer_([&]() -> GpuFlowProblemBufferType {
            std::vector<Scalar> alpha0(numCells);
            std::vector<Scalar> alpha1(numCells);
            std::vector<Scalar> alpha2(numCells);
            const auto& thermalHalfTransBoundary
                = problem.eclTransmissibilities().getThermalHalfTransBoundary();
            for (const auto& [key, value] : thermalHalfTransBoundary) {
                const unsigned cell = key.first;
                const unsigned dir = key.second;
                if (dir == 0) {
                    alpha0[cell] = value;
                } else if (dir == 1) {
                    alpha1[cell] = value;
                } else if (dir == 2) {
                    alpha2[cell] = value;
                } else {
                    OPM_THROW(std::logic_error,
                              "Invalid direction for thermal half transmissibility: "
                                  + std::to_string(dir));
                }
            }
            ThermalGasWaterFlowProblem<Scalar> gpuFlowProblem(
                alpha0, alpha1, alpha2, problem.moduleParams());
            return gpuistl::copy_to_gpu(gpuFlowProblem);
        }())
        , boundaryInfoBuffer_(
              gpuistl::copy_to_gpu<VectorBlockGPU, GpuScalarFluidState, BoundaryInfoGPU>(
                  boundaryInfo, *dynamicGpuFluidSystemPtr_.get()))
    {
    }

    auto domainView()
    {
        return make_view(domainBuffer_);
    }

    auto neighborInfoView()
    {
        return gpuistl::make_view(neighborInfoBuffer_);
    }

    auto& diagMatAddressView()
    {
        return diagMatAddressView_;
    }

    auto& residualView()
    {
        return residualView_;
    }

    auto modelView()
    {
        return gpuistl::make_view(gpuModelBuffer_);
    }

    auto flowProblemView()
    {
        return gpuistl::make_view(gpuFlowProblemBuffer_);
    }

    auto boundaryInfoView()
    {
        return gpuistl::make_view(boundaryInfoBuffer_);
    }

    std::size_t boundaryInfoSize() const
    {
        return boundaryInfoBuffer_.size();
    }

    /*!
     * \brief Copy the GPU residual buffer back to the CPU residual vector.
     *
     * Must be called after all GPU kernels have completed (the synchronous cudaMemcpy
     * inside guarantees the GPU has finished).
     */
    void copyResidualToHost(GlobalEqVector& residual, unsigned numCells)
    {
        auto cpuResidualFromGpu = residualBuffer_.asStdVector();
        std::memcpy(residual.data(), cpuResidualFromGpu.data(), numCells * numEq * sizeof(Scalar));
    }

    /*!
     * \brief Copy the GPU Jacobian non-zero values back to the CPU ISTL matrix.
     *
     * Must be called after all GPU kernels have completed.
     */
    void copyJacobianToHost(SparseMatrixAdapter& jacobian,
                            gpuistl::GpuSparseMatrixWrapper<Scalar>& gpuJacobian)
    {
        auto gpuJacobianNonZeroes = gpuJacobian.getNonZeroValues().asStdVector();
        auto& cpuJacobian = jacobian.istlMatrix();
        const std::size_t totalMatrixSize = cpuJacobian.nonzeroes() * numEq * numEq;
        std::memcpy(&(cpuJacobian[0][0][0][0]),
                    gpuJacobianNonZeroes.data(),
                    totalMatrixSize * sizeof(Scalar));
    }
};

} // namespace Opm

#endif // HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER

#endif // TPFA_LINEARIZER_GPU_PARAMS_HH
