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

#ifndef TPFA_LINEARIZER_STRUCTS_HH
#define TPFA_LINEARIZER_STRUCTS_HH

/*!
 * Structs needed for tpfalinearizer and its gpuparams struct extracted
 * to be defined in one place that both can include.
 */

#include <cassert>
#include <vector>

#include <opm/common/Exceptions.hpp>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/common/utility/gpuistl_if_available.hpp>
#include <opm/common/utility/pointerArithmetic.hpp>

namespace Opm {

template <class Storage = std::vector<int>>
struct FullDomain
{
    Storage cells;
};

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
inline FullDomain<gpuistl::GpuBuffer<int>> copy_to_gpu(FullDomain<> CPUDomain)
{
    if (CPUDomain.cells.size() == 0) {
        OPM_THROW(std::runtime_error, "Cannot copy empty full domain to GPU.");
    }
    return FullDomain<gpuistl::GpuBuffer<int>>{
        gpuistl::GpuBuffer<int>(CPUDomain.cells)
    };
};

inline FullDomain<gpuistl::GpuView<int>> make_view(FullDomain<gpuistl::GpuBuffer<int>>& buffer)
{
    if (buffer.cells.size() == 0) {
        OPM_THROW(std::runtime_error, "Cannot make view of empty full domain buffer.");
    }
    return FullDomain<gpuistl::GpuView<int>>{
        gpuistl::make_view(buffer.cells)
    };
};
#endif

template <class ResidualNBInfoType, class BlockType>
struct NeighborInfoStruct
{
    unsigned int neighbor;
    ResidualNBInfoType res_nbinfo;
    BlockType* matBlockAddress;

    template <class PtrType>
    NeighborInfoStruct(unsigned int n, const ResidualNBInfoType& r, PtrType ptr)
        : neighbor(n)
        , res_nbinfo(r)
        , matBlockAddress(static_cast<BlockType*>(ptr))
    {
    }

    // Add a default constructor
    NeighborInfoStruct()
        : neighbor(0)
        , res_nbinfo()
        , matBlockAddress(nullptr)
    {
    }
};

template <class VectorBlock, class ScalarFluidState>
struct BoundaryConditionData {
    BCType type;
    VectorBlock massRate;
    unsigned pvtRegionIdx;
    unsigned boundaryFaceIndex;
    double faceArea;
    double faceZCoord;
    ScalarFluidState exFluidState;
};

template <class BoundaryConditionData>
struct BoundaryInfo {
    unsigned int cell;
    int dir;
    unsigned int bfIndex;
    BoundaryConditionData bcdata;
};

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
namespace gpuistl {

    template <class MiniMatrixType,
              class GpuMatrixType,
              class CpuMatrixType,
              class MatrixBlockType,
              class ResidualNBInfoType>
    auto copy_to_gpu(const SparseTable<NeighborInfoStruct<ResidualNBInfoType, MatrixBlockType>>&
                         cpuNeighborInfoTable,
                     GpuMatrixType& gpuJacobian,
                     CpuMatrixType& cpuJacobian)
    {
        // Convert the DUNE FieldVectors to MiniMatrix types
        using StructWithMinimatrix = NeighborInfoStruct<ResidualNBInfoType, MiniMatrixType>;
        using Scalar = typename GpuMatrixType::field_type;
        std::vector<StructWithMinimatrix> minimatrices(cpuNeighborInfoTable.dataSize());
        Scalar* gpuBufStart = gpuJacobian.getNonZeroValues().data();
        Scalar* cpuBufStart = &(cpuJacobian[0][0][0][0]);

        // To compute the length of the buffer of the cpuJacobian we here assume we have a blocked
        // BCRS matrix with square blocks and that the blocks are stored as Dune::FieldMatrix
        using CpuBlockType = typename CpuMatrixType::block_type::BaseType;

        const size_t gpuBufferSize
            = gpuJacobian.nonzeroes() * gpuJacobian.blockSize() * gpuJacobian.blockSize();
        const size_t cpuBufferSize
            = cpuJacobian.nonzeroes() * CpuBlockType::rows * CpuBlockType::cols;
        assert(gpuBufferSize == cpuBufferSize);

        const auto& dataStorage = cpuNeighborInfoTable.dataStorage();
        const size_t dataSize = cpuNeighborInfoTable.dataSize();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t idx = 0; idx < dataSize; ++idx) {
            const auto& e = dataStorage[idx];
            Scalar* cpuPtr = &((*e.matBlockAddress)[0][0]);

            // convert the pointer from CPU to GPU pointer based on offset in CPU jacobian
            Scalar* gpuPtr = ComputePtrBasedOnOffsetInOtherBuffer(
                gpuBufStart, gpuBufferSize, cpuBufStart, cpuBufferSize, cpuPtr);

            minimatrices[idx] = StructWithMinimatrix(
                e.neighbor, e.res_nbinfo, reinterpret_cast<MiniMatrixType*>(gpuPtr));
        }

        return SparseTable<StructWithMinimatrix, gpuistl::GpuBuffer>(
            gpuistl::GpuBuffer<StructWithMinimatrix>(minimatrices),
            gpuistl::GpuBuffer<int>(cpuNeighborInfoTable.rowStarts())
        );
    }

    // Handle the BoundaryInfo structs
    template <class GpuVecBlock,
              class GpuFluidState,
              class BoundaryInfoTypeGPU,
              class BoundaryInfoTypeCPU,
              typename GpuFluidSystemPtr>
    auto copy_to_gpu(const std::vector<BoundaryInfoTypeCPU>& cpu_boundary_info,
                     const GpuFluidSystemPtr& dynamicGpuFluidSystemPtr)
    {
        std::vector<BoundaryInfoTypeGPU> gpu_boundary_info;
        for (const auto& info : cpu_boundary_info) {
            gpu_boundary_info.push_back(BoundaryInfoTypeGPU {
                info.cell,
                info.dir,
                info.bfIndex,
                BoundaryConditionData<GpuVecBlock, GpuFluidState> {
                    info.bcdata.type,
                    GpuVecBlock(info.bcdata.massRate),
                    info.bcdata.pvtRegionIdx,
                    info.bcdata.boundaryFaceIndex,
                    info.bcdata.faceArea,
                    info.bcdata.faceZCoord,
                    info.bcdata.exFluidState.withOtherFluidSystem(dynamicGpuFluidSystemPtr)}});
        }

        return gpuistl::GpuBuffer<BoundaryInfoTypeGPU>(gpu_boundary_info);
    }

    // Implemented for residual_, which is a vector of FieldVectors.
    // We then make a GpuBuffer of MiniVectors.
    template <class CPUResidualType, class GpuMiniVector>
    auto copy_to_gpu_residual(CPUResidualType& residual)
    {
        std::vector<GpuMiniVector> vectorOfMiniVectors;
        for (const auto& minivec : residual) {
            vectorOfMiniVectors.emplace_back(minivec);
        }

        return gpuistl::GpuBuffer<GpuMiniVector>(vectorOfMiniVectors);
    }

} // namespace gpuistl
#endif

} // namespace Opm

#endif // TPFA_LINEARIZER_STRUCTS_HH
