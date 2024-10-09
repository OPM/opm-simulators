/*
  Copyright 2023 Equinor ASA

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

#include <config.h> // CMake

// MultisegmentWellContribution includes the cuda runtime if found by CMake
// this leads to inclusion of both amd_hip_vector_types.h and vector_types.h
// which both define vector types like uchar2, short3 and double4.
// Restore the value (if defined) afterwards.
#ifdef HAVE_CUDA
#define HIP_HAVE_CUDA_DEFINED HAVE_CUDA
#endif

#undef HAVE_CUDA

#include <opm/simulators/linalg/bda/rocm/rocsparseWellContributions.hpp>

#ifdef HIP_HAVE_CUDA_DEFINED
#define HAVE_CUDA HIP_HAVE_CUDA_DEFINED
#undef HIP_HAVE_CUDA_DEFINED
#endif

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/MultisegmentWellContribution.hpp>
#include <opm/simulators/linalg/bda/Misc.hpp>
#include <hip/hip_runtime.h>

namespace Opm
{

#ifdef __HIP__
/// HIP kernel to apply the standard wellcontributions
template<class Scalar>
__global__ void stdwell_apply(const Scalar* Cnnzs,
                              const Scalar* Dnnzs,
                              const Scalar* Bnnzs,
                              const unsigned* Ccols,
                              const unsigned* Bcols,
                              const Scalar* x,
                              Scalar* y,
                              const unsigned dim,
                              const unsigned dim_wells,
                              const unsigned *val_pointers)
{
    unsigned wgId = blockIdx.x;
    unsigned wiId = threadIdx.x;
    unsigned valSize = val_pointers[wgId + 1] - val_pointers[wgId];
    unsigned valsPerBlock = dim*dim_wells;
    unsigned numActiveWorkItems = (blockDim.x/valsPerBlock)*valsPerBlock;
    unsigned numBlocksPerWarp = blockDim.x/valsPerBlock;
    unsigned c = wiId % dim;
    unsigned r = (wiId/dim) % dim_wells;
    Scalar temp;

    extern __shared__ Scalar localSum[];
    Scalar* z1 = localSum + gridDim.x;
    Scalar* z2 = z1 + dim_wells;

    localSum[wiId] = 0;
    if (wiId < numActiveWorkItems) {
        unsigned b = wiId/valsPerBlock + val_pointers[wgId];
        while (b < valSize + val_pointers[wgId]) {
            int colIdx = Bcols[b];
            localSum[wiId] += Bnnzs[b*dim*dim_wells + r*dim + c]*x[colIdx*dim + c];
            b += numBlocksPerWarp;
        }

        // merge all blocks in this workgroup into 1 block
        // if numBlocksPerWarp >= 3, should use loop
        // block 1:     block 2:
        //  0  1  2     12 13 14
        //  3  4  5     15 16 17
        //  6  7  8     18 19 20
        //  9 10 11     21 22 23
        // workitem i will hold the sum of workitems i and i + valsPerBlock
        if (wiId < valsPerBlock){
            for (unsigned i = 1; i < numBlocksPerWarp; ++i) {
                localSum[wiId] += localSum[wiId + i*valsPerBlock];
	    }
        }

        if (c == 0 && wiId < valsPerBlock){
            for (unsigned i = dim - 1; i > 0; --i) {
                localSum[wiId] += localSum[wiId + i];
            }
            z1[r] = localSum[wiId];
        }
    }

    __syncthreads();

    if(wiId < dim_wells){
        temp = 0.0;
        for (unsigned i = 0; i < dim_wells; ++i) {
            temp += Dnnzs[wgId*dim_wells*dim_wells + wiId*dim_wells + i]*z1[i];
        }
        z2[wiId] = temp;
    }

    __syncthreads();

    if (wiId < dim*valSize){
        temp = 0.0;
        unsigned bb = wiId/dim + val_pointers[wgId];
        for (unsigned j = 0; j < dim_wells; ++j) {
            temp += Cnnzs[bb*dim*dim_wells + j*dim + c]*z2[j];
        }

        int colIdx = Ccols[bb];
        y[colIdx*dim + c] -= temp;
    }
}
#endif

template<class Scalar>
void WellContributionsRocsparse<Scalar>::
apply_stdwells([[maybe_unused]] Scalar* d_x,
               [[maybe_unused]] Scalar* d_y)
{
#ifdef __HIP__
    unsigned gridDim = this->num_std_wells;
    unsigned blockDim = 64;
    unsigned shared_mem_size = (blockDim + 2 * this->dim_wells) * sizeof(Scalar); // shared memory for localSum, z1 and z2
    // dim3(N) will create a vector {N, 1, 1}
    stdwell_apply<<<dim3(gridDim), dim3(blockDim), shared_mem_size, stream>>>(d_Cnnzs_hip,
                                                                              d_Dnnzs_hip,
                                                                              d_Bnnzs_hip,
                                                                              d_Ccols_hip,
                                                                              d_Bcols_hip,
                                                                              d_x,
                                                                              d_y,
                                                                              this->dim,
                                                                              this->dim_wells,
                                                                              d_val_pointers_hip
    );
    HIP_CHECK(hipStreamSynchronize(stream));
#else
    OPM_THROW(std::logic_error, "Error separate wellcontributions for rocsparse only supported when compiling with hipcc");
#endif
}

template<class Scalar>
void WellContributionsRocsparse<Scalar>::
apply_mswells(Scalar* d_x, Scalar* d_y)
{
    if (h_x.empty()) {
        h_x.resize(this->N);
        h_y.resize(this->N);
    }

    HIP_CHECK(hipMemcpyAsync(h_x.data(), d_x, sizeof(Scalar) * this->N, hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipMemcpyAsync(h_y.data(), d_y, sizeof(Scalar) * this->N, hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));

    // actually apply MultisegmentWells
    for (auto& well : this->multisegments) {
        well->apply(h_x.data(), h_y.data());
    }

    // copy vector y from CPU to GPU
    HIP_CHECK(hipMemcpyAsync(d_y, h_y.data(), sizeof(Scalar) * this->N, hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
}

template<class Scalar>
void WellContributionsRocsparse<Scalar>::
apply(Scalar* d_x, Scalar* d_y)
{
    if (this->num_std_wells > 0) {
        apply_stdwells(d_x, d_y);
    }

    if (this->num_ms_wells > 0) {
        apply_mswells(d_x, d_y);
    }
}

template<class Scalar>
void WellContributionsRocsparse<Scalar>::setStream(hipStream_t stream_)
{
    stream = stream_;
}

template<class Scalar>
void WellContributionsRocsparse<Scalar>::
APIaddMatrix(MatrixType type,
             int* colIndices,
             Scalar* values,
             unsigned int val_size)
{
    if (!this->allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }

    switch (type) {
    case MatrixType::C:
        HIP_CHECK(hipMemcpyAsync(d_Cnnzs_hip + this->num_blocks_so_far * this->dim * this->dim_wells,
                                 values, sizeof(d_Cnnzs_hip) * val_size * this->dim * this->dim_wells,
                                 hipMemcpyHostToDevice, stream));
        HIP_CHECK(hipMemcpyAsync(d_Ccols_hip + this->num_blocks_so_far, colIndices,
                                 sizeof(d_Ccols_hip) * val_size,
                                 hipMemcpyHostToDevice, stream));
        break;

    case MatrixType::D:
        HIP_CHECK(hipMemcpyAsync(d_Dnnzs_hip + this->num_std_wells_so_far * this->dim_wells * this->dim_wells,
                                 values, sizeof(d_Dnnzs_hip) * this->dim_wells * this->dim_wells,
                                 hipMemcpyHostToDevice, stream));
        break;

    case MatrixType::B:
        HIP_CHECK(hipMemcpyAsync(d_Bnnzs_hip + this->num_blocks_so_far * this->dim * this->dim_wells,
                                 values, sizeof(d_Bnnzs_hip) * val_size * this->dim * this->dim_wells,
                                 hipMemcpyHostToDevice, stream));
        HIP_CHECK(hipMemcpyAsync(d_Bcols_hip + this->num_blocks_so_far, colIndices,
                                 sizeof(d_Bcols_hip) * val_size,
                                 hipMemcpyHostToDevice, stream));

        this->val_pointers[this->num_std_wells_so_far] = this->num_blocks_so_far;
        if (this->num_std_wells_so_far == this->num_std_wells - 1) {
            this->val_pointers[this->num_std_wells] = this->num_blocks;
            HIP_CHECK(hipMemcpyAsync(d_val_pointers_hip, this->val_pointers.data(),
                                     sizeof(d_val_pointers_hip) * (this->num_std_wells + 1),
                                     hipMemcpyHostToDevice, stream));
        }
        break;

    default:
        OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributionsRocsparse::addMatrix()");
    }
    HIP_CHECK(hipStreamSynchronize(stream));
}

template<class Scalar>
void WellContributionsRocsparse<Scalar>::APIalloc()
{
    HIP_CHECK(hipMalloc((void**)&d_Cnnzs_hip,
                        sizeof(d_Cnnzs_hip) * this->num_blocks * this->dim * this->dim_wells));
    HIP_CHECK(hipMalloc((void**)&d_Dnnzs_hip,
                        sizeof(d_Dnnzs_hip) * this->num_std_wells * this->dim_wells * this->dim_wells));
    HIP_CHECK(hipMalloc((void**)&d_Bnnzs_hip,
                        sizeof(d_Bnnzs_hip) * this->num_blocks * this->dim * this->dim_wells));
    HIP_CHECK(hipMalloc((void**)&d_Ccols_hip, sizeof(d_Ccols_hip) * this->num_blocks));
    HIP_CHECK(hipMalloc((void**)&d_Bcols_hip, sizeof(d_Bcols_hip) * this->num_blocks));
    HIP_CHECK(hipMalloc((void**)&d_val_pointers_hip,
                        sizeof(d_val_pointers_hip) * (this->num_std_wells + 1)));
}

template class WellContributionsRocsparse<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellContributionsRocsparse<float>;
#endif

} // namespace Opm
