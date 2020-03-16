/*
  Copyright 2020 Equinor ASA

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


#include <cstdlib>
#include <cstring>
#include <config.h> // CMake

#ifdef __NVCC__
#include "opm/simulators/linalg/bda/cuda_header.hpp"
#include <cuda_runtime.h>
#endif

#include "opm/simulators/linalg/bda/WellContributions.hpp"

namespace Opm
{

    // apply WellContributions using y -= C^T * (D^-1 * (B * x))
#if HAVE_CUDA
    __global__ void apply_well_contributions( \
        const double * __restrict__ Cnnzs, \
        const double * __restrict__ Dnnzs, \
        const double * __restrict__ Bnnzs, \
        const int * __restrict__ Ccols, \
        const int * __restrict__ Bcols, \
        const double * __restrict__ x, \
        double * __restrict__ y, \
        const int dim, \
        const int dim_wells, \
        const unsigned int * __restrict__ val_pointers \
        )
    {
        const int idx_b = blockIdx.x;
        const int idx_t = threadIdx.x;
        int idx = idx_b * blockDim.x + idx_t;
        const unsigned int val_size = val_pointers[idx_b+1] - val_pointers[idx_b];

        const int vals_per_block = dim * dim_wells;        // 12
        const int num_active_threads = (32/vals_per_block)*vals_per_block; // 24
        const int num_blocks_per_warp = 32/vals_per_block; // 2
        const int lane = idx_t % 32;
        const int c = lane % dim;                           // col in block
        const int r = (lane / dim) % dim_wells;             // row in block
        const int NUM_THREADS = gridDim.x * blockDim.x;

        extern __shared__ double smem[];
        double * __restrict__ z1 = smem;
        double * __restrict__ z2 = z1 + dim_wells;

        if(idx_t < dim_wells){
            z1[idx_t] = 0.0;
        }

        __syncthreads();

        // z1 = B * x
        if(idx_t < num_active_threads){
            // multiply all blocks with x
            double temp = 0.0;
            int b = idx_t / vals_per_block + val_pointers[idx_b];       // block id, val_size indicates number of blocks
            while(b < val_size + val_pointers[idx_b]){
                int colIdx = Bcols[b];
                temp += Bnnzs[b * dim * dim_wells + r * dim + c] * x[colIdx * dim + c];
                b += num_blocks_per_warp;
            }

            // merge all blocks into 1 dim*dim_wells block
            // since NORNE has only 2 parallel blocks, do not use a loop
            temp += __shfl_down_sync(0x00ffffff, temp, dim * dim_wells);

            b = idx_t / vals_per_block + val_pointers[idx_b];

            // merge all (dim) columns of 1 block, results in a single 1*dim_wells vector, which is used to multiply with invD
            if(idx_t < vals_per_block){
                // should be a loop as well, now only works for dim == 3
                if(c == 0 || c == 2){temp += __shfl_down_sync(0x00000B6D, temp, 2);}  // add col 2 to col 0
                if(c == 0 || c == 1){temp += __shfl_down_sync(0x000006DB, temp, 1);}  // add col 1 to col 0
            }

            // write 1*dim_wells vector to gmem, could be replaced with shfl broadcast to remove z1 altogether
            if(c == 0 && idx_t < vals_per_block){
                z1[r] = temp;
            }
        }

        __syncthreads();

        // z2 = D^-1 * B * x = D^-1 * z1
        if(idx_t < dim_wells){
            double temp = 0.0;
            for(int c = 0; c < dim_wells; ++c){
                temp += Dnnzs[idx_b * dim_wells * dim_wells + idx_t * dim_wells + c] * z1[c];
            }
            z2[idx_t] = temp;
        }

        __syncthreads();

        // y -= C^T * D^-1 * B * x
        // use dim * val_size threads, each block is assigned 'dim' threads
        if(idx_t < dim * val_size){
            double temp = 0.0;
            int b = idx_t / dim + val_pointers[idx_b];
            int cc = idx_t % dim;
            int colIdx = Ccols[b];
            for(unsigned int c = 0; c < dim_wells; ++c){
                temp += Cnnzs[b * dim * dim_wells + c * dim + cc] * z2[c];
            }
            y[colIdx * dim + cc] -= temp;
        }

    }
#endif

        void WellContributions::alloc_all(){
#if HAVE_CUDA
            if(gpu_mode){
                alloc_gpu();
            }else{
                alloc_cpu();
            }
#else
            alloc_cpu();
#endif
            allocated = true;
        }

        void WellContributions::alloc_cpu(){
            Cnnzs = new double[num_blocks * dim * dim_wells];
            Dnnzs = new double[num_wells * dim_wells * dim_wells];
            Bnnzs = new double[num_blocks * dim * dim_wells];
            Ccols = new int[num_blocks];
            Bcols = new int[num_blocks];
            val_pointers = new unsigned int[num_wells + 1];
            z1 = new double[dim_wells];    //        B * x
            z2 = new double[dim_wells];    // D^-1 * B * x
        }

#if HAVE_CUDA
        void WellContributions::alloc_gpu(){
            cudaMalloc((void**)&d_Cnnzs, sizeof(double) * num_blocks * dim * dim_wells);
            cudaMalloc((void**)&d_Dnnzs, sizeof(double) * num_wells * dim_wells * dim_wells);
            cudaMalloc((void**)&d_Bnnzs, sizeof(double) * num_blocks * dim * dim_wells);
            cudaMalloc((void**)&d_Ccols, sizeof(int) * num_blocks);
            cudaMalloc((void**)&d_Bcols, sizeof(int) * num_blocks);
            val_pointers = new unsigned int[num_wells + 1];
            cudaMalloc((void**)&d_val_pointers, sizeof(int) * (num_wells + 1));
            cudaCheckLastError("apply_gpu malloc failed");
        }
#endif

        WellContributions::~WellContributions()
        {
#if HAVE_CUDA
            if(gpu_mode){
                free_gpu();
            }else{
                free_cpu();
            }
#else
            free_cpu();
#endif
        }

        void WellContributions::free_cpu(){
            delete[] Cnnzs;
            delete[] Dnnzs;
            delete[] Bnnzs;
            delete[] Ccols;
            delete[] Bcols;
            delete[] val_pointers;
            delete[] z1;
            delete[] z2;
            //delete[] Mnnzs;
        }

#if HAVE_CUDA
        void WellContributions::free_gpu(){
            cudaFree(d_Cnnzs);
            cudaFree(d_Dnnzs);
            cudaFree(d_Bnnzs);
            cudaFree(d_Ccols);
            cudaFree(d_Bcols);
            delete[] val_pointers;
            cudaFree(d_val_pointers);
            // cudaFree(d_z1);
            // cudaFree(d_z2);
        }
#endif


        void WellContributions::apply(double *x, double *y){
#if HAVE_CUDA
            if (gpu_mode){
                apply_gpu(x, y);
            }else{
                apply_cpu(x, y);
            }
#else
            apply_cpu(x, y);
#endif
        }

        // Apply the WellContributions, similar to StandardWell::apply()
        // y -= (C^T *(D^-1*(   B*x)))
        void WellContributions::apply_cpu(double *x, double *y)
        {
#if 0
            // Mnnzs contains a sparse matrix with a symmetric pattern
            // Mrows would contain 'i*val_size' for every entry i, since every row has the same number of blocks
            // Mcols are the same as Ccols, normally, there is an entry for every block, but since all rows have the same sparsity pattern, we only have to store 1 row
            bool dbg = false;
            for(int i = 0; i < dim*dim*val_size*val_size; ++i){
                if(dbg)printf("Mnnzs[%d]: %.5e\n", i, Mnnzs[i]);
            }
            if(dbg)printf("row_size: %u, val_size: %u\n", row_size, val_size);
            for(int r = 0; r < val_size; ++r){
                for(int c = 0; c < val_size; ++c){
                    int colIdx = Ccols[c];
                    if(dbg)printf("colIdx: %d\n", colIdx);
                    for(int i = 0; i < dim; ++i){
                        double sum = 0.0;
                        for(int j = 0; j < dim; ++j){
                            sum += Mnnzs[r * dim * dim * val_size + c * dim * dim + i * dim + j] * x[colIdx * dim + j];
                        }
                        if(dbg)printf("sum: %f\n", sum);
                        y[colIdx * dim + i] -= sum;
                    }
                }
            }
            if(dbg)exit(0);
#else
            for(int wellID = 0; wellID < num_wells; ++wellID){
                unsigned int val_size = val_pointers[wellID+1] - val_pointers[wellID];

                // B * x
                for (unsigned int i = 0; i < dim_wells; ++i) {
                    z1[i] = 0.0;
                }
                for (unsigned int i = 0; i < val_size; ++i) {
                    unsigned int blockID = i + val_pointers[wellID];
                    int colIdx = Bcols[blockID];
                    for (unsigned int j = 0; j < dim_wells; ++j) {
                        double temp = 0.0;
                        for (unsigned int k = 0; k < dim; ++k) {
                            temp += Bnnzs[blockID * dim * dim_wells + j * dim + k] * x[colIdx * dim + k];
                        }
                        z1[j] += temp;
                    }
                }

                // D^-1 * B * x
                for (unsigned int i = 0; i < dim_wells; ++i) {
                    z2[i] = 0.0;
                }
                for (unsigned int j = 0; j < dim_wells; ++j) {
                    double temp = 0.0;
                    for (unsigned int k = 0; k < dim_wells; ++k) {
                        temp += Dnnzs[wellID * dim_wells * dim_wells + j * dim_wells + k] * z1[k];
                    }
                    z2[j] += temp;
                }

                // C^T * D^-1 * B * x
                for (unsigned int i = 0; i < val_size; ++i) {
                    unsigned int blockID = i + val_pointers[wellID];
                    int colIdx = Ccols[blockID];
                    for (unsigned int j = 0; j < dim; ++j) {
                        double temp = 0.0;
                        for (unsigned int k = 0; k < dim_wells; ++k) {
                            temp += Cnnzs[blockID * dim * dim_wells + j + k * dim] * z2[k];
                        }
                        y[colIdx * dim + j] -= temp;
                    }
                }
            }
#endif
        }


        // Apply the WellContributions, similar to StandardWell::apply()
        // y -= (C^T *(D^-1*(   B*x)))
#if HAVE_CUDA
        void WellContributions::apply_gpu(double *d_x, double *d_y)
        {
            int smem_size = 2 * sizeof(double) * dim_wells;
            apply_well_contributions<<<num_wells, 32, smem_size, stream>>>(d_Cnnzs, d_Dnnzs, d_Bnnzs, d_Ccols, d_Bcols, d_x, d_y, dim, dim_wells, d_val_pointers);
        }
#endif

        void WellContributions::addMatrix(int idx, int *colIndices, double *values, unsigned int val_size)
        {
#if HAVE_CUDA
            if(gpu_mode){
                addMatrix_gpu(idx, colIndices, values, val_size);
            }else{
                addMatrix_cpu(idx, colIndices, values, val_size);
            }
#else
            addMatrix_cpu(idx, colIndices, values, val_size);
#endif
            if(idx == 2){
                num_blocks_so_far += val_size;
            }
            if(idx == 2){
                num_wells_so_far++;
            }
        }


        void WellContributions::addMatrix_cpu(int idx, int *colIndices, double *values, unsigned int val_size)
        {
            switch (idx) {
            case 0:
                memcpy(Cnnzs + num_blocks_so_far * dim * dim_wells, values, sizeof(double) * val_size * dim * dim_wells);
                memcpy(Ccols + num_blocks_so_far, colIndices, sizeof(int) * val_size);
                break;
            case 1:
                memcpy(Dnnzs + num_wells_so_far * dim_wells * dim_wells, values, sizeof(double) * dim_wells * dim_wells);
                break;
            case 2:
                memcpy(Bnnzs + num_blocks_so_far * dim * dim_wells, values, sizeof(double) * val_size * dim * dim_wells);
                memcpy(Bcols + num_blocks_so_far, colIndices, sizeof(int) * val_size);
                val_pointers[num_wells_so_far] = num_blocks_so_far;
                if(num_wells_so_far == num_wells - 1){
                    val_pointers[num_wells] = num_blocks;
                }
                break;
            case 3:
                // store (C*D*B)
                printf("ERROR unsupported matrix ID for WellContributions::addMatrix()\n");
                exit(1);
                // memcpy(Mnnzs, values, sizeof(double) * dim * dim * val_size * val_size);
                // memcpy(Ccols, colIndices, sizeof(int) * val_size);
                break;
            default:
                printf("ERROR unknown matrix ID for WellContributions::addMatrix()\n");
                exit(1);
            }
        }


#if HAVE_CUDA
        void WellContributions::addMatrix_gpu(int idx, int *colIndices, double *values, unsigned int val_size)
        {
            switch (idx) {
            case 0:
                cudaMemcpy(d_Cnnzs + num_blocks_so_far * dim * dim_wells, values, sizeof(double) * val_size * dim * dim_wells, cudaMemcpyHostToDevice);
                cudaMemcpy(d_Ccols + num_blocks_so_far, colIndices, sizeof(int) * val_size, cudaMemcpyHostToDevice);
                break;
            case 1:
                cudaMemcpy(d_Dnnzs + num_wells_so_far * dim_wells * dim_wells, values, sizeof(double) * dim_wells * dim_wells, cudaMemcpyHostToDevice);
                break;
            case 2:
                cudaMemcpy(d_Bnnzs + num_blocks_so_far * dim * dim_wells, values, sizeof(double) * val_size * dim * dim_wells, cudaMemcpyHostToDevice);
                cudaMemcpy(d_Bcols + num_blocks_so_far, colIndices, sizeof(int) * val_size, cudaMemcpyHostToDevice);
                val_pointers[num_wells_so_far] = num_blocks_so_far;
                if(num_wells_so_far == num_wells - 1){
                    val_pointers[num_wells] = num_blocks;
                }
                cudaMemcpy(d_val_pointers, val_pointers, sizeof(int) * (num_wells+1), cudaMemcpyHostToDevice);
                break;
            case 3:
                // store (C*D*B)
                printf("ERROR unsupported matrix ID for WellContributions::addMatrix()\n");
                exit(1);
                break;
            default:
                printf("ERROR unknown matrix ID for WellContributions::addMatrix()\n");
                exit(1);
            }
            cudaCheckLastError("WellContributions::addMatrix() failed");
        }

        void WellContributions::setCudaStream(cudaStream_t stream_)
        {
            this->stream = stream_;
        }

#endif

        void WellContributions::addSizes(unsigned int nnz, unsigned int numEq, unsigned int numWellEq)
        {
            if(allocated){
                std::cerr << "Error cannot add more sizes after allocated in WellContributions" << std::endl;
                exit(1);
            }
            num_blocks += nnz;
            dim = numEq;
            dim_wells = numWellEq;
            num_wells++;
        }

        // Default value
        bool WellContributions::gpu_mode = false;

        // If HAVE_CUDA is false, use_gpu must be false too
        void WellContributions::setMode(bool use_gpu){
            gpu_mode = use_gpu;
        }

} //namespace Opm

