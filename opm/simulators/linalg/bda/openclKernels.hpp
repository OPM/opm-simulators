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

#ifndef OPENCL_HPP
#define OPENCL_HPP

namespace bda
{

    inline const char* axpy_s  = R"(
    __kernel void axpy(
        __global double *in,
        const double a,
        __global double *out,
        const int N)
    {
        unsigned int NUM_THREADS = get_global_size(0);
        int idx = get_global_id(0);

        while(idx < N){
            out[idx] += a * in[idx];
            idx += NUM_THREADS;
        }
    }
    )";


    // returns partial sums, instead of the final dot product
    inline const char* dot_1_s = R"(
    __kernel void dot_1(
        __global double *in1,
        __global double *in2,
        __global double *out,
        const unsigned int N,
        __local double *tmp)
    {
        unsigned int tid = get_local_id(0);
        unsigned int bsize = get_local_size(0);
        unsigned int bid = get_global_id(0) / bsize;
        unsigned int i = get_global_id(0);
        unsigned int NUM_THREADS = get_global_size(0);

        double sum = 0.0;
        while(i < N){
            sum += in1[i] * in2[i];
            i += NUM_THREADS;
        }
        tmp[tid] = sum;

        barrier(CLK_LOCAL_MEM_FENCE);

        // do reduction in shared mem
        for(unsigned int s = get_local_size(0) / 2; s > 0; s >>= 1)
        {
            if (tid < s)
            {
                tmp[tid] += tmp[tid + s];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // write result for this block to global mem
        if (tid == 0) out[get_group_id(0)] = tmp[0];
    }
    )";


    // returns partial sums, instead of the final norm
    // the square root must be computed on CPU
    inline const char* norm_s = R"(
    __kernel void norm(
        __global double *in,
        __global double *out,
        const unsigned int N,
        __local double *tmp)
    {
        unsigned int tid = get_local_id(0);
        unsigned int bsize = get_local_size(0);
        unsigned int bid = get_global_id(0) / bsize;
        unsigned int i = get_global_id(0);
        unsigned int NUM_THREADS = get_global_size(0);

        double local_sum = 0.0;
        while(i < N){
            local_sum += in[i] * in[i];
            i += NUM_THREADS;
        }
        tmp[tid] = local_sum;

        barrier(CLK_LOCAL_MEM_FENCE);

        // do reduction in shared mem
        for(unsigned int s = get_local_size(0) / 2; s > 0; s >>= 1)
        {
            if (tid < s)
            {
                tmp[tid] += tmp[tid + s];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // write result for this block to global mem
        if (tid == 0) out[get_group_id(0)] = tmp[0];
    }
    )";


    // p = (p - omega * v) * beta + r
    inline const char* custom_s = R"(
    __kernel void custom(
        __global double *p,
        __global double *v,
        __global double *r,
        const double omega,
        const double beta,
        const int N)
    {
        const unsigned int NUM_THREADS = get_global_size(0);
        unsigned int idx = get_global_id(0);

        while(idx < N){
            double res = p[idx];
            res -= omega * v[idx];
            res *= beta;
            res += r[idx];
            p[idx] = res;
            idx += NUM_THREADS;
        }
    }
    )";


    // b = mat * x
    // algorithm based on:
    // Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
    // Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
    inline const char* spmv_blocked_s = R"(
    __kernel void spmv_blocked(
        __global const double *vals,
        __global const int *cols,
        __global const int *rows,
        const int Nb,
        __global const double *x,
        __global double *b,
        const unsigned int block_size,
        __local double *tmp)
    {
        const unsigned int warpsize = 32;
        const unsigned int bsize = get_local_size(0);
        const unsigned int idx_b = get_global_id(0) / bsize;
        const unsigned int idx_t = get_local_id(0);
        unsigned int idx = idx_b * bsize + idx_t;
        const unsigned int bs = block_size;
        const unsigned int num_active_threads = (warpsize/bs/bs)*bs*bs;
        const unsigned int num_blocks_per_warp = warpsize/bs/bs;
        const unsigned int NUM_THREADS = get_global_size(0);
        const unsigned int num_warps_in_grid = NUM_THREADS / warpsize;
        unsigned int target_block_row = idx / warpsize;
        const unsigned int lane = idx_t % warpsize;
        const unsigned int c = (lane / bs) % bs;
        const unsigned int r = lane % bs;

        // for 3x3 blocks:
        // num_active_threads: 27
        // num_blocks_per_warp: 3

        while(target_block_row < Nb){
            unsigned int first_block = rows[target_block_row];
            unsigned int last_block = rows[target_block_row+1];
            unsigned int block = first_block + lane / (bs*bs);
            double local_out = 0.0;

            if(lane < num_active_threads){
                for(; block < last_block; block += num_blocks_per_warp){
                    double x_elem = x[cols[block]*bs + c];
                    double A_elem = vals[block*bs*bs + c + r*bs];
                    local_out += x_elem * A_elem;
                }
            }

            // do reduction in shared mem
            tmp[lane] = local_out;
            barrier(CLK_LOCAL_MEM_FENCE);

            for(unsigned int offset = 3; offset <= 24; offset <<= 1)
            {
                if (lane + offset < warpsize)
                {
                    tmp[lane] += tmp[lane + offset];
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }

            if(lane < bs){
                unsigned int row = target_block_row*bs + lane;
                b[row] = tmp[lane];
            }
            target_block_row += num_warps_in_grid;
        }
    }
    )";



    // ILU apply part 1: forward substitution
    // solves L*x=y where L is a lower triangular sparse blocked matrix
    inline const char* ILU_apply1_s = R"(
    __kernel void ILU_apply1(
        __global const double *Lvals,
        __global const unsigned int *Lcols,
        __global const unsigned int *Lrows,
        const unsigned int LrowSize,
        __global const double *y,
        __global double *x,
        __global const unsigned int *nodesPerColorPrefix,
        const unsigned int color,
        const unsigned int block_size,
        __local double *tmp)
    {
        const unsigned int warpsize = 32;
        const unsigned int bs = block_size;
        const unsigned int idx_t = get_local_id(0);
        const unsigned int num_active_threads = (warpsize/bs/bs)*bs*bs;
        const unsigned int num_blocks_per_warp = warpsize/bs/bs;
        const unsigned int NUM_THREADS = get_global_size(0);
        const unsigned int num_warps_in_grid = NUM_THREADS / warpsize;
        unsigned int idx = get_global_id(0);
        unsigned int target_block_row = idx / warpsize;
        idx += nodesPerColorPrefix[color];
        target_block_row += nodesPerColorPrefix[color];
        const unsigned int lane = idx_t % warpsize;
        const unsigned int c = (lane / bs) % bs;
        const unsigned int r = lane % bs;

        while(target_block_row < nodesPerColorPrefix[color+1]){
            const unsigned int first_block = Lrows[target_block_row];
            const unsigned int last_block = Lrows[target_block_row+1];
            unsigned int block = first_block + lane / (bs*bs);
            double local_out = 0.0;
            if(lane < num_active_threads){
                if(lane < bs){
                    local_out = y[target_block_row*bs+lane];
                }
                for(; block < last_block; block += num_blocks_per_warp){
                    const double x_elem = x[Lcols[block]*bs + c];
                    const double A_elem = Lvals[block*bs*bs + c + r*bs];
                    local_out -= x_elem * A_elem;
                }
            }

            // do reduction in shared mem
            tmp[lane] = local_out;
            barrier(CLK_LOCAL_MEM_FENCE);

            for(unsigned int offset = 3; offset <= 24; offset <<= 1)
            {
                if (lane + offset < warpsize)
                {
                    tmp[lane] += tmp[lane + offset];
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }

            if(lane < bs){
                const unsigned int row = target_block_row*bs + lane;
                x[row] = tmp[lane];
            }

            target_block_row += num_warps_in_grid;
        }
    }
    )";


    // ILU apply part 2: backward substitution
    // solves U*x=y where L is a lower triangular sparse blocked matrix
    inline const char* ILU_apply2_s = R"(
    __kernel void ILU_apply2(
        __global const double *Uvals,
        __global const int *Ucols,
        __global const int *Urows,
        const unsigned int UrowSize,
        __global const double *invDiagVals,
        __global double *x,
        __global const unsigned int *nodesPerColorPrefix,
        const unsigned int color,
        const unsigned int block_size,
        __local double *tmp)
    {
        const unsigned int warpsize = 32;
        const unsigned int bs = block_size;
        const unsigned int idx_t = get_local_id(0);
        const unsigned int num_active_threads = (warpsize/bs/bs)*bs*bs;
        const unsigned int num_blocks_per_warp = warpsize/bs/bs;
        const unsigned int NUM_THREADS = get_global_size(0);
        const unsigned int num_warps_in_grid = NUM_THREADS / warpsize;
        unsigned int idx = get_global_id(0);
        unsigned int target_block_row = idx / warpsize;
        idx += nodesPerColorPrefix[color];
        target_block_row += nodesPerColorPrefix[color];
        const unsigned int lane = idx_t % warpsize;
        const unsigned int c = (lane / bs) % bs;
        const unsigned int r = lane % bs;
        const double relaxation = 0.9;

        while(target_block_row < nodesPerColorPrefix[color+1]){
            const unsigned int first_block = Urows[UrowSize-target_block_row-1];
            const unsigned int last_block = Urows[UrowSize-target_block_row];
            unsigned int block = first_block + lane / (bs*bs);
            double local_out = 0.0;
            if(lane < num_active_threads){
                if(lane < bs){
                    const unsigned int row = target_block_row*bs+lane;
                    local_out = x[row];
                }
                for(; block < last_block; block += num_blocks_per_warp){
                    const double x_elem = x[Ucols[block]*bs + c];
                    const double A_elem = Uvals[block*bs*bs + c + r*bs];
                    local_out -= x_elem * A_elem;
                }
            }

            // do reduction in shared mem
            tmp[lane] = local_out;
            barrier(CLK_LOCAL_MEM_FENCE);

            for(unsigned int offset = 3; offset <= 24; offset <<= 1)
            {
                if (lane + offset < warpsize)
                {
                    tmp[lane] += tmp[lane + offset];
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
            local_out = tmp[lane];

            if(lane < bs){
                tmp[lane + bs*idx_t/warpsize] = local_out;
                double sum = 0.0;
                for(int i = 0; i < bs; ++i){
                    sum += invDiagVals[target_block_row*bs*bs + i + lane*bs] * tmp[i + bs*idx_t/warpsize];
                }

                const unsigned int row = target_block_row*bs + lane;
                x[row] = relaxation * sum;
            }

            target_block_row += num_warps_in_grid;
        }
    }
    )";

    inline const char* add_well_contributions_s = R"(
    #pragma OPENCL EXTENSION cl_khr_fp64: enable
    #pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable

    void atomicAdd(volatile __global double *val, const double delta){
        union{
            double f;
            ulong i;
        } old;

        union{
            double f;
            ulong i;
        } new;

        do{
            old.f = *val;
            new.f = old.f + delta;
        } while(atom_cmpxchg((volatile __global ulong *)val, old.i, new.i) != old.i);
    }

    __kernel void add_well_contributions(__global const double *valsC,
                                         __global const double *valsD,
                                         __global const double *valsB,
                                         __global const int *colsC,
                                         __global const int *colsB,
                                         __global const double *x,
                                         __global double *y,
                                         const unsigned int blnc,
                                         const unsigned int blnr,
                                         __global const unsigned int *rowptr,
                                         __local double *localSum,
                                         __local double *z1,
                                         __local double *z2){
        int wgId = get_group_id(0);
        int wiId = get_local_id(0);
        int valSize = rowptr[wgId + 1] - rowptr[wgId];
        int valsPerBlock = blnc*blnr;
        int numActiveWorkItems = (32/valsPerBlock)*valsPerBlock;
        int numBlocksPerWarp = 32/valsPerBlock;
        int c = wiId % blnc;
        int r = (wiId/blnc) % blnr;
        double temp;

        localSum[wiId] = 0;
        if(wiId < numActiveWorkItems){
            int b = wiId/valsPerBlock + rowptr[wgId];
            while(b < valSize + rowptr[wgId]){
                int colIdx = colsB[b];
                localSum[wiId] += valsB[b*blnc*blnr + r*blnc + c]*x[colIdx*blnc + c];
                b += numBlocksPerWarp;
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        int stride = valsPerBlock;
        if(wiId < stride){
            localSum[wiId] += localSum[wiId + stride];
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        if(c == 0 && wiId < valsPerBlock){
            for(stride = 2; stride > 0; stride /= 2){
                localSum[wiId] += localSum[wiId + stride];
            }
            z1[r] = localSum[wiId];
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        if(wiId < blnr){
            temp = 0.0;
            for(unsigned int i = 0; i < blnr; ++i){
                temp += valsD[wgId*blnr*blnr + wiId*blnr + i]*z1[i];
            }
            z2[wiId] = temp;
        }

        barrier(CLK_GLOBAL_MEM_FENCE);

        if(wiId < blnc*valSize){
            temp = 0.0;
            int bb = wiId/blnc + rowptr[wgId];
            int colIdx = colsC[bb];
            for (unsigned int j = 0; j < blnr; ++j){
                temp += valsC[bb*blnc*blnr + j*blnc + c]*z2[j];
            }
            atomicAdd(&y[colIdx*blnc + c], -temp);
        }
    }
    )";

} // end namespace bda

#endif
