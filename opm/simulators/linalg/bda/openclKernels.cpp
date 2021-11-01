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

#include <config.h>
#include <opm/simulators/linalg/bda/openclKernels.hpp>

namespace bda
{

    std::string get_axpy_string() {
        return R"(
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
    }


    // scale vector with scalar
    std::string get_scale_string() {
        return R"(
        __kernel void scale(
            __global double *vec,
            const double a,
            const int N)
        {
            unsigned int NUM_THREADS = get_global_size(0);
            int idx = get_global_id(0);

            while(idx < N){
                vec[idx] *= a;
                idx += NUM_THREADS;
            }
        }
        )";
    }


    // returns partial sums, instead of the final dot product
    std::string get_dot_1_string() {
        return R"(
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
    }


    // returns partial sums, instead of the final norm
    // the square root must be computed on CPU
    std::string get_norm_string() {
        return R"(
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
    }


    // p = (p - omega * v) * beta + r
    std::string get_custom_string() {
        return R"(
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
    }

    std::string get_spmv_blocked_string() {
        return R"(
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
    }



    std::string get_ILU_apply1_string(bool full_matrix) {
        std::string s = R"(
            __kernel void ILU_apply1(
                __global const double *LUvals,
                __global const unsigned int *LUcols,
                __global const unsigned int *LUrows,
                __global const int *diagIndex,
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
                target_block_row += nodesPerColorPrefix[color];
                const unsigned int lane = idx_t % warpsize;
                const unsigned int c = (lane / bs) % bs;
                const unsigned int r = lane % bs;

                while(target_block_row < nodesPerColorPrefix[color+1]){
                    const unsigned int first_block = LUrows[target_block_row];
                    )";
        if (full_matrix) {
            s += "const unsigned int last_block = diagIndex[target_block_row];  ";
        } else {
            s += "const unsigned int last_block = LUrows[target_block_row+1];  ";
        }
        s += R"(
                    unsigned int block = first_block + lane / (bs*bs);
                    double local_out = 0.0;
                    if(lane < num_active_threads){
                        if(lane < bs){
                            local_out = y[target_block_row*bs+lane];
                        }
                        for(; block < last_block; block += num_blocks_per_warp){
                            const double x_elem = x[LUcols[block]*bs + c];
                            const double A_elem = LUvals[block*bs*bs + c + r*bs];
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
        return s;
    }


    std::string get_ILU_apply2_string(bool full_matrix) {
        std::string s = R"(
            __kernel void ILU_apply2(
                __global const double *LUvals,
                __global const int *LUcols,
                __global const int *LUrows,
                __global const int *diagIndex,
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
                unsigned int idx_g = get_global_id(0);
                unsigned int target_block_row = idx_g / warpsize;
                target_block_row += nodesPerColorPrefix[color];
                const unsigned int lane = idx_t % warpsize;
                const unsigned int c = (lane / bs) % bs;
                const unsigned int r = lane % bs;

                while(target_block_row < nodesPerColorPrefix[color+1]){
                    )";
        if (full_matrix) {
            s +=   "const unsigned int first_block = diagIndex[target_block_row] + 1;  ";
        } else {
            s +=   "const unsigned int first_block = LUrows[target_block_row];  ";
        }
        s += R"(
                    const unsigned int last_block = LUrows[target_block_row+1];
                    unsigned int block = first_block + lane / (bs*bs);
                    double local_out = 0.0;
                    if(lane < num_active_threads){
                        if(lane < bs){
                            const unsigned int row = target_block_row*bs+lane;
                            local_out = x[row];
                        }
                        for(; block < last_block; block += num_blocks_per_warp){
                            const double x_elem = x[LUcols[block]*bs + c];
                            const double A_elem = LUvals[block*bs*bs + c + r*bs];
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
                        x[row] = sum;
                    }

                    target_block_row += num_warps_in_grid;
                }
            }
        )";
        return s;
    }

    std::string get_stdwell_apply_string(bool reorder) {
        std::string kernel_name = reorder ? "stdwell_apply" : "stdwell_apply_no_reorder";
        std::string s = "__kernel void " + kernel_name + R"((
                        __global const double *Cnnzs,
                        __global const double *Dnnzs,
                        __global const double *Bnnzs,
                        __global const int *Ccols,
                        __global const int *Bcols,
                        __global const double *x,
                        __global double *y,
                        )";
        if (reorder) {
            s +=     R"(__global const int *toOrder,
                        )";
        }
        s +=         R"(const unsigned int dim,
                        const unsigned int dim_wells,
                        __global const unsigned int *val_pointers,
                        __local double *localSum,
                        __local double *z1,
                        __local double *z2){
                int wgId = get_group_id(0);
                int wiId = get_local_id(0);
                int valSize = val_pointers[wgId + 1] - val_pointers[wgId];
                int valsPerBlock = dim*dim_wells;
                int numActiveWorkItems = (get_local_size(0)/valsPerBlock)*valsPerBlock;
                int numBlocksPerWarp = get_local_size(0)/valsPerBlock;
                int c = wiId % dim;
                int r = (wiId/dim) % dim_wells;
                double temp;

                barrier(CLK_LOCAL_MEM_FENCE);

                localSum[wiId] = 0;
                if(wiId < numActiveWorkItems){
                    int b = wiId/valsPerBlock + val_pointers[wgId];
                    while(b < valSize + val_pointers[wgId]){
                        )";
        if (reorder) {
            s +=       "int colIdx = toOrder[Bcols[b]];  ";
        } else {
            s +=       "int colIdx = Bcols[b];  ";
        }
        s += R"(
                        localSum[wiId] += Bnnzs[b*dim*dim_wells + r*dim + c]*x[colIdx*dim + c];
                        b += numBlocksPerWarp;
                    }

                    if(wiId < valsPerBlock){
                        localSum[wiId] += localSum[wiId + valsPerBlock];
                    }

                    b = wiId/valsPerBlock + val_pointers[wgId];

                    if(c == 0 && wiId < valsPerBlock){
                        for(unsigned int stride = 2; stride > 0; stride >>= 1){
                            localSum[wiId] += localSum[wiId + stride];
                        }
                        z1[r] = localSum[wiId];
                    }
                }

                barrier(CLK_LOCAL_MEM_FENCE);

                if(wiId < dim_wells){
                    temp = 0.0;
                    for(unsigned int i = 0; i < dim_wells; ++i){
                        temp += Dnnzs[wgId*dim_wells*dim_wells + wiId*dim_wells + i]*z1[i];
                    }
                    z2[wiId] = temp;
                }

                barrier(CLK_LOCAL_MEM_FENCE);

                if(wiId < dim*valSize){
                    temp = 0.0;
                    int bb = wiId/dim + val_pointers[wgId];
                    for (unsigned int j = 0; j < dim_wells; ++j){
                        temp += Cnnzs[bb*dim*dim_wells + j*dim + c]*z2[j];
                    }
                    )";
        if (reorder) {
            s +=   "int colIdx = toOrder[Ccols[bb]];  ";
        } else {
            s +=   "int colIdx = Ccols[bb];  ";
        }
        s += R"(
                    y[colIdx*dim + c] -= temp;
                }
            }
            )";
        return s;
    }


    std::string get_ilu_decomp_string() {
        return R"(

        // a = a - (b * c)
        __kernel void block_mult_sub(__global double *a, __local double *b, __global double *c)
        {
            const unsigned int block_size = 3;
            const unsigned int hwarp_size = 16;
            const unsigned int idx_t = get_local_id(0);                   // thread id in work group
            const unsigned int thread_id_in_hwarp = idx_t % hwarp_size;   // thread id in warp (16 threads)
            if(thread_id_in_hwarp < block_size * block_size){
                const unsigned int row = thread_id_in_hwarp / block_size;
                const unsigned int col = thread_id_in_hwarp % block_size;
                double temp = 0.0;
                for (unsigned int k = 0; k < block_size; k++) {
                    temp += b[block_size * row + k] * c[block_size * k + col];
                }
                a[block_size * row + col] -= temp;
            }
        }

        // c = a * b
        __kernel void block_mult(__global double *a, __global double *b, __local double *c) {
            const unsigned int block_size = 3;
            const unsigned int hwarp_size = 16;
            const unsigned int idx_t = get_local_id(0);                   // thread id in work group
            const unsigned int thread_id_in_hwarp = idx_t % hwarp_size;   // thread id in warp (16 threads)
            if(thread_id_in_hwarp < block_size * block_size){
                const unsigned int row = thread_id_in_hwarp / block_size;
                const unsigned int col = thread_id_in_hwarp % block_size;
                double temp = 0.0;
                for (unsigned int k = 0; k < block_size; k++) {
                    temp += a[block_size * row + k] * b[block_size * k + col];
                }
                c[block_size * row + col] = temp;
            }
        }

        // invert 3x3 matrix
        __kernel void inverter(__global double *matrix, __global double *inverse) {
            const unsigned int block_size = 3;
            const unsigned int bs = block_size;                           // rename to shorter name
            const unsigned int hwarp_size = 16;
            const unsigned int idx_t = get_local_id(0);                   // thread id in work group
            const unsigned int thread_id_in_hwarp = idx_t % hwarp_size;   // thread id in warp (16 threads)
            if(thread_id_in_hwarp < bs * bs){
                double t4  = matrix[0] * matrix[4];
                double t6  = matrix[0] * matrix[5];
                double t8  = matrix[1] * matrix[3];
                double t10 = matrix[2] * matrix[3];
                double t12 = matrix[1] * matrix[6];
                double t14 = matrix[2] * matrix[6];

                double det = (t4 * matrix[8] - t6 * matrix[7] - t8 * matrix[8] +
                              t10 * matrix[7] + t12 * matrix[5] - t14 * matrix[4]);
                double t17 = 1.0 / det;

                const unsigned int r = thread_id_in_hwarp / bs;
                const unsigned int c = thread_id_in_hwarp % bs;
                const unsigned int r1 = (r+1) % bs;
                const unsigned int c1 = (c+1) % bs;
                const unsigned int r2 = (r+bs-1) % bs;
                const unsigned int c2 = (c+bs-1) % bs;
                inverse[c*bs+r] = ((matrix[r1*bs+c1] * matrix[r2*bs+c2]) - (matrix[r1*bs+c2] * matrix[r2*bs+c1])) * t17;
            }
        }

        __kernel void ilu_decomp(const unsigned int firstRow,
                                 const unsigned int lastRow,
                                 __global double *LUvals,
                                 __global const int *LUcols,
                                 __global const int *LUrows,
                                 __global double *invDiagVals,
                                 __global int *diagIndex,
                                 const unsigned int Nb,
                                 __local double *pivot){

            const unsigned int bs = 3;
            const unsigned int hwarp_size = 16;
            const unsigned int work_group_size = get_local_size(0);
            const unsigned int work_group_id = get_group_id(0);
            const unsigned int num_groups = get_num_groups(0);
            const unsigned int hwarps_per_group = work_group_size / hwarp_size;
            const unsigned int thread_id_in_group = get_local_id(0);      // thread id in work group
            const unsigned int thread_id_in_hwarp = thread_id_in_group % hwarp_size;     // thread id in hwarp (16 threads)
            const unsigned int hwarp_id_in_group = thread_id_in_group / hwarp_size;
            const unsigned int lmem_offset = hwarp_id_in_group * bs * bs;  // each workgroup gets some lmem, but the workitems have to share it
                                                                           // every workitem in a hwarp has the same lmem_offset

            // go through all rows
            for (int i = firstRow + work_group_id * hwarps_per_group + hwarp_id_in_group; i < lastRow; i += num_groups * hwarps_per_group)
            {
                int iRowStart = LUrows[i];
                int iRowEnd = LUrows[i + 1];

                // go through all elements of the row
                for (int ij = iRowStart; ij < iRowEnd; ij++) {
                    int j = LUcols[ij];

                    if (j < i) {
                        // calculate the pivot of this row
                        block_mult(LUvals + ij * bs * bs, invDiagVals + j * bs * bs, pivot + lmem_offset);

                        // copy pivot
                        if (thread_id_in_hwarp < bs * bs) {
                            LUvals[ij * bs * bs + thread_id_in_hwarp] = pivot[lmem_offset + thread_id_in_hwarp];
                        }

                        int jRowEnd = LUrows[j + 1];
                        int jk = diagIndex[j] + 1;
                        int ik = ij + 1;
                        // substract that row scaled by the pivot from this row.
                        while (ik < iRowEnd && jk < jRowEnd) {
                            if (LUcols[ik] == LUcols[jk]) {
                                block_mult_sub(LUvals + ik * bs * bs, pivot + lmem_offset, LUvals + jk * bs * bs);
                                ik++;
                                jk++;
                            } else {
                                if (LUcols[ik] < LUcols[jk])
                                { ik++; }
                                else
                                { jk++; }
                            }
                        }
                    }
                }

                // store the inverse in the diagonal
                inverter(LUvals + diagIndex[i] * bs * bs, invDiagVals + i * bs * bs);

                // copy inverse
                if (thread_id_in_hwarp < bs * bs) {
                    LUvals[diagIndex[i] * bs * bs + thread_id_in_hwarp] = invDiagVals[i * bs * bs + thread_id_in_hwarp];
                }

            }
        }
        )";
    }

} // end namespace bda

