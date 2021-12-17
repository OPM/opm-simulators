// This file is auto-generated. Do not edit!

#include "openclKernels.hpp"

namespace Opm{

namespace Accelerator{

const std::string OpenclKernels::ILU_apply1_fm_str = R"( 
/// ILU apply part 1: forward substitution.
/// Solves L*x=y where L is a lower triangular sparse blocked matrix.
/// Here, L is inside a normal, square matrix.
/// In this case, diagIndex indicates where the rows of L end.
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
        const unsigned int last_block = diagIndex[target_block_row];
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

const std::string OpenclKernels::ILU_apply2_fm_str = R"( 
/// ILU apply part 2: backward substitution.
/// Solves U*x=y where U is an upper triangular sparse blocked matrix.
/// Here, U is inside a normal, square matrix.
/// In this case diagIndex indicates where the rows of U start.
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
    const unsigned int lane = idx_t % warpsize;
    const unsigned int c = (lane / bs) % bs;
    const unsigned int r = lane % bs;

    target_block_row += nodesPerColorPrefix[color];

    while(target_block_row < nodesPerColorPrefix[color+1]){
        const unsigned int first_block = diagIndex[target_block_row] + 1;
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

const std::string OpenclKernels::ILU_decomp_str = R"( 
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
__kernel void block_mult(__global double *a, __global double *b, __local double *c)
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
            temp += a[block_size * row + k] * b[block_size * k + col];
        }
        c[block_size * row + col] = temp;
    }
}

// invert 3x3 matrix
__kernel void inverter(__global double *matrix, __global double *inverse)
{
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

/// Exact ilu decomposition kernel
/// The kernel takes a full BSR matrix and performs inplace ILU decomposition
__kernel void ilu_decomp(const unsigned int firstRow,
                         const unsigned int lastRow,
                         __global double *LUvals,
                         __global const int *LUcols,
                         __global const int *LUrows,
                         __global double *invDiagVals,
                         __global int *diagIndex,
                         const unsigned int Nb,
                         __local double *pivot)
{
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
                // subtract that row scaled by the pivot from this row.
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

const std::string OpenclKernels::add_coarse_pressure_correction_str = R"( 
/// add the coarse pressure solution back to the finer, complete solution
/// every workitem handles one blockrow
__kernel void add_coarse_pressure_correction(
    __global const double *coarse_x,
    __global double *fine_x,
    const unsigned int pressure_idx,
    const unsigned int Nb)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    const unsigned int block_size = 3;
    unsigned int target_block_row = get_global_id(0);

    while(target_block_row < Nb){
        fine_x[target_block_row * block_size + pressure_idx] += coarse_x[target_block_row];
        target_block_row += NUM_THREADS;
    }
}
)"; 

const std::string OpenclKernels::axpy_str = R"( 
/// axpy kernel: a = a + alpha * b
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

const std::string OpenclKernels::custom_str = R"( 
/// Custom kernel: combines some bicgstab vector operations into 1
/// p = (p - omega * v) * beta + r
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

const std::string OpenclKernels::dot_1_str = R"( 
/// returns partial sums, instead of the final dot product
/// partial sums are added on CPU
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

const std::string OpenclKernels::full_to_pressure_restriction_str = R"( 
/// transform blocked vector to scalar vector using pressure-weights
/// every workitem handles one blockrow
__kernel void full_to_pressure_restriction(
    __global const double *fine_y,
    __global const double *weights,
    __global double *coarse_y,
    const unsigned int Nb)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    const unsigned int block_size = 3;
    unsigned int target_block_row = get_global_id(0);

    while(target_block_row < Nb){
        double sum = 0.0;
        unsigned int idx = block_size * target_block_row;
        for (unsigned int i = 0; i < block_size; ++i) {
            sum += fine_y[idx + i] * weights[idx + i];
        }
        coarse_y[target_block_row] = sum;
        target_block_row += NUM_THREADS;
    }
}
)"; 

const std::string OpenclKernels::norm_str = R"( 
/// returns partial sums, instead of the final norm
/// the square root must be computed on CPU
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

const std::string OpenclKernels::prolongate_vector_str = R"( 
/// prolongate vector during amg cycle
/// every workitem handles one row
__kernel void prolongate_vector(
    __global const double *in,
    __global double *out,
    __global const int *cols,
    const unsigned int N)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    unsigned int row = get_global_id(0);

    while(row < N){
        out[row] += in[cols[row]];
        row += NUM_THREADS;
    }
}
)"; 

const std::string OpenclKernels::residual_str = R"( 
/// res = rhs - mat * x
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void residual(
    __global const double *vals,
    __global const int *cols,
    __global const int *rows,
    const int N,
    __global const double *x,
    __global const double *rhs,
    __global double *out,
    __local double *tmp)
{
    const unsigned int bsize = get_local_size(0);
    const unsigned int idx_b = get_global_id(0) / bsize;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int num_workgroups = get_num_groups(0);

    int row = idx_b;

    while (row < N) {
        int rowStart = rows[row];
        int rowEnd = rows[row+1];
        int rowLength = rowEnd - rowStart;
        double local_sum = 0.0;
        for (int j = rowStart + idx_t; j < rowEnd; j += bsize) {
            int col = cols[j];
            local_sum += vals[j] * x[col];
        }

        tmp[idx_t] = local_sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        int offset = bsize / 2;
        while(offset > 0) {
            if (idx_t < offset) {
                tmp[idx_t] += tmp[idx_t + offset];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            offset = offset / 2;
        }

        if (idx_t == 0) {
            out[row] = rhs[row] - tmp[idx_t];
        }

        row += num_workgroups;
    }
}
)"; 

const std::string OpenclKernels::residual_blocked_str = R"( 
/// res = rhs - mat * x
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void residual_blocked(
    __global const double *vals,
    __global const int *cols,
    __global const int *rows,
    const int Nb,
    __global const double *x,
    __global const double *rhs,
    __global double *out,
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
            out[row] = rhs[row] - tmp[lane];
        }
        target_block_row += num_warps_in_grid;
    }
}
)"; 

const std::string OpenclKernels::scale_str = R"( 
/// scale vector with scalar: a = a * alpha
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

const std::string OpenclKernels::spmv_str = R"( 
/// b = mat * x
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void spmv(
    __global const double *vals,
    __global const int *cols,
    __global const int *rows,
    const int N,
    __global const double *x,
    __global double *out,
    __local double *tmp)
{
    const unsigned int bsize = get_local_size(0);
    const unsigned int idx_b = get_global_id(0) / bsize;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int num_workgroups = get_num_groups(0);

    int row = idx_b;

    while (row < N) {
        int rowStart = rows[row];
        int rowEnd = rows[row+1];
        int rowLength = rowEnd - rowStart;
        double local_sum = 0.0;
        for (int j = rowStart + idx_t; j < rowEnd; j += bsize) {
            int col = cols[j];
            local_sum += vals[j] * x[col];
        }

        tmp[idx_t] = local_sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        int offset = bsize / 2;
        while(offset > 0) {
            if (idx_t < offset) {
                tmp[idx_t] += tmp[idx_t + offset];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            offset = offset / 2;
        }

        if (idx_t == 0) {
            out[row] = tmp[idx_t];
        }

        row += num_workgroups;
    }
}
)"; 

const std::string OpenclKernels::spmv_blocked_str = R"( 
/// b = mat * x
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void spmv_blocked(
    __global const double *vals,
    __global const int *cols,
    __global const int *rows,
    const int Nb,
    __global const double *x,
    __global double *out,
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
            out[row] = tmp[lane];
        }
        target_block_row += num_warps_in_grid;
    }
}
)"; 

const std::string OpenclKernels::spmv_noreset_str = R"( 
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void spmv_noreset(
    __global const double *vals,
    __global const int *cols,
    __global const int *rows,
    const int N,
    __global const double *x,
    __global double *out,
    __local double *tmp)
{
    const unsigned int bsize = get_local_size(0);
    const unsigned int idx_b = get_global_id(0) / bsize;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int num_workgroups = get_num_groups(0);

    int row = idx_b;

    while (row < N) {
        int rowStart = rows[row];
        int rowEnd = rows[row+1];
        int rowLength = rowEnd - rowStart;
        double local_sum = 0.0;
        for (int j = rowStart + idx_t; j < rowEnd; j += bsize) {
            int col = cols[j];
            local_sum += vals[j] * x[col];
        }

        tmp[idx_t] = local_sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        int offset = bsize / 2;
        while(offset > 0) {
            if (idx_t < offset) {
                tmp[idx_t] += tmp[idx_t + offset];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            offset = offset / 2;
        }

        if (idx_t == 0) {
            out[row] += tmp[idx_t];
        }

        row += num_workgroups;
    }
}
)"; 

const std::string OpenclKernels::stdwell_apply_str = R"( 
/// In this kernel there is reordering: the B/Ccols do not correspond with the x/y vector
/// the x/y vector is reordered, using toOrder to address that
__kernel void stdwell_apply(
            __global const double *Cnnzs,
            __global const double *Dnnzs,
            __global const double *Bnnzs,
            __global const int *Ccols,
            __global const int *Bcols,
            __global const double *x,
            __global double *y,
            __global const int *toOrder,
            const unsigned int dim,
            const unsigned int dim_wells,
            __global const unsigned int *val_pointers,
            __local double *localSum,
            __local double *z1,
            __local double *z2)
{
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
            int colIdx = toOrder[Bcols[b]];
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

        int colIdx = toOrder[Ccols[bb]];
        y[colIdx*dim + c] -= temp;
    }
}
)"; 

const std::string OpenclKernels::stdwell_apply_no_reorder_str = R"( 
/// Applies sdtwells without reordering
__kernel void stdwell_apply_no_reorder(
            __global const double *Cnnzs,
            __global const double *Dnnzs,
            __global const double *Bnnzs,
            __global const int *Ccols,
            __global const int *Bcols,
            __global const double *x,
            __global double *y,
            const unsigned int dim,
            const unsigned int dim_wells,
            __global const unsigned int *val_pointers,
            __local double *localSum,
            __local double *z1,
            __local double *z2)
{
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
            int colIdx = Bcols[b];
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

        int colIdx = Ccols[bb];
        y[colIdx*dim + c] -= temp;
    }
}
)"; 

const std::string OpenclKernels::vmul_str = R"( 
/// multiply vector with another vector and a scalar, element-wise
/// add result to a third vector
__kernel void vmul(
    const double alpha,
    __global double const *in1,
    __global double const *in2,
    __global double *out,
    const int N)
{
    unsigned int NUM_THREADS = get_global_size(0);
    int idx = get_global_id(0);

    while(idx < N){
        out[idx] += alpha * in1[idx] * in2[idx];
        idx += NUM_THREADS;
    }
}
)"; 

}
}