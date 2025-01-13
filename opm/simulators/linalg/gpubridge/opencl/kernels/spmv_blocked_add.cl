/// b = mat * x
/// algorithm based on:
/// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
/// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
__kernel void spmv_blocked_add(
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
            out[row] += x[row];
        }
        target_block_row += num_warps_in_grid;
    }
}
