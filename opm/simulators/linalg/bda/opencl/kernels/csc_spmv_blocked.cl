#pragma OPENCL EXTENSION cl_khr_fp64: enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable

void AtomicAdd(__global double *val, double delta){
    union{
        double f;
        ulong i;
    } old;
    union {
        double f;
        ulong i;
    } new;

    do{
        old.f = *val;
        new.f = old.f + delta;
    } while(atom_cmpxchg((volatile __global ulong *)val, old.i, new.i) != old.i);
}

__kernel void csc_spmv_blocked(__global const double *vals,
                               __global const int *cptrs,
                               __global const int *rinds,
                               __global const double *x,
                               __global double *out,
                               __local double *tmp,
                               const int Nb,
                               const int block_size)
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
    unsigned int target_block_col = idx / warpsize;
    const unsigned int lane = idx_t % warpsize;
    const unsigned int c = (lane / bs) % bs;
    const unsigned int r = lane % bs;

    while(target_block_col < Nb){
        unsigned int first_block = cptrs[target_block_col];
        unsigned int last_block = cptrs[target_block_col + 1];
        unsigned int block = first_block + lane / (bs*bs);
        double local_out = 0.0;

        for(; block < last_block; block += num_blocks_per_warp){
            if(lane < num_active_threads){
                double x_elem = x[target_block_col * bs + c];
                double A_elem = vals[block * bs * bs + r * bs + c];
                local_out = A_elem * x_elem;
            }

            tmp[lane] = local_out;
            barrier(CLK_LOCAL_MEM_FENCE);

            for(unsigned int offset = bs; offset <= bs * (bs - 1); offset += bs){
                if(c == 0){
                    tmp[lane] += tmp[lane + offset];
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }

            if(lane < num_active_threads && c == 0){
                AtomicAdd(out + rinds[block] * bs + r, tmp[lane]);
                /* out[rinds[block] * bs + r] += tmp[lane]; */
            }
        }

        target_block_col += num_warps_in_grid;
    }
}
