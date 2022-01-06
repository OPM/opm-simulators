__kernel void block_sub(__global double *mat1, __global double *mat2, __global double *result)
{
    const unsigned int bs = 3;
    const unsigned int warpsize = 32;
    const unsigned int num_active_threads = (warpsize / bs / bs) * bs * bs;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int lane = idx_t % warpsize;

    if(lane < num_active_threads){
        const unsigned int row = lane % bs;
        const unsigned int col = (lane / bs) % bs;
        result[bs * row + col] = mat1[bs * row + col] - mat2[bs * row + col];
    }
}

__kernel void block_mult_sub_isai(__global double *a, __global double *b, __global double *c)
{
    const unsigned int bs = 3;
    const unsigned int warpsize = 32;
    const unsigned int num_active_threads = (warpsize / bs / bs) * bs * bs;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int lane = idx_t % warpsize;

    if(lane < num_active_threads){
        const unsigned int row = lane % bs;
        const unsigned int col = (lane / bs) % bs;
        double temp = 0.0;

        for (unsigned int k = 0; k < bs; k++) {
            temp += b[bs * row + k] * c[bs * k + col];
        }

        a[bs * row + col] -= temp;
    }
}

__kernel void isaiL(__global const int *diagIndex,
                    __global const int *colPtr,
                    __global const int *mapping,
                    __global const int *nvc,
                    __global const int *luIdxs,
                    __global const int *xxIdxs,
                    __global const int *dxIdxs,
                    __global const double *LU,
                    __global double *invL,
                    const unsigned int Nb)
{
    const unsigned int warpsize = 32;
    const unsigned int idx_b = get_group_id(0);
    const unsigned int idx_t = get_local_id(0);
    const unsigned int idx = get_global_id(0);
    const unsigned int bs = 3;
    const unsigned int num_threads = get_global_size(0);
    const unsigned int num_warps_in_grid = num_threads / warpsize;
    const unsigned int num_active_threads = (warpsize / bs / bs) * bs * bs;
    const unsigned int num_blocks_per_warp = warpsize / bs / bs;
    const unsigned int lane = idx_t % warpsize;
    const unsigned int c = (lane / bs) % bs;
    const unsigned int r = lane % bs;
    unsigned int tcol = idx / warpsize;

    while(tcol < Nb){
        const unsigned int frow = diagIndex[tcol] + 1;
        const unsigned int lrow = colPtr[tcol + 1];
        const unsigned int nx = lrow - frow;

        if(lane < num_active_threads){
            for(unsigned int xid = 0; xid < nx; xid++){
                unsigned int xpos = mapping[frow + xid];
                block_sub(invL + xpos * bs * bs, LU + xpos * bs * bs, invL + xpos * bs * bs);
            }

            for(unsigned int v = nvc[tcol]; v < nvc[tcol + 1]; v++){
                block_mult_sub_isai(invL + xxIdxs[v] * bs * bs, LU + luIdxs[v] * bs * bs, invL + dxIdxs[v] * bs * bs);
            }
        }

        tcol += num_warps_in_grid;
    }
}
