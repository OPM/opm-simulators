__kernel void block_add(__global double *mat1, __global double *mat2, __global double *result)
{
    const unsigned int bs = 3;
    const unsigned int warpsize = 32;
    const unsigned int num_active_threads = (warpsize / bs / bs) * bs * bs;
    const unsigned int idx_t = get_local_id(0);
    const unsigned int lane = idx_t % warpsize;

    if(lane < num_active_threads){
        const unsigned int row = lane % bs;
        const unsigned int col = (lane / bs) % bs;
        result[bs * row + col] = mat1[bs * row + col] + mat2[bs * row + col];
    }
}

__kernel void block_mult_isai(__global double *mat1, __global double *mat2, __global double *result)
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
            temp += mat1[bs * row + k] * mat2[bs * k + col];
        }

        result[bs * row + col] = temp;
    }
}

__kernel void isaiU(__global const int *diagIndex,
                    __global const int *colPtr,
                    __global const int *rowIndices,
                    __global const int *mapping,
                    __global const int *nvc,
                    __global const int *luIdxs,
                    __global const int *xxIdxs,
                    __global const int *dxIdxs,
                    __global const double *LU,
                    __global const double *invDiagVals,
                    __global double *invU,
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
        const unsigned int frow = colPtr[tcol];
        const unsigned int lrow = diagIndex[tcol];
        const unsigned int nx = lrow - frow + 1;

        if(lane < num_active_threads){
            block_add(invU + lrow * bs * bs, invDiagVals + tcol * bs * bs, invU + lrow * bs * bs);

            for(unsigned int v = nvc[tcol]; v < nvc[tcol + 1]; v++){
                block_mult_sub_isai(invU + xxIdxs[v] * bs * bs, LU + luIdxs[v] * bs * bs, invU + dxIdxs[v] * bs * bs);
            }

            for(unsigned int xid = 1; xid < nx; xid++){
                unsigned int xpos = mapping[lrow - xid];
                block_mult_isai(invDiagVals + rowIndices[lrow - xid] * bs * bs, invU + xpos * bs * bs, invU + xpos * bs * bs);
            }
        }

        tcol += num_warps_in_grid;
    }
}
