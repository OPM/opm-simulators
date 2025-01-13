/// Applies sdtwells
__kernel void stdwell_apply(
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

        // merge all blocks in this workgroup into 1 block
        // if numBlocksPerWarp >= 3, should use loop
        // block 1:     block 2:
        //  0  1  2     12 13 14
        //  3  4  5     15 16 17
        //  6  7  8     18 19 20
        //  9 10 11     21 22 23
        // workitem i will hold the sum of workitems i and i + valsPerBlock
        if(wiId < valsPerBlock){
            for (int i = 1; i < numBlocksPerWarp; ++i) {
                localSum[wiId] += localSum[wiId + i*valsPerBlock];
	    }
        }

        if(c == 0 && wiId < valsPerBlock){
            for(unsigned int i = dim - 1; i > 0; --i){
                localSum[wiId] += localSum[wiId + i];
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
