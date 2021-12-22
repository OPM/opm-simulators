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
