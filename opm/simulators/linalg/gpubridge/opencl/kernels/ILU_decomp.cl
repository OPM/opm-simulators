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
                         __global const unsigned *rowIndices,
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
        const unsigned row = rowIndices[i];
        int iRowStart = LUrows[row];
        int iRowEnd = LUrows[row + 1];

        // go through all elements of the row
        for (int ij = iRowStart; ij < iRowEnd; ij++) {
            int j = LUcols[ij];

            if (j < row) {
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
        inverter(LUvals + diagIndex[row] * bs * bs, invDiagVals + row * bs * bs);

        // copy inverse
        if (thread_id_in_hwarp < bs * bs) {
            LUvals[diagIndex[row] * bs * bs + thread_id_in_hwarp] = invDiagVals[row * bs * bs + thread_id_in_hwarp];
        }
    }
}
