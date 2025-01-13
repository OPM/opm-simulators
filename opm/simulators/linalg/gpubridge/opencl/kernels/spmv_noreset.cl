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
