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
