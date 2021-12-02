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
