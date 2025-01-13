/// multiply vector with another vector and a scalar, element-wise
/// add result to a third vector
__kernel void vmul(
    const double alpha,
    __global double const *in1,
    __global double const *in2,
    __global double *out,
    const int N)
{
    unsigned int NUM_THREADS = get_global_size(0);
    int idx = get_global_id(0);

    while(idx < N){
        out[idx] += alpha * in1[idx] * in2[idx];
        idx += NUM_THREADS;
    }
}
