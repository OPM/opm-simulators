/// scale vector with scalar: a = a * alpha
__kernel void scale(
    __global double *vec,
    const double a,
    const int N)
{
    unsigned int NUM_THREADS = get_global_size(0);
    int idx = get_global_id(0);

    while(idx < N){
        vec[idx] *= a;
        idx += NUM_THREADS;
    }
}
