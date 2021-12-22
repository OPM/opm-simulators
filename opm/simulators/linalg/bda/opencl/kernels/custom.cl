/// Custom kernel: combines some bicgstab vector operations into 1
/// p = (p - omega * v) * beta + r
__kernel void custom(
    __global double *p,
    __global double *v,
    __global double *r,
    const double omega,
    const double beta,
    const int N)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    unsigned int idx = get_global_id(0);

    while(idx < N){
        double res = p[idx];
        res -= omega * v[idx];
        res *= beta;
        res += r[idx];
        p[idx] = res;
        idx += NUM_THREADS;
    }
}
