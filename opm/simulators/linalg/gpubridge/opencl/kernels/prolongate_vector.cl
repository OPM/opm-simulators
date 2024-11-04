/// prolongate vector during amg cycle
/// every workitem handles one row
__kernel void prolongate_vector(
    __global const double *in,
    __global double *out,
    __global const int *cols,
    const unsigned int N)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    unsigned int row = get_global_id(0);

    while(row < N){
        out[row] += in[cols[row]];
        row += NUM_THREADS;
    }
}
