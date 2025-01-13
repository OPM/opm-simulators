/// transform blocked vector to scalar vector using pressure-weights
/// every workitem handles one blockrow
__kernel void full_to_pressure_restriction(
    __global const double *fine_y,
    __global const double *weights,
    __global double *coarse_y,
    const unsigned int Nb)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    const unsigned int block_size = 3;
    unsigned int target_block_row = get_global_id(0);

    while(target_block_row < Nb){
        double sum = 0.0;
        unsigned int idx = block_size * target_block_row;
        for (unsigned int i = 0; i < block_size; ++i) {
            sum += fine_y[idx + i] * weights[idx + i];
        }
        coarse_y[target_block_row] = sum;
        target_block_row += NUM_THREADS;
    }
}
