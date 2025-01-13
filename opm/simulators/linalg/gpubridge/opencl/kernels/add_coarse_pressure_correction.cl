/// add the coarse pressure solution back to the finer, complete solution
/// every workitem handles one blockrow
__kernel void add_coarse_pressure_correction(
    __global const double *coarse_x,
    __global double *fine_x,
    const unsigned int pressure_idx,
    const unsigned int Nb)
{
    const unsigned int NUM_THREADS = get_global_size(0);
    const unsigned int block_size = 3;
    unsigned int target_block_row = get_global_id(0);

    while(target_block_row < Nb){
        fine_x[target_block_row * block_size + pressure_idx] += coarse_x[target_block_row];
        target_block_row += NUM_THREADS;
    }
}
