
# Case Study Simulation Configurations
We use the following command to run flow with 16 processes in parallel for the CPU benchmark, where we specify that the linear solver is described in a json file called "cpu_ilu0.json".
```bash
    mpirun -n 16 flow \
        /path/to/DATA/file \
        --output-dir=/path/to/output \
        --linear-solver=/path/to/cpu_ilu0.json \
        --threads-per-process=1 \
        --newton-min-iterations=1 \
        --matrix-add-well-contributions=true
```
To recreate the CPU runs use the following "cpu_ilu0.json".
```json
{
    "tol": "0.01",
    "maxiter": "200",
    "verbosity": "0",
    "solver": "bicgstab",
    "preconditioner": {
        "type": "ILU0",
        "ilulevel": "0"
    }
}
```

We use the following command to run the GPU simulations, where we only use one process. To save simulation time we use 16 threads in parallel, this does not effect the timing for the linear solver which exists entirely on the GPU.
```bash
    mpirun -n 1 flow \
        /path/to/DATA/file \
        --output-dir=/path/to/output \
        --linear-solver=/path/to/gpu_ilu0.json \
        --threads-per-process=16 \
        --newton-min-iterations=1 \
        --matrix-add-well-contributions=true
```
The "gpu_ilu0.json" file contains arguments specifying that we want to use the split matrix storage format and automatically tune the blocksizes of the GPU kernels.
```json
{
    "tol": "0.01",
    "maxiter": "200",
    "verbosity": "0",
    "solver": "gpubicgstab",
    "preconditioner": {
        "type": "OPMGPUILU0",
        "ilulevel": "0"
    }
}
```
