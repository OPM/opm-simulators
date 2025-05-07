---
title: gpu-ISTL - Extending OPM Flow with GPU Linear Solvers
tags:
  - Reservoir Simulation
  - OPM
  - GPU Computing
  - High Performance Computing
  - Linear Solvers
  - CUDA
  - HIP
authors:
  - name: Kjetil Olsen Lye
    orcid: 0000-0001-6914-1010
    affiliation: 1
  - name: Tobias Meyer Andersen
    orcid: 0009-0000-8913-3034
    affiliation: 1
  - given-names: Atgeirr Flø Rasmussen
    orcid: 0000-0002-7932-3835
    affiliation: 1
  - name: Jakob Torben
    orcid: 0009-0002-9036-3739
    affiliation: 1
affiliations:
 - name: Mathematics and Cybernetics, SINTEF Digital, Oslo, Norway
   index: 1
date: 18 October 2024
bibliography: paper.bib
---

# Summary

The `gpu-ISTL` framework provides `Dune-ISTL` compatible sparse linear operations on the GPU for the OPM Flow simulator [@OPMFLOW]. OPM Flow is an open-source fully implicit simulator for subsurface reservoir simulation of industrial complexity, where efficient linear solvers are critical for computational performance. It is written in C++ and relies on the Distributed and Unified Numerics Environment (Dune) [@BASTIAN202175] for several numerical algorithms. In particular, templated linear-solver algorithms provided by the Dune Iterative Solver Template Library (`Dune-ISTL`) [@dune-istl] are essential in the simulator. To GPU-accelerate the simulator with minimal effort, `gpu-ISTL` provides classes that can instantiate `Dune-ISTL` algorithms with types that automatically runs the algorithms provided on the GPU. Both AMD and Nvidia GPUs are supported in the same source code by utilizing hipify-perl [@hipify] to generate HIP code [@HIP] from the existing CUDA code [@cuda] when compiling OPM Flow.

# Statement of Need

Simulating multiphase and multicomponent fluid behavior within complex geological formations is crucial for modern geoenergy operations. The repeated solution of large, sparse systems of linear equations is typically the most computationally expensive part of these simulations. Therefore, it is essential for a reservoir simulator like OPM Flow, which is used both commercially and for research, to have flexible linear algebra libraries that optimally utilize contemporary computer hardware. To fill this need, we introduce `gpu-ISTL`, a generic `Dune-ISTL` compatible GPU-accelerated linear solver library and employ it within OPM Flow.

A current and prominent trend in computer architecture is hardware specialization [@computer_arch_hennessy_patterson], where modern computational units are increasingly tailored to specific workloads. Having to write tailored implementations of each numerical method for every possible hardware provider is generally undesirable, especially in a large codebase like OPM Flow, which exceeds half a million lines of code. The `gpu-ISTL` library leverages the parallel computational resources of GPUs to meet industry performance demands. It also facilitates further research into the development of efficient GPU-accelerated linear solvers, without requiring intrusive and substantial changes to the OPM Flow codebase.

Our new library enhances research in numerical methods for reservoir simulation by leveraging OPM Flow's extensive infrastructure for complex reservoirs. It enables investigations into GPU-based linear solvers, preconditioners, autotuning and mixed-precision computations, see for instance [@andersen_ilu_dilu_rsc]. By integrating `gpu-ISTL` into the OPM Flow linear solver subsystem, development remains synchronized with the simulator's rapidly evolving codebase, ensuring compatibility and continuity in development efforts.

Several standalone libraries GPU-accelerate linear algebra operations. Some also support at least both Nvidia and AMD cards, as well as the iterative solvers required by OPM Flow. Examples of open-source libraries fulfilling all these criteria are `ViennaCL` [@ViennaCL] and `PETSc` [@mills2021], though those are not `Dune-ISTL` compatible. The Bandicoot library [@curtin_bandicoot_2023], an extension to the Armadillo [@armadillo] library, is another worthy mention of a GPU library that effectively computes linear algebra operations, although it does not feature the iterative solvers that are vital to OPM.

# The `gpu-ISTL` Components

The `gpu-ISTL` library encompasses numerical algorithms, Dune adapters, and essential linear algebra components such as `GpuVector` and `GpuSparseMatrix` that leverage optimized libraries from Nvidia and AMD for efficient mathematical operations. Preconditioners such as Jacobi, ILU(0), and DILU are implemented with custom GPU kernels to enhance the performance of the bi-conjugate gradient stabilized (BiCGSTAB) method within `gpu-ISTL`. Adapter classes like `SolverAdapter` and `PreconditionerAdapter` allow mixing GPU and CPU solvers or preconditioners for ease of development, testing and validation. Moreover, `GpuOwnerOverlapCopy` extends the MPI functionality in Dune to support multi-GPU simulations, including CUDA-aware MPI for Nvidia cards to accelerate inter-process memory transfers. The library also provides capabilities for autotuning thread-block sizes in user-implemented GPU kernels. Figure 1 contains a simplified class diagram with of some of the components that constitute `gpu-ISTL`.

![Class diagram showing a simplified view of how `gpu-ISTL` is implemented. The colored backgrounds indicate which namespace the classes belong to. Colors on individual classes correspond to what conceptual part of the simulator they are a part of. White represents linear solvers, red represents the classes enabling multiprocess simulations, orange represents the implementation linear solver preconditioners, and green represents general linear algebra functionality.](figures/mpiparallel_gpuistl.png)

Furthermore, Figure 2 provides a sequence diagram for the call stack of a single linear solve in OPM Flow using the gpu-ISTL framework. We note that the sequence diagram is for the MPI free version. In Figure 2 we denote the (indirect) calls to specific GPU kernels as an arrow from the corresponding `gpu-ISTL` classes to the right-hand red box representing the GPU.

![A sequence diagram for a single linear solver in OPM Flow using the `gpu-ISTL` framework. For readability, components from OPM are marked in green, components from DUNE in yellow and components from `gpu-ISTL` are marked in red. The classes GpuVector, GpuSparseMatrix and GpuILU0 of `gpu-ISTL` will (through some indirection) call GPU kernels, which is visualized as a big red block on the right-hand side of the diagram.](figures/timeloopgpuistlproper.png)

OPM Flow's build system utilizes hipify-perl [@hipify], which translates CUDA [@cuda] code to HIP [@HIP] code if one wants to compile for AMD GPUs. Incorporating the translation in the build system ensures that the hardware of the two largest GPU vendors are supported without incurring any extra overhead for the developers of the CUDA code in the repository.

# Performance Case Study of $CO_2$ Storage Simulation

To demonstrate the effectiveness of `gpu-ISTL`, we conduct a scaling study on geological carbon storage cases derived from the 11th Society of Petroleum Engineers (SPE) Comparative Solutions Project (CSP) [@Nordbotten2024]. Specifically, we simulate Case C from the SPE11 CSP using Pyopmspe11 [@pyopmspe11] with successively finer discretizations.

We simulate using the ILU0 preconditioned BiCGSTAB as the linear solver on three different hardware configurations, tracking the time spent in the linear solver normalized by the number of linear iterations. Specifically, we utilize 16 cores of an AMD 7950X CPU for the first run, an AMD Radeon RX 7900XTX GPU for the second run, and an Nvidia RTX 4070Ti GPU for the third run. Figure 3 shows the speedup of the GPU compared to the normalized CPU performance. The results demonstrate that GPU runtimes scale better than CPU runtimes as problem sizes increase. This highlights the performance benefits of using the advanced linear solvers with preconditioners implemented in `gpu-ISTL`. Furthermore, it underscores how `gpu-ISTL` and OPM Flow can serve as a robust platform for exploring and testing novel preconditioners within an industry-relevant environment, offering rapid evaluation on both synthetic and real-life cases.

![Speedup of the mean time per linear iteration across a 2000-year simulation case derived from Case C of the 11th SPE Comparative Solution Project compared to the CPU implementation. A speedup of 5.6 and 4.8 is achieved on the largest case for the RTX 4070Ti and Radeon RX 7900XTX respectively.](figures/lin_its_small.svg)

# Case Study Configuration
The cases generated by Pyopmspe11 can be found in [@spe11c_dataset], in the ZIP-file "SPE11C_LONG.zip". The file names include the sidelength of the domain in terms of cells to distinguish between the cases. The total number of cells in the domain is then the cube of that number, although Pyopmspe11 will add some extra cells as padding on the domain.

To run OPM Flow simulations for the case study on CPU, we use the following command to run OPM Flow with 16 processes in parallel, where we specify that the linear solver is described in a JSON file called "cpu_ilu0.json". The "\<S\>" represents the sidelength in the file path and name.
```bash
    mpirun -n 16 flow \
        /path/to/SPE11C_LONG/SPE11C<S>CUBE_LONG/deck/SPE11C<S>CUBE_LONG.DATA \
        --output-dir=/path/to/output \
        --linear-solver=/path/to/cpu_ilu0.json \
        --threads-per-process=1 \
        --newton-min-iterations=1 \
        --matrix-add-well-contributions=true
```
To recreate the CPU runs, use the following "cpu_ilu0.json".
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

We use the command below to run the GPU simulations, where we only use one process, because we run it on a single GPU. To save simulation time we use 16 threads in parallel, which does not affect the timing for the linear solver which exists entirely on the GPU. The "\<S\>" represents the sidelength in the file path and name.
```bash
    mpirun -n 1 flow \
        /path/to/SPE11C_LONG/SPE11C<S>CUBE_LONG/deck/SPE11C<S>CUBE_LONG.DATA \
        --output-dir=/path/to/output \
        --linear-solver=/path/to/gpu_ilu0.json \
        --threads-per-process=16 \
        --newton-min-iterations=1 \
        --matrix-add-well-contributions=true
```
The "gpu_ilu0.json" file enables the use of `gpu-ISTL` by specifying that we want to use the `gpubicgstab` linear solver, and the `OPMGPUILU0` preconditioner.
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
# Acknowledgements

Development of `gpu-ISTL` in OPM Flow from 2021 to 2024 was part of the EU project *HPC, Big Data, and Artificial Intelligence convergent platform* (ACROSS). Development in 2024 has been financed by Equinor.

# References

