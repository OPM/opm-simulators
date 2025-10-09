/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/linalgparameters.hh>
#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void FlowLinearSolverParameters::init(bool cprRequestedInDataFile)
{
    // TODO: these parameters have undocumented non-trivial dependencies
    relaxed_linear_solver_reduction_ = Parameters::Get<Parameters::RelaxedLinearSolverReduction>();
    linear_solver_restart_ = Parameters::Get<Parameters::LinearSolverRestart>();
    linear_solver_verbosity_ = Parameters::Get<Parameters::LinearSolverVerbosity>();
    ilu_relaxation_ = Parameters::Get<Parameters::IluRelaxation>();
    ilu_fillin_level_ = Parameters::Get<Parameters::IluFillinLevel>();
    ilu_milu_ = convertString2Milu(Parameters::Get<Parameters::MiluVariant>());
    ilu_redblack_ = Parameters::Get<Parameters::IluRedblack>();
    ilu_reorder_sphere_ = Parameters::Get<Parameters::IluReorderSpheres>();
    newton_use_gmres_ = Parameters::Get<Parameters::UseGmres>();
    ignoreConvergenceFailure_ = Parameters::Get<Parameters::LinearSolverIgnoreConvergenceFailure>();
    scale_linear_system_ = Parameters::Get<Parameters::ScaleLinearSystem>();
    linsolver_ = Parameters::Get<Parameters::LinearSolver>();
    linear_solver_print_json_definition_ = Parameters::Get<Parameters::LinearSolverPrintJsonDefinition>();
    cpr_reuse_setup_  = Parameters::Get<Parameters::CprReuseSetup>();
    cpr_reuse_interval_  = Parameters::Get<Parameters::CprReuseInterval>();
    gpu_aware_mpi_ = Parameters::Get<Parameters::GpuAwareMpi>();

    if (!Parameters::IsSet<Parameters::LinearSolver>() && cprRequestedInDataFile) {
        linsolver_ = "cpr";
    } else {
        linsolver_ = Parameters::Get<Parameters::LinearSolver>();
    }

    // if local solver from nldd, use nldd local linear solver parameters
    if (is_nldd_local_solver_) {
        linsolver_ = Parameters::Get<Parameters::NlddLocalLinearSolver>();
        linear_solver_reduction_ = Parameters::Get<Parameters::NlddLocalLinearSolverReduction>();
        linear_solver_maxiter_ = Parameters::Get<Parameters::NlddLocalLinearSolverMaxIter>();
    }
    else {
        linsolver_ = Parameters::Get<Parameters::LinearSolver>();
        linear_solver_reduction_ = Parameters::Get<Parameters::LinearSolverReduction>();
        linear_solver_maxiter_ = Parameters::Get<Parameters::LinearSolverMaxIter>();
    }

    accelerator_mode_ = Parameters::Get<Parameters::AcceleratorMode>();
    cpr_weights_thread_parallel_ = Parameters::Get<Parameters::CprWeightsThreadParallel>();
    gpu_device_id_ = Parameters::Get<Parameters::GpuDeviceId>();
    opencl_platform_id_ = Parameters::Get<Parameters::OpenclPlatformId>();
    opencl_ilu_parallel_ = Parameters::Get<Parameters::OpenclIluParallel>();

    // NLDD local solvers must use CPU accelerator because extractMatrix and other
    // operations in the NLDD implementation are not GPU-compatible yet.
    if (is_nldd_local_solver_) {
        linear_solver_accelerator_ = Parameters::LinearSolverAcceleratorType::CPU;
    } else {
        linear_solver_accelerator_ = Parameters::linearSolverAcceleratorTypeFromCLI();
    }

    if (linear_solver_accelerator_ == Parameters::LinearSolverAcceleratorType::GPU) {
        if (!Parameters::IsSet<Parameters::LinearSolver>()) {
            // Since well operator is not implemented for GPU we can not use default cprw
            linsolver_ = "cpr_trueimpes";
        }
    }
}

void FlowLinearSolverParameters::registerParameters()
{
    Parameters::Register<Parameters::LinearSolverReduction>
        ("The minimum reduction of the residual which the linear solver must achieve");
    Parameters::Register<Parameters::NlddLocalLinearSolverReduction>
        ("The minimum reduction of the residual which the NLDD local linear solver must achieve");
    Parameters::Register<Parameters::RelaxedLinearSolverReduction>
        ("The minimum reduction of the residual which the linear solver need to "
         "achieve for the solution to be accepted");
    Parameters::Register<Parameters::LinearSolverMaxIter>
        ("The maximum number of iterations of the linear solver");
    Parameters::Register<Parameters::NlddLocalLinearSolverMaxIter>
        ("The maximum number of iterations of the NLDD local linear solver");
    Parameters::Register<Parameters::LinearSolverRestart>
        ("The number of iterations after which GMRES is restarted");
    Parameters::Register<Parameters::LinearSolverVerbosity>
        ("The verbosity level of the linear solver (0: off, 2: all)");
    Parameters::Register<Parameters::IluRelaxation>
        ("The relaxation factor of the linear solver's ILU preconditioner");
    Parameters::Register<Parameters::IluFillinLevel>
        ("The fill-in level of the linear solver's ILU preconditioner");
    Parameters::Register<Parameters::MiluVariant>
        ("Specify which variant of the modified-ILU preconditioner ought to be used. "
         "Possible variants are: ilu (default, plain ILU), "
         "milu_1 (lump diagonal with dropped row entries), "
         "milu_2 (lump diagonal with the sum of the absolute values of the dropped row entries), "
         "milu_3 (if diagonal is positive add sum of dropped row entries, otherwise subtract them), "
         "milu_4 (if diagonal is positive add sum of dropped row entries, otherwise do nothing)");
    Parameters::Register<Parameters::IluRedblack>
        ("Use red-black partitioning for the ILU preconditioner");
    Parameters::Register<Parameters::IluReorderSpheres>
        ("Whether to reorder the entries of the matrix in the red-black "
         "ILU preconditioner in spheres starting at an edge. "
         "If false the original ordering is preserved in each color. "
         "Otherwise why try to ensure D4 ordering (in a 2D structured grid, "
         "the diagonal elements are consecutive).");
    Parameters::Register<Parameters::UseGmres>
        ("Use GMRES as the linear solver");
    Parameters::Register<Parameters::LinearSolverIgnoreConvergenceFailure>
        ("Continue with the simulation like nothing happened "
         "after the linear solver did not converge");
    Parameters::Register<Parameters::ScaleLinearSystem>
        ("Scale linear system according to equation scale and primary variable types");
    Parameters::Register<Parameters::LinearSolver>
        ("Configuration of solver. Valid options are: cprw (default), "
         "ilu0, dilu, cpr (an alias for cprw), cpr_quasiimpes, "
         "cpr_trueimpes, cpr_trueimpesanalytic, amg or hybrid (experimental). "
         "Alternatively, you can request a configuration to be read from a "
         "JSON file by giving the filename here, ending with '.json.'");
    Parameters::Register<Parameters::NlddLocalLinearSolver>
        ("Configuration of NLDD local linear solver. Valid options are: ilu0 (default), "
            "dilu, cpr_quasiimpes and amg. "
            "Alternatively, you can request a configuration to be read from a "
            "JSON file by giving the filename here, ending with '.json.'");
    Parameters::Register<Parameters::LinearSolverPrintJsonDefinition>
        ("Write the JSON definition of the linear solver setup to the DBG file.");
    Parameters::Register<Parameters::CprReuseSetup>
        ("Reuse preconditioner setup. Valid options are "
         "0: recreate the preconditioner for every linear solve, "
         "1: recreate once every timestep, "
         "2: recreate if last linear solve took more than 10 iterations, "
         "3: never recreate, "
         "4: recreated every CprReuseInterval");
    Parameters::Register<Parameters::CprReuseInterval>
        ("Reuse preconditioner interval. Used when CprReuseSetup is set to 4, "
         "then the preconditioner will be fully recreated instead of reused "
         "every N linear solve, where N is this parameter.");
    Parameters::Register<Parameters::AcceleratorMode>
        ("Choose a linear solver, usage: "
         "'--accelerator-mode=[none|cusparse|opencl|amgcl|rocalution|rocsparse]'");
    Parameters::Register<Parameters::GpuDeviceId>
        ("Choose device ID for cusparseSolver or openclSolver, "
         "use 'nvidia-smi' or 'clinfo' to determine valid IDs");
    Parameters::Register<Parameters::OpenclPlatformId>
        ("Choose platform ID for openclSolver, use 'clinfo' "
         "to determine valid platform IDs");
    Parameters::Register<Parameters::OpenclIluParallel>
        ("Parallelize ILU decomposition and application on GPU");
    Parameters::Register<Parameters::LinearSolverAccelerator>
        ("Choose the backend for the linear solver, usage: "
         "'--linear-solver-accelerator=[cpu|gpu]'.");
    Parameters::Register<Parameters::GpuAwareMpi>
        ("MPI communication use GPU aware MPI in the sense that "
            "it will use GPU direct communication. Setting this to true "
            " will require that the MPI implementation "
            "supports GPU direct communication. "
            "If you are unsure, set this to false. "
            "Usage: --gpu-aware-mpi=[true|false]. ");
    Parameters::Register<Parameters::VerifyGpuAwareMpi>
        ("Verify that the MPI implementation supports GPU aware MPI. "
            "If this is set to true *and* --gpu-aware-mpi=true, the simulation will fail if the "
            "MPI implementation does not support GPU aware MPI. "
            "Note that the verification is not exhaustive, "
            "and some configurations will not verify, but will work in practice. "
            "Usage: --verify-gpu-aware-mpi=[true|false]. ");
    Parameters::Register<Parameters::CprWeightsThreadParallel>
        ("Enable OpenMP thread parallelization of CPR weight calculation. "
            "This can improve performance for large models but is disabled by default."
            "Usage: --cpr-weights-thread-parallel=[true|false]. ");

    Parameters::SetDefault<Parameters::LinearSolverVerbosity>(0);
}

void FlowLinearSolverParameters::reset()
{
    relaxed_linear_solver_reduction_ = 1e-2;
    linear_solver_reduction_  = 1e-2;
    linear_solver_maxiter_    = 200;
    linear_solver_restart_    = 40;
    linear_solver_verbosity_  = 0;
    ilu_relaxation_           = 0.9;
    ilu_fillin_level_         = 0;
    ilu_milu_                 = MILU_VARIANT::ILU;
    ilu_redblack_             = false;
    ilu_reorder_sphere_       = false;
    newton_use_gmres_         = false;
    ignoreConvergenceFailure_ = false;
    scale_linear_system_      = false;
    is_nldd_local_solver_     = false;
    linsolver_                = "cprw";
    linear_solver_print_json_definition_ = true;
    cpr_reuse_setup_          = 4;
    cpr_reuse_interval_       = 30;
    accelerator_mode_         = "none";
    gpu_device_id_            = 0;
    opencl_platform_id_       = 0;
    opencl_ilu_parallel_      = true;
    linear_solver_accelerator_ = Parameters::LinearSolverAcceleratorType::CPU;
    gpu_aware_mpi_              = false;
    verify_gpu_aware_mpi_       = false;
    cpr_weights_thread_parallel_ = false;
}

} // namespace Opm
