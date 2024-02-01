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

#ifndef OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
#define OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED

#include <opm/simulators/linalg/MILU.hpp>

#include <opm/simulators/linalg/linalgproperties.hh>
#include <opm/models/utils/parametersystem.hh>

namespace Opm {
template <class TypeTag>
class ISTLSolverBda;
template <class TypeTag>
class ISTLSolver;
}



namespace Opm::Properties {

namespace TTag {
struct FlowIstlSolverParams {};
}

template<class TypeTag, class MyTypeTag>
struct LinearSolverReduction {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct RelaxedLinearSolverReduction {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct LinearSolverMaxIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverRestart {
    using type = UndefinedProperty;
};
//
// LinearSolverVerbosity defined in opm-models
//
template<class TypeTag, class MyTypeTag>
struct IluRelaxation {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct IluFillinLevel {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MiluVariant {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct IluRedblack {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct IluReorderSpheres {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct UseGmres {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverIgnoreConvergenceFailure{
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ScaleLinearSystem {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolver {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverPrintJsonDefinition {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct CprReuseSetup {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct CprReuseInterval {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct AcceleratorMode {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct BdaDeviceId {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct OpenclPlatformId {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct OpenclIluParallel {
    using type = UndefinedProperty;
};
template<class TypeTag>
struct LinearSolverReduction<TypeTag, TTag::FlowIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
template<class TypeTag>
struct RelaxedLinearSolverReduction<TypeTag, TTag::FlowIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
template<class TypeTag>
struct LinearSolverMaxIter<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 200;
};
template<class TypeTag>
struct LinearSolverRestart<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 40;
};
template<class TypeTag>
struct LinearSolverVerbosity<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct IluRelaxation<TypeTag, TTag::FlowIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.9;
};
template<class TypeTag>
struct IluFillinLevel<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct MiluVariant<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = "ILU";
};
template<class TypeTag>
struct IluRedblack<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct IluReorderSpheres<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct UseGmres<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct LinearSolverIgnoreConvergenceFailure<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct ScaleLinearSystem<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = "ilu0";
};
template<class TypeTag>
struct LinearSolverPrintJsonDefinition<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = true;
};
template<class TypeTag>
struct CprReuseSetup<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 4;
};
template<class TypeTag>
struct CprReuseInterval<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 30;
};
template<class TypeTag>
struct AcceleratorMode<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = "none";
};
template<class TypeTag>
struct BdaDeviceId<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct OpenclPlatformId<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct OpenclIluParallel<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = true; // note: false should only be used in debug
};

// Set the backend to be used.
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowIstlSolverParams> {
#if COMPILE_BDA_BRIDGE
    using type = ISTLSolverBda<TypeTag>;
#else
    using type = ISTLSolver<TypeTag>;
#endif
};
} // namespace Opm::Properties

namespace Opm
{


    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct FlowLinearSolverParameters
    {
        double linear_solver_reduction_;
        double relaxed_linear_solver_reduction_;
        int    linear_solver_maxiter_;
        int    linear_solver_restart_;
        int    linear_solver_verbosity_;
        double ilu_relaxation_;
        int    ilu_fillin_level_;
        MILU_VARIANT   ilu_milu_;
        bool   ilu_redblack_;
        bool   ilu_reorder_sphere_;
        bool   newton_use_gmres_;
        bool   ignoreConvergenceFailure_;
        bool scale_linear_system_;
        std::string linsolver_;
        bool linear_solver_print_json_definition_;
        int cpr_reuse_setup_;
        int cpr_reuse_interval_;
        std::string accelerator_mode_;
        int bda_device_id_;
        int opencl_platform_id_;
        bool opencl_ilu_parallel_;

        template <class TypeTag>
        void init(bool cprRequestedInDataFile)
        {
            // TODO: these parameters have undocumented non-trivial dependencies
            linear_solver_reduction_ = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
            relaxed_linear_solver_reduction_ = EWOMS_GET_PARAM(TypeTag, double, RelaxedLinearSolverReduction);
            linear_solver_maxiter_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
            linear_solver_restart_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverRestart);
            linear_solver_verbosity_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
            ilu_relaxation_ = EWOMS_GET_PARAM(TypeTag, double, IluRelaxation);
            ilu_fillin_level_ = EWOMS_GET_PARAM(TypeTag, int, IluFillinLevel);
            ilu_milu_ = convertString2Milu(EWOMS_GET_PARAM(TypeTag, std::string, MiluVariant));
            ilu_redblack_ = EWOMS_GET_PARAM(TypeTag, bool, IluRedblack);
            ilu_reorder_sphere_ = EWOMS_GET_PARAM(TypeTag, bool, IluReorderSpheres);
            newton_use_gmres_ = EWOMS_GET_PARAM(TypeTag, bool, UseGmres);
            ignoreConvergenceFailure_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure);
            scale_linear_system_ = EWOMS_GET_PARAM(TypeTag, bool, ScaleLinearSystem);
            linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolver);
            linear_solver_print_json_definition_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverPrintJsonDefinition);
            cpr_reuse_setup_  =  EWOMS_GET_PARAM(TypeTag, int, CprReuseSetup);
            cpr_reuse_interval_  =  EWOMS_GET_PARAM(TypeTag, int, CprReuseInterval);

            if (!EWOMS_PARAM_IS_SET(TypeTag, std::string, LinearSolver) && cprRequestedInDataFile) {
                linsolver_ = "cpr";
            } else {
                linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolver);
            }

            accelerator_mode_ = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
            bda_device_id_ = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
            opencl_platform_id_ = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
            opencl_ilu_parallel_ = EWOMS_GET_PARAM(TypeTag, bool, OpenclIluParallel);
        }

        template <class TypeTag>
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverReduction, "The minimum reduction of the residual which the linear solver must achieve");
            EWOMS_REGISTER_PARAM(TypeTag, double, RelaxedLinearSolverReduction, "The minimum reduction of the residual which the linear solver need to achieve for the solution to be accepted");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIter, "The maximum number of iterations of the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverRestart, "The number of iterations after which GMRES is restarted");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity, "The verbosity level of the linear solver (0: off, 2: all)");
            EWOMS_REGISTER_PARAM(TypeTag, double, IluRelaxation, "The relaxation factor of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, int, IluFillinLevel, "The fill-in level of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, MiluVariant, "Specify which variant of the modified-ILU preconditioner ought to be used. Possible variants are: ILU (default, plain ILU), MILU_1 (lump diagonal with dropped row entries), MILU_2 (lump diagonal with the sum of the absolute values of the dropped row  entries), MILU_3 (if diagonal is positive add sum of dropped row entrires. Otherwise subtract them), MILU_4 (if diagonal is positive add sum of dropped row entrires. Otherwise do nothing");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluRedblack, "Use red-black partitioning for the ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluReorderSpheres, "Whether to reorder the entries of the matrix in the red-black ILU preconditioner in spheres starting at an edge. If false the original ordering is preserved in each color. Otherwise why try to ensure D4 ordering (in a 2D structured grid, the diagonal elements are consecutive).");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseGmres, "Use GMRES as the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure, "Continue with the simulation like nothing happened after the linear solver did not converge");
            EWOMS_REGISTER_PARAM(TypeTag, bool, ScaleLinearSystem, "Scale linear system according to equation scale and primary variable types");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, LinearSolver, "Configuration of solver. Valid options are: ilu0 (default), dilu, cprw, cpr (an alias for cprw), cpr_quasiimpes, cpr_trueimpes, cpr_trueimpesanalytic, amg or hybrid (experimental). Alternatively, you can request a configuration to be read from a JSON file by giving the filename here, ending with '.json.'");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverPrintJsonDefinition, "Write the JSON definition of the linear solver setup to the DBG file.");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprReuseSetup, "Reuse preconditioner setup. Valid options are 0: recreate the preconditioner for every linear solve, 1: recreate once every timestep, 2: recreate if last linear solve took more than 10 iterations, 3: never recreate, 4: recreated every CprReuseInterval");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprReuseInterval, "Reuse preconditioner interval. Used when CprReuseSetup is set to 4, then the preconditioner will be fully recreated instead of reused every N linear solve, where N is this parameter.");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, AcceleratorMode, "Choose a linear solver, usage: '--accelerator-mode=[none|cusparse|opencl|amgcl|rocalution|rocsparse]'");
            EWOMS_REGISTER_PARAM(TypeTag, int, BdaDeviceId, "Choose device ID for cusparseSolver or openclSolver, use 'nvidia-smi' or 'clinfo' to determine valid IDs");
            EWOMS_REGISTER_PARAM(TypeTag, int, OpenclPlatformId, "Choose platform ID for openclSolver, use 'clinfo' to determine valid platform IDs");
            EWOMS_REGISTER_PARAM(TypeTag, bool, OpenclIluParallel, "Parallelize ILU decomposition and application on GPU");
        }

        FlowLinearSolverParameters() { reset(); }

        // set default values
        void reset()
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
            linsolver_                = "ilu0";
            linear_solver_print_json_definition_ = true;
            cpr_reuse_setup_          = 4;
            cpr_reuse_interval_       = 30;
            accelerator_mode_         = "none";
            bda_device_id_            = 0;
            opencl_platform_id_       = 0;
            opencl_ilu_parallel_      = true;
        }
    };


} // namespace Opm




#endif // OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
