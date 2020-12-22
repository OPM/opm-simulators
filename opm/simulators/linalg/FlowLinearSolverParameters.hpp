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

#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>

#include <opm/simulators/linalg/linalgproperties.hh>
#include <opm/models/utils/parametersystem.hh>

#include <array>
#include <memory>

namespace Opm {
template <class TypeTag>
class ISTLSolverEbos;
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
struct IluRelaxation {
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
template<class TypeTag, class MyTypeTag>
struct FlowLinearSolverVerbosity {
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
struct LinearSolverRequireFullSparsityPattern {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverIgnoreConvergenceFailure{
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct PreconditionerAddWellContributions {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ScaleLinearSystem {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct CprMaxEllIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct CprEllSolvetype {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct CprReuseSetup {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct Linsolver {
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
struct OpenclIluReorder {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct FpgaBitstream {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct LinearSolverReduction<TypeTag, TTag::FlowIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
template<class TypeTag>
struct IluRelaxation<TypeTag, TTag::FlowIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.9;
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
struct FlowLinearSolverVerbosity<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
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
struct LinearSolverRequireFullSparsityPattern<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct LinearSolverIgnoreConvergenceFailure<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowIstlSolverParams> {
    using type = Opm::ISTLSolverEbos<TypeTag>;
};
template<class TypeTag>
struct PreconditionerAddWellContributions<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct ScaleLinearSystem<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct CprMaxEllIter<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 20;
};
template<class TypeTag>
struct CprEllSolvetype<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct CprReuseSetup<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr int value = 3;
};
template<class TypeTag>
struct Linsolver<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = "ilu0";
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
struct OpenclIluReorder<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = ""; // note: default value is chosen depending on the solver used
};
template<class TypeTag>
struct FpgaBitstream<TypeTag, TTag::FlowIstlSolverParams> {
    static constexpr auto value = "";
};

} // namespace Opm::Properties

namespace Opm
{


    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct FlowLinearSolverParameters
    {
        double linear_solver_reduction_;
        double ilu_relaxation_;
        int    linear_solver_maxiter_;
        int    linear_solver_restart_;
        int    linear_solver_verbosity_;
        int    ilu_fillin_level_;
        Opm::MILU_VARIANT   ilu_milu_;
        bool   ilu_redblack_;
        bool   ilu_reorder_sphere_;
        bool   newton_use_gmres_;
        bool   require_full_sparsity_pattern_;
        bool   ignoreConvergenceFailure_;
        bool scale_linear_system_;
        std::string linsolver_;
        std::string accelerator_mode_;
        int bda_device_id_;
        int opencl_platform_id_;
        int cpr_max_ell_iter_ = 20;
        int cpr_reuse_setup_ = 0;
        std::string opencl_ilu_reorder_;
        std::string fpga_bitstream_;

        template <class TypeTag>
        void init()
        {
            // TODO: these parameters have undocumented non-trivial dependencies
            linear_solver_reduction_ = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
            ilu_relaxation_ = EWOMS_GET_PARAM(TypeTag, double, IluRelaxation);
            linear_solver_maxiter_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
            linear_solver_restart_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverRestart);
            linear_solver_verbosity_ = EWOMS_GET_PARAM(TypeTag, int, FlowLinearSolverVerbosity);
            ilu_fillin_level_ = EWOMS_GET_PARAM(TypeTag, int, IluFillinLevel);
            ilu_milu_ = convertString2Milu(EWOMS_GET_PARAM(TypeTag, std::string, MiluVariant));
            ilu_redblack_ = EWOMS_GET_PARAM(TypeTag, bool, IluRedblack);
            ilu_reorder_sphere_ = EWOMS_GET_PARAM(TypeTag, bool, IluReorderSpheres);
            newton_use_gmres_ = EWOMS_GET_PARAM(TypeTag, bool, UseGmres);
            require_full_sparsity_pattern_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverRequireFullSparsityPattern);
            ignoreConvergenceFailure_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure);
            scale_linear_system_ = EWOMS_GET_PARAM(TypeTag, bool, ScaleLinearSystem);
            cpr_max_ell_iter_  =  EWOMS_GET_PARAM(TypeTag, int, CprMaxEllIter);
            cpr_reuse_setup_  =  EWOMS_GET_PARAM(TypeTag, int, CprReuseSetup);
            linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, Linsolver);
            accelerator_mode_ = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
            bda_device_id_ = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
            opencl_platform_id_ = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
            opencl_ilu_reorder_ = EWOMS_GET_PARAM(TypeTag, std::string, OpenclIluReorder);
            fpga_bitstream_ = EWOMS_GET_PARAM(TypeTag, std::string, FpgaBitstream);
        }

        template <class TypeTag>
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverReduction, "The minimum reduction of the residual which the linear solver must achieve");
            EWOMS_REGISTER_PARAM(TypeTag, double, IluRelaxation, "The relaxation factor of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIter, "The maximum number of iterations of the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverRestart, "The number of iterations after which GMRES is restarted");
            EWOMS_REGISTER_PARAM(TypeTag, int, FlowLinearSolverVerbosity, "The verbosity level of the linear solver (0: off, 2: all)");
            EWOMS_REGISTER_PARAM(TypeTag, int, IluFillinLevel, "The fill-in level of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, MiluVariant, "Specify which variant of the modified-ILU preconditioner ought to be used. Possible variants are: ILU (default, plain ILU), MILU_1 (lump diagonal with dropped row entries), MILU_2 (lump diagonal with the sum of the absolute values of the dropped row  entries), MILU_3 (if diagonal is positive add sum of dropped row entrires. Otherwise substract them), MILU_4 (if diagonal is positive add sum of dropped row entrires. Otherwise do nothing");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluRedblack, "Use red-black partioning for the ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluReorderSpheres, "Whether to reorder the entries of the matrix in the red-black ILU preconditioner in spheres starting at an edge. If false the original ordering is preserved in each color. Otherwise why try to ensure D4 ordering (in a 2D structured grid, the diagonal elements are consecutive).");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseGmres, "Use GMRES as the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverRequireFullSparsityPattern, "Produce the full sparsity pattern for the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure, "Continue with the simulation like nothing happened after the linear solver did not converge");
            EWOMS_REGISTER_PARAM(TypeTag, bool, ScaleLinearSystem, "Scale linear system according to equation scale and primary variable types");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprMaxEllIter, "MaxIterations of the elliptic pressure part of the cpr solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprReuseSetup, "Reuse preconditioner setup. Valid options are 0: recreate the preconditioner for every linear solve, 1: recreate once every timestep, 2: recreate if last linear solve took more than 10 iterations, 3: never recreate");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, Linsolver, "Configuration of solver. Valid options are: ilu0 (default), cpr (an alias for cpr_trueimpes), cpr_quasiimpes, cpr_trueimpes or amg. Alternatively, you can request a configuration to be read from a JSON file by giving the filename here, ending with '.json.'");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, AcceleratorMode, "Use GPU (cusparseSolver or openclSolver) or FPGA (fpgaSolver) as the linear solver, usage: '--accelerator-mode=[none|cusparse|opencl|fpga]'");
            EWOMS_REGISTER_PARAM(TypeTag, int, BdaDeviceId, "Choose device ID for cusparseSolver or openclSolver, use 'nvidia-smi' or 'clinfo' to determine valid IDs");
            EWOMS_REGISTER_PARAM(TypeTag, int, OpenclPlatformId, "Choose platform ID for openclSolver, use 'clinfo' to determine valid platform IDs");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, OpenclIluReorder, "Choose the reordering strategy for ILU for openclSolver and fpgaSolver, usage: '--opencl-ilu-reorder=[level_scheduling|graph_coloring], level_scheduling behaves like Dune and cusparse, graph_coloring is more aggressive and likely to be faster, but is random-based and generally increases the number of linear solves and linear iterations significantly.");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, FpgaBitstream, "Specify the bitstream file for fpgaSolver (including path), usage: '--fpga-bitstream=<filename>'");
        }

        FlowLinearSolverParameters() { reset(); }

        // set default values
        void reset()
        {
            newton_use_gmres_        = false;
            linear_solver_reduction_ = 1e-2;
            linear_solver_maxiter_   = 150;
            linear_solver_restart_   = 40;
            linear_solver_verbosity_ = 0;
            require_full_sparsity_pattern_ = false;
            ignoreConvergenceFailure_ = false;
            ilu_fillin_level_         = 0;
            ilu_relaxation_           = 0.9;
            ilu_milu_                 = MILU_VARIANT::ILU;
            ilu_redblack_             = false;
            ilu_reorder_sphere_       = true;
            accelerator_mode_         = "none";
            bda_device_id_            = 0;
            opencl_platform_id_       = 0;
            opencl_ilu_reorder_       = "";  // note: the default value is chosen depending on the solver used
            fpga_bitstream_           = "";
        }
    };


} // namespace Opm




#endif // OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
