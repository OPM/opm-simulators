/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
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

#include <ewoms/common/parametersystem.hh>

#include <array>
#include <memory>

namespace Opm {
template <class TypeTag>
class ISTLSolverEbos;
}


BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowIstlSolverParams);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(LinearSolverReduction);
NEW_PROP_TAG(IluRelaxation);
NEW_PROP_TAG(LinearSolverMaxIter);
NEW_PROP_TAG(LinearSolverRestart);
NEW_PROP_TAG(FlowLinearSolverVerbosity);
NEW_PROP_TAG(IluFillinLevel);
NEW_PROP_TAG(MiluVariant);
NEW_PROP_TAG(IluRedblack);
NEW_PROP_TAG(IluReorderSpheres);
NEW_PROP_TAG(UseGmres);
NEW_PROP_TAG(LinearSolverRequireFullSparsityPattern);
NEW_PROP_TAG(LinearSolverIgnoreConvergenceFailure);
NEW_PROP_TAG(UseAmg);
NEW_PROP_TAG(UseCpr);
NEW_PROP_TAG(LinearSolverBackend);
NEW_PROP_TAG(PreconditionerAddWellContributions);
NEW_PROP_TAG(SystemStrategy);
NEW_PROP_TAG(ScaleLinearSystem);
NEW_PROP_TAG(CprSolverVerbose);
NEW_PROP_TAG(CprUseDrs);
NEW_PROP_TAG(CprMaxEllIter);
NEW_PROP_TAG(CprEllSolvetype);
NEW_PROP_TAG(CprReuseSetup);
NEW_PROP_TAG(LinearSolverConfigurationJsonFile);

SET_SCALAR_PROP(FlowIstlSolverParams, LinearSolverReduction, 1e-2);
SET_SCALAR_PROP(FlowIstlSolverParams, IluRelaxation, 0.9);
SET_INT_PROP(FlowIstlSolverParams, LinearSolverMaxIter, 200);
SET_INT_PROP(FlowIstlSolverParams, LinearSolverRestart, 40);
SET_INT_PROP(FlowIstlSolverParams, FlowLinearSolverVerbosity, 0);
SET_INT_PROP(FlowIstlSolverParams, IluFillinLevel, 0);
SET_STRING_PROP(FlowIstlSolverParams, MiluVariant, "ILU");
SET_BOOL_PROP(FlowIstlSolverParams, IluRedblack, false);
SET_BOOL_PROP(FlowIstlSolverParams, IluReorderSpheres, false);
SET_BOOL_PROP(FlowIstlSolverParams, UseGmres, false);
SET_BOOL_PROP(FlowIstlSolverParams, LinearSolverRequireFullSparsityPattern, false);
SET_BOOL_PROP(FlowIstlSolverParams, LinearSolverIgnoreConvergenceFailure, false);
SET_BOOL_PROP(FlowIstlSolverParams, UseAmg, false);
SET_BOOL_PROP(FlowIstlSolverParams, UseCpr, false);
SET_TYPE_PROP(FlowIstlSolverParams, LinearSolverBackend, Opm::ISTLSolverEbos<TypeTag>);
SET_BOOL_PROP(FlowIstlSolverParams, PreconditionerAddWellContributions, false);
SET_STRING_PROP(FlowIstlSolverParams, SystemStrategy, "none");
SET_BOOL_PROP(FlowIstlSolverParams, ScaleLinearSystem, false);
SET_INT_PROP(FlowIstlSolverParams, CprSolverVerbose, 0);
SET_BOOL_PROP(FlowIstlSolverParams, CprUseDrs, false);
SET_INT_PROP(FlowIstlSolverParams, CprMaxEllIter, 20);
SET_INT_PROP(FlowIstlSolverParams, CprEllSolvetype, 0);
SET_INT_PROP(FlowIstlSolverParams, CprReuseSetup, 0);
SET_STRING_PROP(FlowIstlSolverParams, LinearSolverConfigurationJsonFile, "none");



END_PROPERTIES

namespace Opm
{

    /**
     * \brief Parameters used to configure the CPRPreconditioner.
     */
    struct CPRParameter
    {
        double cpr_relax_;
        double cpr_solver_tol_;
        int cpr_ilu_n_;
        MILU_VARIANT cpr_ilu_milu_;
        bool cpr_ilu_redblack_;
        bool cpr_ilu_reorder_sphere_;
        bool cpr_use_drs_;
        int cpr_max_ell_iter_;
        int cpr_ell_solvetype_;
        bool cpr_use_amg_;
        bool cpr_use_bicgstab_;
        int cpr_solver_verbose_;
        bool cpr_pressure_aggregation_;
        int cpr_reuse_setup_;
        CPRParameter() { reset(); }

        void reset()
        {
            cpr_solver_tol_           = 1e-2;
            cpr_ilu_n_                = 0;
            cpr_ilu_milu_             = MILU_VARIANT::ILU;
            cpr_ilu_redblack_         = false;
            cpr_ilu_reorder_sphere_   = true;
            cpr_max_ell_iter_         = 25;
            cpr_ell_solvetype_        = 0;
            cpr_use_drs_              = false;
            cpr_use_amg_              = true;
            cpr_use_bicgstab_         = true;
            cpr_solver_verbose_       = 0;
            cpr_pressure_aggregation_ = false;
            cpr_reuse_setup_          = 0;
        }
    };




    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct FlowLinearSolverParameters
        : public CPRParameter
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
        bool   linear_solver_use_amg_;
        bool   use_cpr_;
        std::string system_strategy_;
        bool scale_linear_system_;
        std::string linear_solver_configuration_json_file_;

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
            linear_solver_use_amg_ = EWOMS_GET_PARAM(TypeTag, bool, UseAmg);
            use_cpr_ = EWOMS_GET_PARAM(TypeTag, bool, UseCpr);
            system_strategy_ = EWOMS_GET_PARAM(TypeTag, std::string, SystemStrategy);
            scale_linear_system_ = EWOMS_GET_PARAM(TypeTag, bool, ScaleLinearSystem);
            cpr_solver_verbose_  =  EWOMS_GET_PARAM(TypeTag, int, CprSolverVerbose);
            cpr_use_drs_  =  EWOMS_GET_PARAM(TypeTag, bool, CprUseDrs);
            cpr_max_ell_iter_  =  EWOMS_GET_PARAM(TypeTag, int, CprMaxEllIter);
            cpr_ell_solvetype_  =  EWOMS_GET_PARAM(TypeTag, int, CprEllSolvetype);
            cpr_reuse_setup_  =  EWOMS_GET_PARAM(TypeTag, int, CprReuseSetup);
            linear_solver_configuration_json_file_ = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolverConfigurationJsonFile);
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
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseAmg, "Use AMG as the linear solver's preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseCpr, "Use CPR as the linear solver's preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, SystemStrategy, "Strategy for reformulating and scaling linear system (none: no scaling -- should not be used with CPR, original: use weights that are equivalent to no scaling -- should not be used with CPR, simple: form pressure equation as simple sum of conservation equations, quasiimpes: form pressure equation based on diagonal block, trueimpes: form pressure equation based on linearization of accumulation term)");
            EWOMS_REGISTER_PARAM(TypeTag, bool, ScaleLinearSystem, "Scale linear system according to equation scale and primary variable types");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprSolverVerbose, "Verbosity of cpr solver (0: silent, 1: print summary of inner linear solver, 2: print extensive information about inner linear solve, including setup information)");
            EWOMS_REGISTER_PARAM(TypeTag, bool, CprUseDrs, "Use dynamic row sum using weights");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprMaxEllIter, "MaxIterations of the elliptic pressure part of the cpr solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprEllSolvetype, "Solver type of elliptic pressure solve (0: bicgstab, 1: cg, 2: only amg preconditioner)");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprReuseSetup, "Reuse Amg Setup");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, LinearSolverConfigurationJsonFile, "Filename of JSON configuration for flexible linear solver system.");
        }

        FlowLinearSolverParameters() { reset(); }

        // set default values
        void reset()
        {
            use_cpr_     = false;
            newton_use_gmres_        = false;
            linear_solver_reduction_ = 1e-2;
            linear_solver_maxiter_   = 150;
            linear_solver_restart_   = 40;
            linear_solver_verbosity_ = 0;
            require_full_sparsity_pattern_ = false;
            ignoreConvergenceFailure_ = false;
            linear_solver_use_amg_    = false;
            ilu_fillin_level_         = 0;
            ilu_relaxation_           = 0.9;
            ilu_milu_                 = MILU_VARIANT::ILU;
            ilu_redblack_             = false;
            ilu_reorder_sphere_       = true;
        }
    };


} // namespace Opm




#endif // OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
