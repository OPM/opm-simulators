/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS
  Copyright 2015 IRIS AS

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

#ifndef OPM_NEWTONITERATIONBLACKOILINTERLEAVED_HEADER_INCLUDED
#define OPM_NEWTONITERATIONBLACKOILINTERLEAVED_HEADER_INCLUDED

#include <opm/autodiff/CPRPreconditioner.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/autodiff/ParallelOverlappingILU0.hpp>

#include <ewoms/common/parametersystem.hh>

#include <array>
#include <memory>

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

END_PROPERTIES

namespace Opm
{
    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class
    struct NewtonIterationBlackoilInterleavedParameters
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
        }

        NewtonIterationBlackoilInterleavedParameters() { reset(); }
        // read values from parameter class
        NewtonIterationBlackoilInterleavedParameters( const ParameterGroup& param )
            : CPRParameter(param)
        {
            // set default parameters
            reset();

            // read parameters (using previsouly set default values)
            newton_use_gmres_        = param.getDefault("newton_use_gmres", newton_use_gmres_ );
            linear_solver_reduction_ = param.getDefault("linear_solver_reduction", linear_solver_reduction_ );
            linear_solver_maxiter_   = param.getDefault("linear_solver_maxiter", linear_solver_maxiter_);
            linear_solver_restart_   = param.getDefault("linear_solver_restart", linear_solver_restart_);
            linear_solver_verbosity_ = param.getDefault("linear_solver_verbosity", linear_solver_verbosity_);
            require_full_sparsity_pattern_ = param.getDefault("require_full_sparsity_pattern", require_full_sparsity_pattern_);
            ignoreConvergenceFailure_ = param.getDefault("linear_solver_ignoreconvergencefailure", ignoreConvergenceFailure_);
            linear_solver_use_amg_    = param.getDefault("linear_solver_use_amg", linear_solver_use_amg_ );
            ilu_relaxation_           = param.getDefault("ilu_relaxation", ilu_relaxation_ );
            ilu_fillin_level_         = param.getDefault("ilu_fillin_level",  ilu_fillin_level_ );
            ilu_redblack_             = param.getDefault("ilu_redblack", cpr_ilu_redblack_);
            ilu_reorder_sphere_       = param.getDefault("ilu_reorder_sphere", cpr_ilu_reorder_sphere_);
            std::string milu("ILU");
            ilu_milu_ = convertString2Milu(param.getDefault("ilu_milu", milu));

            // Check whether to use cpr approach
            const std::string cprSolver = "cpr";
            use_cpr_ = ( param.getDefault("solver_approach", std::string()) == cprSolver );
        }

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


    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables).
    class NewtonIterationBlackoilInterleaved : public NewtonIterationBlackoilInterface
    {
    public:

        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of the linear solvers
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        NewtonIterationBlackoilInterleaved(const ParameterGroup& param,
                                           const boost::any& parallelInformation=boost::any());

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const;

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        virtual int iterations () const { return iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        virtual const boost::any& parallelInformation() const;

    private:
        // max number of equations supported, increase if necessary
        static const int maxNumberEquations_ = 6 ;

        mutable std::array< std::unique_ptr< NewtonIterationBlackoilInterface >, maxNumberEquations_+1 > newtonIncrementDoublePrecision_;
        mutable std::array< std::unique_ptr< NewtonIterationBlackoilInterface >, maxNumberEquations_+1 > newtonIncrementSinglePrecision_;
        NewtonIterationBlackoilInterleavedParameters parameters_;
        boost::any parallelInformation_;
        mutable int iterations_;
    };

} // namespace Opm


#endif // OPM_NEWTONITERATIONBLACKOILINTERLEAVED_HEADER_INCLUDED
