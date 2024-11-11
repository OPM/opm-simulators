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

#include <opm/simulators/linalg/linalgparameters.hh>
#include <opm/simulators/linalg/linalgproperties.hh>

namespace Opm {

template <class TypeTag>
class ISTLSolverGpuBridge;

template <class TypeTag>
class ISTLSolver;

}

namespace Opm::Properties {

namespace TTag {

struct FlowIstlSolverParams {};

}

// Set the backend to be used.
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowIstlSolverParams>
{
#if COMPILE_GPU_BRIDGE
    using type = ISTLSolverGpuBridge<TypeTag>;
#else
    using type = ISTLSolver<TypeTag>;
#endif
};

}

namespace Opm::Parameters {

struct LinearSolverReduction { static constexpr double value = 1e-2; };
struct NlddLocalLinearSolverReduction { static constexpr double value = 1e-2; };
struct RelaxedLinearSolverReduction { static constexpr double value = 1e-2; };
struct IluRelaxation { static constexpr double value = 0.9; };
struct LinearSolverMaxIter { static constexpr int value = 200; };
struct NlddLocalLinearSolverMaxIter { static constexpr int value = 200; };
struct LinearSolverRestart { static constexpr int value = 40; };
struct IluFillinLevel { static constexpr int value = 0; };
struct MiluVariant { static constexpr auto value = "ILU"; };
struct IluRedblack { static constexpr bool value = false; };
struct IluReorderSpheres { static constexpr bool value = false; };
struct UseGmres { static constexpr bool value = false; };
struct LinearSolverIgnoreConvergenceFailure { static constexpr bool value = false; };
struct ScaleLinearSystem { static constexpr bool value = false; };
struct LinearSolver { static constexpr auto value = "cprw"; };
struct NlddLocalLinearSolver { static constexpr auto value = "ilu0"; };
struct LinearSolverPrintJsonDefinition { static constexpr auto value = true; };
struct CprReuseSetup { static constexpr int value = 4; };
struct CprReuseInterval { static constexpr int value = 30; };
struct AcceleratorMode { static constexpr auto value = "none"; };
struct GpuDeviceId { static constexpr int value = 0; };
struct OpenclPlatformId { static constexpr int value = 0; };
struct OpenclIluParallel { static constexpr bool value = true; }; // note: false should only be used in debug

} // namespace Opm::Parameters

namespace Opm {

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
    bool is_nldd_local_solver_;
    std::string linsolver_;
    bool linear_solver_print_json_definition_;
    int cpr_reuse_setup_;
    int cpr_reuse_interval_;
    std::string accelerator_mode_;
    int gpu_device_id_;
    int opencl_platform_id_;
    bool opencl_ilu_parallel_;

    FlowLinearSolverParameters() { reset(); }

    void init(bool cprRequestedInDataFile);

    static void registerParameters();

    // set default values
    void reset();
};

} // namespace Opm

#endif // OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
