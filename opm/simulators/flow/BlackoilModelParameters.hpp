/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
#define OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED

#include <opm/simulators/flow/SubDomain.hpp>

#include <string>

namespace Opm::Parameters {

template<class Scalar>
struct DbhpMaxRel { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct DwellFractionMax { static constexpr Scalar value = 0.2; };

struct EclDeckFileName { static constexpr auto value = ""; };

template<class Scalar>
struct InjMultOscThreshold { static constexpr Scalar value = 0.1; };

template<class Scalar>
struct InjMultDampMult { static constexpr Scalar value = 0.9; };

template<class Scalar>
struct InjMultMinDampFactor { static constexpr Scalar value = 0.05; };

template<class Scalar>
struct MaxResidualAllowed { static constexpr Scalar value = 1e7; };

template<class Scalar>
struct RelaxedMaxPvFraction { static constexpr Scalar value = 0.03; };

template<class Scalar>
struct ToleranceMb { static constexpr Scalar value = 1e-7; };

template<class Scalar>
struct ToleranceMbRelaxed { static constexpr Scalar value = 1e-6; };

template<class Scalar>
struct ToleranceEnergyBalance { static constexpr Scalar value = 1e-7; };

template<class Scalar>
struct ToleranceEnergyBalanceRelaxed { static constexpr Scalar value = 1e-6; };

template<class Scalar>
struct ToleranceCnv { static constexpr Scalar value = 1e-2; };

template<class Scalar>
struct ToleranceCnvRelaxed { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct ToleranceCnvEnergy { static constexpr Scalar value = 1e-2; };

template<class Scalar>
struct ToleranceCnvEnergyRelaxed { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct ToleranceWells { static constexpr Scalar value = 1e-4; };

template<class Scalar>
struct ToleranceWellControl { static constexpr Scalar value = 1e-7; };

struct MaxWelleqIter { static constexpr int value = 30; };

template<class Scalar>
struct MaxSinglePrecisionDays { static constexpr Scalar value = 20.0; };

struct MinStrictCnvIter { static constexpr int value = -1; };
struct MinStrictMbIter { static constexpr int value = -1; };
struct SolveWelleqInitially { static constexpr bool value = true; };
struct PreSolveNetwork { static constexpr bool value = true; };
struct UpdateEquationsScaling { static constexpr bool value = false; };
struct UseUpdateStabilization { static constexpr bool value = true; };
struct MatrixAddWellContributions { static constexpr bool value = false; };

struct UseMultisegmentWell { static constexpr bool value = true; };

template<class Scalar>
struct TolerancePressureMsWells { static constexpr Scalar value = 0.01*1e5; };

template<class Scalar>
struct MaxPressureChangeMsWells { static constexpr Scalar value = 10*1e5; };

struct MaxNewtonIterationsWithInnerWellIterations { static constexpr int value = 8; };
struct MaxInnerIterMsWells { static constexpr int value = 100; };
struct MaxInnerIterWells { static constexpr int value = 50; };
struct MaxWellStatusSwitchInInnerIterWells { static constexpr int value = 99; };
struct ShutUnsolvableWells { static constexpr bool value = true; };
struct AlternativeWellRateInit { static constexpr bool value = true; };
struct StrictOuterIterWells { static constexpr int value = 6; };
struct StrictInnerIterWells { static constexpr int value = 40; };

template<class Scalar>
struct RegularizationFactorWells { static constexpr Scalar value = 100.0; };

struct EnableWellOperabilityCheck { static constexpr bool value = true; };
struct EnableWellOperabilityCheckIter { static constexpr bool value = false; };
struct DebugEmitCellPartition { static constexpr bool value = false; };

template<class Scalar>
struct RelaxedWellFlowTol { static constexpr Scalar value = 1e-3; };

template<class Scalar>
struct RelaxedPressureTolMsw { static constexpr Scalar value = 1e4; };

struct MaximumNumberOfWellSwitches { static constexpr int value = 3; };
struct MaximumNumberOfGroupSwitches { static constexpr int value = 3; };
struct UseAverageDensityMsWells { static constexpr bool value = false; };
struct LocalWellSolveControlSwitching { static constexpr bool value = true; };
struct UseImplicitIpr { static constexpr bool value = true; };
struct CheckGroupConstraintsInnerWellIterations { static constexpr bool value = true; };

// Network solver parameters
struct NetworkMaxStrictOuterIterations { static constexpr int value = 10; };
struct NetworkMaxOuterIterations { static constexpr int value = 10; };
struct NetworkMaxSubIterations { static constexpr int value = 20; };
template<class Scalar>
struct NetworkPressureUpdateDampingFactor { static constexpr Scalar value = 0.1; };
template<class Scalar>
struct NetworkMaxPressureUpdateInBars { static constexpr Scalar value = 5.0; };
struct NonlinearSolver { static constexpr auto value = "newton"; };
struct LocalSolveApproach { static constexpr auto value = "gauss-seidel"; };
struct MaxLocalSolveIterations { static constexpr int value = 20; };
struct NewtonMinIterations { static constexpr int value = 2; };

struct WellGroupConstraintsMaxIterations { static constexpr int value = 1; };
template<class Scalar>
struct LocalToleranceScalingMb { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct LocalToleranceScalingCnv { static constexpr Scalar value = 0.1; };
struct NlddNumInitialNewtonIter { static constexpr int value = 1; };
template<class Scalar>
struct NlddRelativeMobilityChangeTol { static constexpr Scalar value = 0.1; };
struct NumLocalDomains { static constexpr int value = 0; };

template<class Scalar>
struct LocalDomainsPartitioningImbalance { static constexpr Scalar value = 1.03; };

struct LocalDomainsPartitioningMethod { static constexpr auto value = "zoltan"; };
struct LocalDomainsPartitionWellNeighborLevels { static constexpr int value = 1; };
struct LocalDomainsOrderingMeasure { static constexpr auto value = "maxpressure"; };

struct ConvergenceMonitoring { static constexpr bool value = false; };
struct ConvergenceMonitoringCutOff { static constexpr int value = 6; };
template<class Scalar>
struct ConvergenceMonitoringDecayFactor { static constexpr Scalar value = 0.75; };


template<class Scalar>
struct NupcolGroupRateTolerance { static constexpr Scalar value = 0.001; };

} // namespace Opm::Parameters

namespace Opm {

/// Solver parameters for the BlackoilModel.
template <class Scalar>
struct BlackoilModelParameters
{
public:
    /// Max relative change in bhp in single iteration.
    Scalar dbhp_max_rel_;
    /// Max absolute change in well volume fraction in single iteration.
    Scalar dwell_fraction_max_;
    /// Injectivity multiplier oscillation threshold
    Scalar inj_mult_osc_threshold_;
    /// Injectivity multiplier dampening multiplier
    Scalar inj_mult_damp_mult_;
    /// Minimum damping factor for injectivity multipliers
    Scalar inj_mult_min_damp_factor_;
    /// Absolute max limit for residuals.
    Scalar max_residual_allowed_;
    //// Max allowed pore volume faction where CNV is violated. Below the
    //// relaxed tolerance tolerance_cnv_relaxed_ is used.
    Scalar relaxed_max_pv_fraction_;
    /// Relative mass balance tolerance (total mass balance error).
    Scalar tolerance_mb_;
    /// Relaxed mass balance tolerance (can be used when iter >= min_strict_mb_iter_).
    Scalar tolerance_mb_relaxed_;
    /// Relative energy balance tolerance (total energy balance error).
    Scalar tolerance_energy_balance_;
    /// Relaxed energy balance tolerance (can be used when iter >= min_strict_mb_iter_).
    Scalar tolerance_energy_balance_relaxed_;
    /// Local convergence tolerance (max of local saturation errors).
    Scalar tolerance_cnv_;
    /// Relaxed local convergence tolerance (can be used when iter >= min_strict_cnv_iter_ && cnvViolatedPV < relaxed_max_pv_fraction_).
    Scalar tolerance_cnv_relaxed_;
    /// Local energy convergence tolerance (max of local energy errors).
    Scalar tolerance_cnv_energy_;
    /// Relaxed local energy convergence tolerance (can be used when iter >= min_strict_cnv_iter_ && cnvViolatedPV < relaxed_max_pv_fraction_).
    Scalar tolerance_cnv_energy_relaxed_;
    /// Well convergence tolerance.
    Scalar tolerance_wells_;
    /// Tolerance for the well control equations
    //  TODO: it might need to distinguish between rate control and pressure control later
    Scalar tolerance_well_control_;
    /// Tolerance for the pressure equations for multisegment wells
    Scalar tolerance_pressure_ms_wells_;
    /// Relaxed tolerance for for the well flow residual
    Scalar relaxed_tolerance_flow_well_;

    /// Relaxed tolerance for the MSW pressure solution
    Scalar relaxed_tolerance_pressure_ms_well_;

    /// Maximum pressure change over an iteratio for ms wells
    Scalar max_pressure_change_ms_wells_;

    /// Maximum inner iteration number for ms wells
    int max_inner_iter_ms_wells_;

    /// Strict inner iteration number for wells
    int strict_inner_iter_wells_;

    /// Newton iteration where wells are stricly convergent
    int strict_outer_iter_wells_;

    /// Regularization factor for wells
    Scalar regularization_factor_wells_;

    /// Maximum newton iterations with inner well iterations
    int max_niter_inner_well_iter_;

    /// Whether to shut unsolvable well
    bool shut_unsolvable_wells_;

    /// Maximum inner iteration number for standard wells
    int max_inner_iter_wells_;

    /// Maximum iteration number of the well equation solution
    int max_welleq_iter_;

    /// Tolerance for time step in seconds where single precision can be used
    /// for solving for the Jacobian
    Scalar maxSinglePrecisionTimeStep_;

    /// Minimum number of Newton iterations before we can use relaxed CNV convergence criterion
    int min_strict_cnv_iter_;

    /// Minimum number of Newton iterations before we can use relaxed MB convergence criterion
    int min_strict_mb_iter_;

    /// Solve well equation initially
    bool solve_welleq_initially_;

    /// Pre solve and iterate network model
    bool pre_solve_network_;

    /// Update scaling factors for mass balance equations
    bool update_equations_scaling_;

    /// Try to detect oscillation or stagnation
    bool use_update_stabilization_;

    /// Whether to use MultisegmentWell to handle multisegment wells
    /// it is something temporary before the multisegment well model is considered to be
    /// well developed and tested.
    /// if it is false, we will handle multisegment wells as standard wells, which will be
    /// the default behavoir for the moment. Later, we might set it to be true by default if necessary
    bool use_multisegment_well_;

    /// The file name of the deck
    std::string deck_file_name_;

    /// Whether to add influences of wells between cells to the matrix and preconditioner matrix
    bool matrix_add_well_contributions_;

    /// Whether to check well operability
    bool check_well_operability_;
    /// Whether to check well operability during iterations
    bool check_well_operability_iter_;

    /// Maximum number of times a well can switch to the same control
    int max_number_of_well_switches_;

    /// Maximum number of times group can switch to the same control
    int max_number_of_group_switches_;

    /// Whether to approximate segment densities by averaging over segment and its outlet
    bool use_average_density_ms_wells_;

    /// Whether to allow control switching during local well solutions
    bool local_well_solver_control_switching_;

    /// Whether to use implicit IPR for thp stability checks and solution search
    bool use_implicit_ipr_;

    /// Whether to allow checking/changing to group controls during inner well iterations
    bool check_group_constraints_inner_well_iterations_; 

    /// Maximum number of iterations in the network solver before relaxing tolerance
    int network_max_strict_outer_iterations_;

    /// Maximum number of iterations in the network solver before giving up
    int network_max_outer_iterations_;

    /// Maximum number of sub-iterations to update network pressures (within a single well/group control update)
    int network_max_sub_iterations_;

    /// Damping factor in the inner network pressure update iterations
    Scalar network_pressure_update_damping_factor_;

    /// Maximum pressure update in the inner network pressure update iterations
    Scalar network_max_pressure_update_in_bars_;

    /// Maximum number of iterations in the well/group switch algorithm
    int well_group_constraints_max_iterations_;

    /// Maximum number of status switches (open<->shut> in local well iterations
    int max_well_status_switch_;

    /// Nonlinear solver type: newton or nldd
    std::string nonlinear_solver_;
  
    /// 'jacobi' and 'gauss-seidel' supported
    DomainSolveApproach local_solve_approach_{DomainSolveApproach::Jacobi};

    /// Maximum number of Newton iterations per time step
    int newton_max_iter_;

    /// Minimum number of Newton iterations per time step
    int newton_min_iter_;

    int max_local_solve_iterations_;

    Scalar local_tolerance_scaling_mb_;
    Scalar local_tolerance_scaling_cnv_;

    int nldd_num_initial_newton_iter_{1};
    /// Threshold for single cell relative mobility change in NLDD
    Scalar nldd_relative_mobility_change_tol_;
    int num_local_domains_{0};
    Scalar local_domains_partition_imbalance_{1.03};
    std::string local_domains_partition_method_;
    int local_domains_partition_well_neighbor_levels_{1};
    DomainOrderingMeasure local_domains_ordering_{DomainOrderingMeasure::MaxPressure};

    bool write_partitions_{false};

    /// Struct holding convergence monitor params
    struct ConvergenceMonitorParams
    {
        /// Whether to enable convergence monitoring
        bool enabled_;
        /// Cut-off limit for convergence monitoring
        int cutoff_;
        /// Decay factor used in convergence monitoring
        Scalar decay_factor_;
    };

    ConvergenceMonitorParams monitor_params_; //!< Convergence monitoring parameters

    // Relative tolerance of group rates (VREP, REIN)
    // If violated the nupcol wellstate is updated
    Scalar nupcol_group_rate_tolerance_;

    /// Construct from user parameters or defaults.
    BlackoilModelParameters();

    static void registerParameters();
};

} // namespace Opm

#endif // OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
