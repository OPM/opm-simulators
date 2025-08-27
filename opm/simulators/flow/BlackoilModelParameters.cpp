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

#include <config.h>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/nonlinear/newtonmethodparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <algorithm>
#include <stdexcept>

namespace Opm {

template<class Scalar>
BlackoilModelParameters<Scalar>::BlackoilModelParameters()
{
    dbhp_max_rel_ = Parameters::Get<Parameters::DbhpMaxRel<Scalar>>();
    dwell_fraction_max_ = Parameters::Get<Parameters::DwellFractionMax<Scalar>>();
    inj_mult_osc_threshold_ = Parameters::Get<Parameters::InjMultOscThreshold<Scalar>>();
    inj_mult_damp_mult_ = Parameters::Get<Parameters::InjMultDampMult<Scalar>>();
    inj_mult_min_damp_factor_ = Parameters::Get<Parameters::InjMultMinDampFactor<Scalar>>();
    max_residual_allowed_ = Parameters::Get<Parameters::MaxResidualAllowed<Scalar>>();
    relaxed_max_pv_fraction_ = Parameters::Get<Parameters::RelaxedMaxPvFraction<Scalar>>();
    tolerance_mb_ = Parameters::Get<Parameters::ToleranceMb<Scalar>>();
    tolerance_mb_relaxed_ = std::max(tolerance_mb_, Parameters::Get<Parameters::ToleranceMbRelaxed<Scalar>>());
    tolerance_energy_balance_ = Parameters::Get<Parameters::ToleranceEnergyBalance<Scalar>>();
    tolerance_energy_balance_relaxed_ = std::max(tolerance_energy_balance_, Parameters::Get<Parameters::ToleranceEnergyBalanceRelaxed<Scalar>>());
    tolerance_cnv_ = Parameters::Get<Parameters::ToleranceCnv<Scalar>>();
    tolerance_cnv_relaxed_ = std::max(tolerance_cnv_, Parameters::Get<Parameters::ToleranceCnvRelaxed<Scalar>>());
    tolerance_cnv_energy_ = Parameters::Get<Parameters::ToleranceCnvEnergy<Scalar>>();
    tolerance_cnv_energy_relaxed_ = std::max(tolerance_cnv_energy_, Parameters::Get<Parameters::ToleranceCnvEnergyRelaxed<Scalar>>());
    tolerance_wells_ = Parameters::Get<Parameters::ToleranceWells<Scalar>>();
    tolerance_well_control_ = Parameters::Get<Parameters::ToleranceWellControl<Scalar>>();
    max_welleq_iter_ = Parameters::Get<Parameters::MaxWelleqIter>();
    use_multisegment_well_ = Parameters::Get<Parameters::UseMultisegmentWell>();
    tolerance_pressure_ms_wells_ = Parameters::Get<Parameters::TolerancePressureMsWells<Scalar>>();
    relaxed_tolerance_flow_well_ = Parameters::Get<Parameters::RelaxedWellFlowTol<Scalar>>();
    relaxed_tolerance_pressure_ms_well_ = Parameters::Get<Parameters::RelaxedPressureTolMsw<Scalar>>();
    max_pressure_change_ms_wells_ = Parameters::Get<Parameters::MaxPressureChangeMsWells<Scalar>>();
    max_inner_iter_ms_wells_ = Parameters::Get<Parameters::MaxInnerIterMsWells>();
    strict_inner_iter_wells_ = Parameters::Get<Parameters::StrictInnerIterWells>();
    strict_outer_iter_wells_ = Parameters::Get<Parameters::StrictOuterIterWells>();
    regularization_factor_wells_ = Parameters::Get<Parameters::RegularizationFactorWells<Scalar>>();
    max_niter_inner_well_iter_ = Parameters::Get<Parameters::MaxNewtonIterationsWithInnerWellIterations>();
    shut_unsolvable_wells_ = Parameters::Get<Parameters::ShutUnsolvableWells>();
    max_inner_iter_wells_ = Parameters::Get<Parameters::MaxInnerIterWells>();
    max_well_status_switch_ = Parameters::Get<Parameters::MaxWellStatusSwitchInInnerIterWells>();
    maxSinglePrecisionTimeStep_ = Parameters::Get<Parameters::MaxSinglePrecisionDays<Scalar>>() * 24 * 60 * 60;
    min_strict_cnv_iter_ = Parameters::Get<Parameters::MinStrictCnvIter>();
    min_strict_mb_iter_ = Parameters::Get<Parameters::MinStrictMbIter>();
    solve_welleq_initially_ = Parameters::Get<Parameters::SolveWelleqInitially>();
    pre_solve_network_ = Parameters::Get<Parameters::PreSolveNetwork>();
    update_equations_scaling_ = Parameters::Get<Parameters::UpdateEquationsScaling>();
    use_update_stabilization_ = Parameters::Get<Parameters::UseUpdateStabilization>();
    matrix_add_well_contributions_ = Parameters::Get<Parameters::MatrixAddWellContributions>();
    check_well_operability_ = Parameters::Get<Parameters::EnableWellOperabilityCheck>();
    check_well_operability_iter_ = Parameters::Get<Parameters::EnableWellOperabilityCheckIter>();
    max_number_of_well_switches_ = Parameters::Get<Parameters::MaximumNumberOfWellSwitches>();
    max_number_of_group_switches_ = Parameters::Get<Parameters::MaximumNumberOfGroupSwitches>();
    use_average_density_ms_wells_ = Parameters::Get<Parameters::UseAverageDensityMsWells>();
    local_well_solver_control_switching_ = Parameters::Get<Parameters::LocalWellSolveControlSwitching>();
    use_implicit_ipr_ = Parameters::Get<Parameters::UseImplicitIpr>();
    check_group_constraints_inner_well_iterations_ = Parameters::Get<Parameters::CheckGroupConstraintsInnerWellIterations>();
    nonlinear_solver_ = Parameters::Get<Parameters::NonlinearSolver>();
    const auto approach = Parameters::Get<Parameters::LocalSolveApproach>();
    if (approach == "jacobi") {
        local_solve_approach_ = DomainSolveApproach::Jacobi;
    } else if (approach == "gauss-seidel") {
        local_solve_approach_ = DomainSolveApproach::GaussSeidel;
    } else {
        throw std::runtime_error("Invalid domain solver approach '" + approach + "' specified.");
    }

    max_local_solve_iterations_ = Parameters::Get<Parameters::MaxLocalSolveIterations>();
    local_tolerance_scaling_mb_ = Parameters::Get<Parameters::LocalToleranceScalingMb<Scalar>>();
    local_tolerance_scaling_cnv_ = Parameters::Get<Parameters::LocalToleranceScalingCnv<Scalar>>();
    newton_max_iter_ = Parameters::Get<Parameters::NewtonMaxIterations>();
    newton_min_iter_ = Parameters::Get<Parameters::NewtonMinIterations>();
    nldd_num_initial_newton_iter_ = Parameters::Get<Parameters::NlddNumInitialNewtonIter>();
    nldd_relative_mobility_change_tol_ = Parameters::Get<Parameters::NlddRelativeMobilityChangeTol<Scalar>>();
    num_local_domains_ = Parameters::Get<Parameters::NumLocalDomains>();
    local_domains_partition_imbalance_ = std::max(Scalar{1.0}, Parameters::Get<Parameters::LocalDomainsPartitioningImbalance<Scalar>>());
    local_domains_partition_method_ = Parameters::Get<Parameters::LocalDomainsPartitioningMethod>();
    local_domains_partition_well_neighbor_levels_ = Parameters::Get<Parameters::LocalDomainsPartitionWellNeighborLevels>();
    deck_file_name_ = Parameters::Get<Parameters::EclDeckFileName>();
    network_max_strict_outer_iterations_ = Parameters::Get<Parameters::NetworkMaxStrictOuterIterations>();
    network_max_outer_iterations_ = Parameters::Get<Parameters::NetworkMaxOuterIterations>();
    network_max_sub_iterations_ = Parameters::Get<Parameters::NetworkMaxSubIterations>();
    network_pressure_update_damping_factor_ = Parameters::Get<Parameters::NetworkPressureUpdateDampingFactor<Scalar>>();
    network_max_pressure_update_in_bars_ = Parameters::Get<Parameters::NetworkMaxPressureUpdateInBars<Scalar>>();
    local_domains_ordering_ = domainOrderingMeasureFromString(Parameters::Get<Parameters::LocalDomainsOrderingMeasure>());
    write_partitions_ = Parameters::Get<Parameters::DebugEmitCellPartition>();

    monitor_params_.enabled_ = Parameters::Get<Parameters::ConvergenceMonitoring>();
    monitor_params_.cutoff_ = Parameters::Get<Parameters::ConvergenceMonitoringCutOff>();
    monitor_params_.decay_factor_ = Parameters::Get<Parameters::ConvergenceMonitoringDecayFactor<Scalar>>();

    nupcol_group_rate_tolerance_ = Parameters::Get<Parameters::NupcolGroupRateTolerance<Scalar>>();
    well_group_constraints_max_iterations_ = Parameters::Get<Parameters::WellGroupConstraintsMaxIterations>();
}

template<class Scalar>
void BlackoilModelParameters<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::DbhpMaxRel<Scalar>>
        ("Maximum relative change of the bottom-hole pressure in a single iteration");
    Parameters::Register<Parameters::DwellFractionMax<Scalar>>
        ("Maximum absolute change of a well's volume fraction in a single iteration");
    Parameters::Register<Parameters::InjMultOscThreshold<Scalar>>
        ("Injection multiplier oscillation threshold (used for multiplier dampening)");
    Parameters::Register<Parameters::InjMultDampMult<Scalar>>
        ("Injection multiplier dampening factor (dampening multiplied by this each time oscillation is detected)");
    Parameters::Register<Parameters::InjMultMinDampFactor<Scalar>>
        ("Minimum injection multiplier dampening factor (maximum dampening level)");
    Parameters::Register<Parameters::MaxResidualAllowed<Scalar>>
        ("Absolute maximum tolerated for residuals without cutting the time step size");
    Parameters::Register<Parameters::RelaxedMaxPvFraction<Scalar>>
        ("The fraction of the pore volume of the reservoir "
         "where the volumetric error (CNV) may be violated "
         "during strict Newton iterations.");
    Parameters::Register<Parameters::ToleranceMb<Scalar>>
        ("Tolerated mass balance error relative to total mass present");
    Parameters::Register<Parameters::ToleranceMbRelaxed<Scalar>>
        ("Relaxed tolerated mass balance error that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceEnergyBalance<Scalar>>
        ("Tolerated energy balance error relative to (scaled) total energy present");
    Parameters::Register<Parameters::ToleranceEnergyBalanceRelaxed<Scalar>>
        ("Relaxed tolerated energy balance error that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceCnv<Scalar>>
        ("Local convergence tolerance (Maximum of local saturation errors)");
    Parameters::Register<Parameters::ToleranceCnvRelaxed<Scalar>>
        ("Relaxed local convergence tolerance that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceCnvEnergy<Scalar>>
        ("Local energy convergence tolerance (Maximum of local energy errors)");
    Parameters::Register<Parameters::ToleranceCnvEnergyRelaxed<Scalar>>
        ("Relaxed local energy convergence tolerance that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceWells<Scalar>>
        ("Well convergence tolerance");
    Parameters::Register<Parameters::ToleranceWellControl<Scalar>>
        ("Tolerance for the well control equations");
    Parameters::Register<Parameters::MaxWelleqIter>
        ("Maximum number of iterations to determine solution the well equations");
    Parameters::Register<Parameters::UseMultisegmentWell>
        ("Use the well model for multi-segment wells instead of the "
         "one for single-segment wells");
    Parameters::Register<Parameters::TolerancePressureMsWells<Scalar>>
        ("Tolerance for the pressure equations for multi-segment wells");
    Parameters::Register<Parameters::RelaxedWellFlowTol<Scalar>>
        ("Relaxed tolerance for the well flow residual");
    Parameters::Register<Parameters::RelaxedPressureTolMsw<Scalar>>
        ("Relaxed tolerance for the MSW pressure solution");
    Parameters::Register<Parameters::MaxPressureChangeMsWells<Scalar>>
        ("Maximum relative pressure change for a single iteration "
         "of the multi-segment well model");
    Parameters::Register<Parameters::MaxInnerIterMsWells>
        ("Maximum number of inner iterations for multi-segment wells");
    Parameters::Register<Parameters::StrictInnerIterWells>
        ("Number of inner well iterations with strict tolerance");
    Parameters::Register<Parameters::StrictOuterIterWells>
        ("Number of newton iterations for which wells are checked with strict tolerance");
    Parameters::Register<Parameters::MaxNewtonIterationsWithInnerWellIterations>
        ("Maximum newton iterations with inner well iterations");
    Parameters::Register<Parameters::ShutUnsolvableWells>
        ("Shut unsolvable wells");
    Parameters::Register<Parameters::MaxInnerIterWells>
        ("Maximum number of inner iterations for standard wells");
    Parameters::Register<Parameters::MaxWellStatusSwitchInInnerIterWells>
        ("Maximum number of status switching (shut<->open) in inner iterations for wells");
    Parameters::Register<Parameters::AlternativeWellRateInit>
        ("Use alternative well rate initialization procedure");
    Parameters::Register<Parameters::RegularizationFactorWells<Scalar>>
        ("Regularization factor for wells");
    Parameters::Register<Parameters::MaxSinglePrecisionDays<Scalar>>
        ("Maximum time step size where single precision floating point "
         "arithmetic can be used solving for the linear systems of equations");
    Parameters::Register<Parameters::MinStrictCnvIter>
        ("Minimum number of Newton iterations before relaxed tolerances "
         "can be used for the CNV convergence criterion");
    Parameters::Register<Parameters::MinStrictMbIter>
        ("Minimum number of Newton iterations before relaxed tolerances "
         "can be used for the MB convergence criterion. "
         "Default -1 means that the relaxed tolerance is used when maximum "
         "number of Newton iterations are reached.");
    Parameters::Register<Parameters::SolveWelleqInitially>
        ("Fully solve the well equations before each iteration of the reservoir model");
    Parameters::Register<Parameters::PreSolveNetwork>
        ("Pre solve and iterate the network model at start-up");
    Parameters::Register<Parameters::UpdateEquationsScaling>
        ("Update scaling factors for mass balance equations during the run");
    Parameters::Register<Parameters::UseUpdateStabilization>
        ("Try to detect and correct oscillations or stagnation during the Newton method");
    Parameters::Register<Parameters::MatrixAddWellContributions>
        ("Explicitly specify the influences of wells between cells in "
         "the Jacobian and preconditioner matrices");
    Parameters::Register<Parameters::EnableWellOperabilityCheck>
        ("Enable the well operability checking");
    Parameters::Register<Parameters::EnableWellOperabilityCheckIter>
        ("Enable the well operability checking during iterations");
    Parameters::Register<Parameters::MaximumNumberOfWellSwitches>
        ("Maximum number of times a well can switch to the same control");
    Parameters::Register<Parameters::MaximumNumberOfGroupSwitches>
        ("Maximum number of times a group can switch to the same control");
    Parameters::Register<Parameters::UseAverageDensityMsWells>
        ("Approximate segment densitities by averaging over segment and its outlet");
    Parameters::Register<Parameters::LocalWellSolveControlSwitching>
        ("Allow control switching during local well solutions");
    Parameters::Register<Parameters::UseImplicitIpr>
        ("Compute implict IPR for stability checks and stable solution search");
    Parameters::Register<Parameters::CheckGroupConstraintsInnerWellIterations>
        ("Allow checking of group constraints during inner well iterations");        
    Parameters::Register<Parameters::NetworkMaxStrictOuterIterations>
        ("Maximum outer iterations in network solver before relaxing tolerance");
    Parameters::Register<Parameters::NetworkMaxOuterIterations>
        ("Maximum outer number of iterations in the network solver before giving up");
    Parameters::Register<Parameters::NetworkMaxSubIterations>
        ("Maximum number of sub-iterations to update network pressures (within a single well/group control update)");
    Parameters::Register<Parameters::NetworkPressureUpdateDampingFactor<Scalar>>
        ("Damping factor in the inner network pressure update iterations");
    Parameters::Register<Parameters::NetworkMaxPressureUpdateInBars<Scalar>>
        ("Maximum pressure update in the inner network pressure update iterations");
    Parameters::Register<Parameters::NonlinearSolver>
        ("Choose nonlinear solver. Valid choices are newton or nldd.");
    Parameters::Register<Parameters::LocalSolveApproach>
        ("Choose local solve approach. Valid choices are jacobi and gauss-seidel");
    Parameters::SetDefault<Parameters::NewtonMaxIterations>(20);
    Parameters::Register<Parameters::NewtonMinIterations>
        ("The minimum number of Newton iterations per time step");
    Parameters::Register<Parameters::MaxLocalSolveIterations>
        ("Max iterations for local solves with NLDD nonlinear solver.");
    Parameters::Register<Parameters::LocalToleranceScalingMb<Scalar>>
        ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
    Parameters::Register<Parameters::LocalToleranceScalingCnv<Scalar>>
        ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
    Parameters::Register<Parameters::NlddNumInitialNewtonIter>
        ("Number of initial global Newton iterations when running the NLDD nonlinear solver.");
    Parameters::Register<Parameters::NlddRelativeMobilityChangeTol<Scalar>>
        ("Threshold for single cell relative mobility change in the NLDD solver");
    Parameters::Register<Parameters::NumLocalDomains>
        ("Number of local domains for NLDD nonlinear solver.");
    Parameters::Register<Parameters::LocalDomainsPartitioningImbalance<Scalar>>
        ("Subdomain partitioning imbalance tolerance. 1.03 is 3 percent imbalance.");
    Parameters::Register<Parameters::LocalDomainsPartitioningMethod>
        ("Subdomain partitioning method. Allowed values are "
         "'zoltan', "
         "'simple', "
         "and the name of a partition file ending with '.partition'.");
    Parameters::Register<Parameters::LocalDomainsPartitionWellNeighborLevels>
        ("Number of neighbor levels around wells to include in the same domain during NLDD partitioning");
    Parameters::Register<Parameters::LocalDomainsOrderingMeasure>
        ("Subdomain ordering measure. Allowed values are "
         "'maxpressure', "
         "'averagepressure' "
         "and  'residual'.");
    Parameters::Register<Parameters::DebugEmitCellPartition>
        ("Whether or not to emit cell partitions as a debugging aid.");

    Parameters::Register<Parameters::ConvergenceMonitoring>
        ("Enable convergence monitoring");
    Parameters::Register<Parameters::ConvergenceMonitoringCutOff>
        ("Cut off limit for convergence monitoring");
    Parameters::Register<Parameters::ConvergenceMonitoringDecayFactor<Scalar>>
        ("Decay factor for convergence monitoring");

    Parameters::Register<Parameters::NupcolGroupRateTolerance<Scalar>>
        ("Tolerance for acceptable changes in VREP/RAIN group rates");

    Parameters::Hide<Parameters::DebugEmitCellPartition>();

    Parameters::Register<Parameters::WellGroupConstraintsMaxIterations>
    ("Maximum number of iterations in the well/group switching algorithm");

    // if openMP is available, use two threads per mpi rank by default
#if _OPENMP
    Parameters::SetDefault<Parameters::ThreadsPerProcess>(2);
#endif
}

template struct BlackoilModelParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct BlackoilModelParameters<float>;
#endif

} // namespace Opm
