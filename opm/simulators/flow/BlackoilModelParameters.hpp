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

#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/flow/SubDomain.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace Opm::Parameters {

template<class Scalar>
struct DbhpMaxRel { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct DwellFractionMax { static constexpr Scalar value = 0.2; };

struct EclDeckFileName { static constexpr auto value = ""; };

template<class Scalar>
struct MaxResidualAllowed { static constexpr Scalar value = 1e7; };

template<class Scalar>
struct RelaxedMaxPvFraction { static constexpr Scalar value = 0.03; };

template<class Scalar>
struct ToleranceMb { static constexpr Scalar value = 1e-7; };

template<class Scalar>
struct ToleranceMbRelaxed { static constexpr Scalar value = 1e-6; };

template<class Scalar>
struct ToleranceCnv { static constexpr Scalar value = 1e-2; };

template<class Scalar>
struct ToleranceCnvRelaxed { static constexpr Scalar value = 1.0; };

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
struct UseAverageDensityMsWells { static constexpr bool value = false; };
struct LocalWellSolveControlSwitching { static constexpr bool value = true; };
struct UseImplicitIpr { static constexpr bool value = true; };

// Network solver parameters
struct NetworkMaxStrictIterations { static constexpr int value = 100; };
struct NetworkMaxIterations { static constexpr int value = 200; };
struct NonlinearSolver { static constexpr auto value = "newton"; };
struct LocalSolveApproach { static constexpr auto value = "gauss-seidel"; };
struct MaxLocalSolveIterations { static constexpr int value = 20; };

template<class Scalar>
struct LocalToleranceScalingMb { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct LocalToleranceScalingCnv { static constexpr Scalar value = 0.1; };
struct NlddNumInitialNewtonIter { static constexpr int value = 1; };
struct NumLocalDomains { static constexpr int value = 0; };

template<class Scalar>
struct LocalDomainsPartitioningImbalance { static constexpr Scalar value = 1.03; };

struct LocalDomainsPartitioningMethod { static constexpr auto value = "zoltan"; };
struct LocalDomainsOrderingMeasure { static constexpr auto value = "maxpressure"; };

} // namespace Opm::Parameters

namespace Opm {

/// Solver parameters for the BlackoilModel.
template <class TypeTag>
struct BlackoilModelParameters
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    /// Max relative change in bhp in single iteration.
    Scalar dbhp_max_rel_;
    /// Max absolute change in well volume fraction in single iteration.
    Scalar dwell_fraction_max_;
    /// Absolute max limit for residuals.
    Scalar max_residual_allowed_;
    //// Max allowed pore volume faction where CNV is violated. Below the
    //// relaxed tolerance tolerance_cnv_relaxed_ is used.
    Scalar relaxed_max_pv_fraction_;
    /// Relative mass balance tolerance (total mass balance error).
    Scalar tolerance_mb_;
    /// Relaxed mass balance tolerance (can be used when iter >= min_strict_mb_iter_).
    Scalar tolerance_mb_relaxed_;
    /// Local convergence tolerance (max of local saturation errors).
    Scalar tolerance_cnv_;
    /// Relaxed local convergence tolerance (can be used when iter >= min_strict_cnv_iter_ && cnvViolatedPV < relaxed_max_pv_fraction_).
    Scalar tolerance_cnv_relaxed_;
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

    /// Update scaling factors for mass balance equations
    bool update_equations_scaling_;

    /// Try to detect oscillation or stagnation.
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

    /// Maximum number of times a well can switch to the same controt
    int max_number_of_well_switches_;

    /// Whether to approximate segment densities by averaging over segment and its outlet
    bool use_average_density_ms_wells_;

    /// Whether to allow control switching during local well solutions
    bool local_well_solver_control_switching_;

    /// Whether to use implicit IPR for thp stability checks and solution search
    bool use_implicit_ipr_;

    /// Maximum number of iterations in the network solver before relaxing tolerance
    int network_max_strict_iterations_;

    /// Maximum number of iterations in the network solver before giving up
    int network_max_iterations_;

    /// Nonlinear solver type: newton or nldd.
    std::string nonlinear_solver_;
    /// 'jacobi' and 'gauss-seidel' supported.
    DomainSolveApproach local_solve_approach_{DomainSolveApproach::Jacobi};

    int max_local_solve_iterations_;

    Scalar local_tolerance_scaling_mb_;
    Scalar local_tolerance_scaling_cnv_;

    int nldd_num_initial_newton_iter_{1};
    int num_local_domains_{0};
    Scalar local_domain_partition_imbalance_{1.03};
    std::string local_domain_partition_method_;
    DomainOrderingMeasure local_domain_ordering_{DomainOrderingMeasure::MaxPressure};

    bool write_partitions_{false};

    /// Construct from user parameters or defaults.
    BlackoilModelParameters()
    {
        dbhp_max_rel_ = Parameters::Get<Parameters::DbhpMaxRel<Scalar>>();
        dwell_fraction_max_ = Parameters::Get<Parameters::DwellFractionMax<Scalar>>();
        max_residual_allowed_ = Parameters::Get<Parameters::MaxResidualAllowed<Scalar>>();
        relaxed_max_pv_fraction_ = Parameters::Get<Parameters::RelaxedMaxPvFraction<Scalar>>();
        tolerance_mb_ = Parameters::Get<Parameters::ToleranceMb<Scalar>>();
        tolerance_mb_relaxed_ = std::max(tolerance_mb_, Parameters::Get<Parameters::ToleranceMbRelaxed<Scalar>>());
        tolerance_cnv_ = Parameters::Get<Parameters::ToleranceCnv<Scalar>>();
        tolerance_cnv_relaxed_ = std::max(tolerance_cnv_, Parameters::Get<Parameters::ToleranceCnvRelaxed<Scalar>>());
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
        maxSinglePrecisionTimeStep_ = Parameters::Get<Parameters::MaxSinglePrecisionDays<Scalar>>() * 24 * 60 * 60;
        min_strict_cnv_iter_ = Parameters::Get<Parameters::MinStrictCnvIter>();
        min_strict_mb_iter_ = Parameters::Get<Parameters::MinStrictMbIter>();
        solve_welleq_initially_ = Parameters::Get<Parameters::SolveWelleqInitially>();
        update_equations_scaling_ = Parameters::Get<Parameters::UpdateEquationsScaling>();
        use_update_stabilization_ = Parameters::Get<Parameters::UseUpdateStabilization>();
        matrix_add_well_contributions_ = Parameters::Get<Parameters::MatrixAddWellContributions>();
        check_well_operability_ = Parameters::Get<Parameters::EnableWellOperabilityCheck>();
        check_well_operability_iter_ = Parameters::Get<Parameters::EnableWellOperabilityCheckIter>();
        max_number_of_well_switches_ = Parameters::Get<Parameters::MaximumNumberOfWellSwitches>();
        use_average_density_ms_wells_ = Parameters::Get<Parameters::UseAverageDensityMsWells>();
        local_well_solver_control_switching_ = Parameters::Get<Parameters::LocalWellSolveControlSwitching>();
        use_implicit_ipr_ = Parameters::Get<Parameters::UseImplicitIpr>();
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
        nldd_num_initial_newton_iter_ = Parameters::Get<Parameters::NlddNumInitialNewtonIter>();
        num_local_domains_ = Parameters::Get<Parameters::NumLocalDomains>();
        local_domain_partition_imbalance_ = std::max(Scalar{1.0}, Parameters::Get<Parameters::LocalDomainsPartitioningImbalance<Scalar>>());
        local_domain_partition_method_ = Parameters::Get<Parameters::LocalDomainsPartitioningMethod>();
        deck_file_name_ = Parameters::Get<Parameters::EclDeckFileName>();
        network_max_strict_iterations_ = Parameters::Get<Parameters::NetworkMaxStrictIterations>();
        network_max_iterations_ = Parameters::Get<Parameters::NetworkMaxIterations>();
        local_domain_ordering_ = domainOrderingMeasureFromString(Parameters::Get<Parameters::LocalDomainsOrderingMeasure>());
        write_partitions_ = Parameters::Get<Parameters::DebugEmitCellPartition>();
    }

    static void registerParameters()
    {
        Parameters::Register<Parameters::DbhpMaxRel<Scalar>>
            ("Maximum relative change of the bottom-hole pressure in a single iteration");
        Parameters::Register<Parameters::DwellFractionMax<Scalar>>
            ("Maximum absolute change of a well's volume fraction in a single iteration");
        Parameters::Register<Parameters::MaxResidualAllowed<Scalar>>
            ("Absolute maximum tolerated for residuals without cutting the time step size");
        Parameters::Register<Parameters::RelaxedMaxPvFraction<Scalar>>
            ("The fraction of the pore volume of the reservoir "
             "where the volumetric error (CNV) may be voilated "
             "during strict Newton iterations.");
        Parameters::Register<Parameters::ToleranceMb<Scalar>>
            ("Tolerated mass balance error relative to total mass present");
        Parameters::Register<Parameters::ToleranceMbRelaxed<Scalar>>
            ("Relaxed tolerated mass balance error that applies for iterations "
             "after the iterations with the strict tolerance");
        Parameters::Register<Parameters::ToleranceCnv<Scalar>>
            ("Local convergence tolerance (Maximum of local saturation errors)");
        Parameters::Register<Parameters::ToleranceCnvRelaxed<Scalar>>
            ("Relaxed local convergence tolerance that applies for iterations "
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
        Parameters::Register<Parameters::UseAverageDensityMsWells>
            ("Approximate segment densitities by averaging over segment and its outlet");
        Parameters::Register<Parameters::LocalWellSolveControlSwitching>
            ("Allow control switching during local well solutions");
        Parameters::Register<Parameters::UseImplicitIpr>
            ("Compute implict IPR for stability checks and stable solution search");
        Parameters::Register<Parameters::NetworkMaxStrictIterations>
            ("Maximum iterations in network solver before relaxing tolerance");
        Parameters::Register<Parameters::NetworkMaxIterations>
            ("Maximum number of iterations in the network solver before giving up");
        Parameters::Register<Parameters::NonlinearSolver>
            ("Choose nonlinear solver. Valid choices are newton or nldd.");
        Parameters::Register<Parameters::LocalSolveApproach>
            ("Choose local solve approach. Valid choices are jacobi and gauss-seidel");
        Parameters::Register<Parameters::MaxLocalSolveIterations>
            ("Max iterations for local solves with NLDD nonlinear solver.");
        Parameters::Register<Parameters::LocalToleranceScalingMb<Scalar>>
            ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
        Parameters::Register<Parameters::LocalToleranceScalingCnv<Scalar>>
            ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
        Parameters::Register<Parameters::NlddNumInitialNewtonIter>
            ("Number of initial global Newton iterations when running the NLDD nonlinear solver.");
        Parameters::Register<Parameters::NumLocalDomains>
            ("Number of local domains for NLDD nonlinear solver.");
        Parameters::Register<Parameters::LocalDomainsPartitioningImbalance<Scalar>>
            ("Subdomain partitioning imbalance tolerance. 1.03 is 3 percent imbalance.");
        Parameters::Register<Parameters::LocalDomainsPartitioningMethod>
            ("Subdomain partitioning method. Allowed values are "
             "'zoltan', "
             "'simple', "
             "and the name of a partition file ending with '.partition'.");
        Parameters::Register<Parameters::LocalDomainsOrderingMeasure>
            ("Subdomain ordering measure. Allowed values are "
             "'maxpressure', "
             "'averagepressure' "
             "and  'residual'.");
        Parameters::Register<Parameters::DebugEmitCellPartition>
            ("Whether or not to emit cell partitions as a debugging aid.");


        Parameters::Hide<Parameters::DebugEmitCellPartition>();

        // if openMP is available, determine the number threads per process automatically.
#if _OPENMP
        Parameters::SetDefault<Parameters::ThreadsPerProcess>(-1);
#endif
    }
};

} // namespace Opm

#endif // OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
