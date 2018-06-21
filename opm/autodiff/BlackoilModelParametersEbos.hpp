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

#ifndef OPM_BLACKOILMODELPARAMETERS_EBOS_HEADER_INCLUDED
#define OPM_BLACKOILMODELPARAMETERS_EBOS_HEADER_INCLUDED

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <string>

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowModelParameters);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclDeckFileName);

NEW_PROP_TAG(FlowDpMaxRel);
NEW_PROP_TAG(FlowDsMax);
NEW_PROP_TAG(FlowDrMaxRel);
NEW_PROP_TAG(FlowDbphMaxRel);
NEW_PROP_TAG(FlowDWellFractionMax);
NEW_PROP_TAG(FlowMaxResidualAllowed);
NEW_PROP_TAG(FlowToleranceMb);
NEW_PROP_TAG(FlowToleranceCnv);
NEW_PROP_TAG(FlowToleranceCnvRelaxed);
NEW_PROP_TAG(FlowToleranceWells);
NEW_PROP_TAG(FlowToleranceWellControl);
NEW_PROP_TAG(FlowMaxWelleqIter);
NEW_PROP_TAG(FlowUseMultisegmentWell);
NEW_PROP_TAG(FlowMaxSinglePrecisionDays);
NEW_PROP_TAG(FlowMaxStrictIter);
NEW_PROP_TAG(FlowSolveWelleqInitially);
NEW_PROP_TAG(FlowUpdateEquationsScaling);
NEW_PROP_TAG(FlowUseUpdateStabilization);
NEW_PROP_TAG(FlowMatrixAddWellContributions);
NEW_PROP_TAG(FlowPreconditionerAddWellContributions);

// parameters for multisegment wells
NEW_PROP_TAG(FlowTolerancePressureMsWells);
NEW_PROP_TAG(FlowMaxPressureChangeMsWells);
NEW_PROP_TAG(FlowUseInnerIterationsMsWells);
NEW_PROP_TAG(FlowMaxInnerIterMsWells);

SET_SCALAR_PROP(FlowModelParameters, FlowDpMaxRel, 0.3);
SET_SCALAR_PROP(FlowModelParameters, FlowDsMax, 0.2);
SET_SCALAR_PROP(FlowModelParameters, FlowDrMaxRel, 1e9);
SET_SCALAR_PROP(FlowModelParameters, FlowDbphMaxRel, 1.0);
SET_SCALAR_PROP(FlowModelParameters, FlowDWellFractionMax, 0.2);
SET_SCALAR_PROP(FlowModelParameters, FlowMaxResidualAllowed, 1e7);
SET_SCALAR_PROP(FlowModelParameters, FlowToleranceMb, 1e-5);
SET_SCALAR_PROP(FlowModelParameters, FlowToleranceCnv,1e-2);
SET_SCALAR_PROP(FlowModelParameters, FlowToleranceCnvRelaxed, 1e9);
SET_SCALAR_PROP(FlowModelParameters, FlowToleranceWells, 1e-4);
SET_SCALAR_PROP(FlowModelParameters, FlowToleranceWellControl, 1e-7);
SET_INT_PROP(FlowModelParameters, FlowMaxWelleqIter, 15);
SET_BOOL_PROP(FlowModelParameters, FlowUseMultisegmentWell, false);
SET_SCALAR_PROP(FlowModelParameters, FlowMaxSinglePrecisionDays, 20.0);
SET_INT_PROP(FlowModelParameters, FlowMaxStrictIter, 8);
SET_BOOL_PROP(FlowModelParameters, FlowSolveWelleqInitially, true);
SET_BOOL_PROP(FlowModelParameters, FlowUpdateEquationsScaling, false);
SET_BOOL_PROP(FlowModelParameters, FlowUseUpdateStabilization, true);
SET_BOOL_PROP(FlowModelParameters, FlowMatrixAddWellContributions, false);
SET_BOOL_PROP(FlowModelParameters, FlowPreconditionerAddWellContributions, false);
SET_SCALAR_PROP(FlowModelParameters, FlowTolerancePressureMsWells, 0.01 *1e5);
SET_SCALAR_PROP(FlowModelParameters, FlowMaxPressureChangeMsWells, 2.0 *1e5);
SET_BOOL_PROP(FlowModelParameters, FlowUseInnerIterationsMsWells, true);
SET_INT_PROP(FlowModelParameters, FlowMaxInnerIterMsWells, 10);

END_PROPERTIES

namespace Opm
{
    /// Solver parameters for the BlackoilModel.
    template <class TypeTag>
    struct BlackoilModelParametersEbos
    {
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    public:
        /// Max relative change in pressure in single iteration.
        double dp_max_rel_;
        /// Max absolute change in saturation in single iteration.
        double ds_max_;
        /// Max relative change in gas-oil or oil-gas ratio in single iteration.
        double dr_max_rel_;
        /// Max relative change in bhp in single iteration.
        double dbhp_max_rel_;
        /// Max absolute change in well volume fraction in single iteration.
        double dwell_fraction_max_;
        /// Absolute max limit for residuals.
        double max_residual_allowed_;
        /// Relative mass balance tolerance (total mass balance error).
        double tolerance_mb_;
        /// Local convergence tolerance (max of local saturation errors).
        double tolerance_cnv_;
        /// Relaxed local convergence tolerance (used when iter >= max_strict_iter_).
        double tolerance_cnv_relaxed_;
        /// Well convergence tolerance.
        double tolerance_wells_;
        /// Tolerance for the well control equations
        //  TODO: it might need to distinguish between rate control and pressure control later
        double tolerance_well_control_;
        /// Tolerance for the pressure equations for multisegment wells
        double tolerance_pressure_ms_wells_;
        /// Maximum pressure change over an iteratio for ms wells
        double max_pressure_change_ms_wells_;

        /// Whether to use inner iterations for ms wells
        bool use_inner_iterations_ms_wells_;

        /// Maximum inner iteration number for ms wells
        int max_inner_iter_ms_wells_;

        /// Maximum iteration number of the well equation solution
        int max_welleq_iter_;

        /// Tolerance for time step in seconds where single precision can be used
        /// for solving for the Jacobian
        double maxSinglePrecisionTimeStep_;

        /// Maximum number of Newton iterations before we give up on the CNV convergence criterion
        int max_strict_iter_;

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

        // Whether to add influences of wells between cells to the matrix and preconditioner matrix
        bool matrix_add_well_contributions_;

        // Whether to add influences of wells between cells to the preconditioner matrix only
        bool preconditioner_add_well_contributions_;

        /// Construct from user parameters or defaults.
        BlackoilModelParametersEbos()
        {
            dp_max_rel_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowDpMaxRel);
            ds_max_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowDsMax);
            dr_max_rel_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowDrMaxRel);
            dbhp_max_rel_=  EWOMS_GET_PARAM(TypeTag, Scalar, FlowDbphMaxRel);
            dwell_fraction_max_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowDWellFractionMax);
            max_residual_allowed_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowMaxResidualAllowed);
            tolerance_mb_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowToleranceMb);
            tolerance_cnv_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowToleranceCnv);
            tolerance_cnv_relaxed_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowToleranceCnvRelaxed);
            tolerance_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowToleranceWells);
            tolerance_well_control_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowToleranceWellControl);
            max_welleq_iter_ = EWOMS_GET_PARAM(TypeTag, int, FlowMaxWelleqIter);
            use_multisegment_well_ = EWOMS_GET_PARAM(TypeTag, bool, FlowUseMultisegmentWell);
            tolerance_pressure_ms_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowTolerancePressureMsWells);
            max_pressure_change_ms_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowMaxPressureChangeMsWells);
            use_inner_iterations_ms_wells_ = EWOMS_GET_PARAM(TypeTag, bool, FlowUseInnerIterationsMsWells);
            max_inner_iter_ms_wells_ = EWOMS_GET_PARAM(TypeTag, int, FlowMaxInnerIterMsWells);
            maxSinglePrecisionTimeStep_ = EWOMS_GET_PARAM(TypeTag, Scalar, FlowMaxSinglePrecisionDays) *24*60*60;
            max_strict_iter_ = EWOMS_GET_PARAM(TypeTag, int, FlowMaxStrictIter);
            solve_welleq_initially_ = EWOMS_GET_PARAM(TypeTag, bool, FlowSolveWelleqInitially);
            update_equations_scaling_ = EWOMS_GET_PARAM(TypeTag, bool, FlowUpdateEquationsScaling);
            use_update_stabilization_ = EWOMS_GET_PARAM(TypeTag, bool, FlowUseUpdateStabilization);
            matrix_add_well_contributions_ = EWOMS_GET_PARAM(TypeTag, bool, FlowMatrixAddWellContributions);
            preconditioner_add_well_contributions_ = EWOMS_GET_PARAM(TypeTag, bool, FlowPreconditionerAddWellContributions);

            deck_file_name_ = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        }

        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowDpMaxRel, "Maximum relative change of pressure in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowDsMax, "Maximum absolute change of any saturation in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowDrMaxRel, "Maximum relative change of the gas-in-oil or oil-in-gas ratio in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowDbphMaxRel, "Maximum relative change of the bottom-hole pressure in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowDWellFractionMax, "Maximum absolute change of a well's volume fraction in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowMaxResidualAllowed, "Absolute maximum tolerated for residuals without cutting the time step size");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowToleranceMb, "Tolerated mass balance error relative to total mass present");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowToleranceCnv, "Local convergence tolerance (Maximum of local saturation errors)");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowToleranceCnvRelaxed, "Relaxed local convergence tolerance that applies for iterations after the iterations with the strict tolerance");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowToleranceWells, "Well convergence tolerance");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowToleranceWellControl, "Tolerance for the well control equations");
            EWOMS_REGISTER_PARAM(TypeTag, int, FlowMaxWelleqIter, "Maximum number of iterations to determine solution the  well equations");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowUseMultisegmentWell, "Use the well model for multi-segment wells instead of the one for single-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowTolerancePressureMsWells, "Tolerance for the pressure equations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowMaxPressureChangeMsWells, "Maximum relative pressure change for a single iteration of the multi-segment well model");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowUseInnerIterationsMsWells, "Use nested iterations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, int, FlowMaxInnerIterMsWells, "Maximum number of inner iterations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlowMaxSinglePrecisionDays, "Maximum time step size where single precision floating point arithmetic can be used solving for the linear systems of equations");
            EWOMS_REGISTER_PARAM(TypeTag, int, FlowMaxStrictIter, "Maximum number of Newton iterations before relaxed tolerances are used for the CNV convergence criterion");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowSolveWelleqInitially, "Fully solve the well equations before each iteration of the reservoir model");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowUpdateEquationsScaling, "Update scaling factors for mass balance equations during the run");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowUseUpdateStabilization, "Try to detect and correct oscillations or stagnation during the Newton method");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowMatrixAddWellContributions, "Explicitly specify the influences of wells between cells in the Jacobian and preconditioner matrices");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FlowPreconditionerAddWellContributions, "Explicitly specify the influences of wells between cells for the preconditioner matrix only");
        }
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELPARAMETERS_EBOS_HEADER_INCLUDED
