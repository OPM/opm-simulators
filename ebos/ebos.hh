// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief The common settings for all ebos variants.
 */
#ifndef EBOS_HH
#define EBOS_HH

#include "eclproblem.hh"

#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/linalg/ISTLSolverEbos.hpp>

#include <ewoms/common/start.hh>

namespace Ewoms {
template <class TypeTag>
class EbosProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EbosTypeTag, INHERITS_FROM(BlackOilModel, EclBaseProblem, FlowModelParameters));

// Set the problem class
SET_TYPE_PROP(EbosTypeTag, Problem, Ewoms::EbosProblem<TypeTag>);

// Enable experimental features for ebos: ebos is the research simulator of the OPM
// project. If you're looking for a more stable "production quality" simulator, consider
// using `flow`
SET_BOOL_PROP(EbosTypeTag, EnableExperiments, true);

// use flow's well model for now
SET_TYPE_PROP(EbosTypeTag, EclWellModel, Opm::BlackoilWellModel<TypeTag>);

// currently, ebos uses the non-multisegment well model by default to avoid
// regressions. the --use-multisegment-well=true|false command line parameter is still
// available in ebos, but hidden from view.
SET_BOOL_PROP(EbosTypeTag, UseMultisegmentWell, false);

// set some properties that are only required by the well model
SET_BOOL_PROP(EbosTypeTag, MatrixAddWellContributions, true);
SET_BOOL_PROP(EbosTypeTag, EnableTerminalOutput, false);
// flow's well model only works with surface volumes
SET_BOOL_PROP(EbosTypeTag, BlackoilConserveSurfaceVolume, true);
// the values for the residual are for the whole cell instead of for a cubic meter of the cell
SET_BOOL_PROP(EbosTypeTag, UseVolumetricResidual, false);

// by default use flow's aquifer model for now
SET_TYPE_PROP(EbosTypeTag, EclAquiferModel, Opm::BlackoilAquiferModel<TypeTag>);

// use flow's linear solver backend for now
SET_TAG_PROP(EbosTypeTag, LinearSolverSplice, FlowIstlSolver);

// the default for the allowed volumetric error for oil per second
SET_SCALAR_PROP(EbosTypeTag, NewtonTolerance, 1e-1);

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
SET_SCALAR_PROP(EbosTypeTag, EclNewtonRelaxedVolumeFraction, 0.05);

// the maximum volumetric error of a cell in the relaxed region
SET_SCALAR_PROP(EbosTypeTag, EclNewtonRelaxedTolerance, 1e6*GET_PROP_VALUE(TypeTag, NewtonTolerance));

// the tolerated amount of "incorrect" amount of oil per time step for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
SET_SCALAR_PROP(EbosTypeTag, EclNewtonSumTolerance, 1e-5);

// make all Newton iterations strict, i.e., the volumetric Newton tolerance must be
// always be upheld in the majority of the spatial domain. In this context, "majority"
// means 1 - EclNewtonRelaxedVolumeFraction.
SET_INT_PROP(EbosTypeTag, EclNewtonStrictIterations, 100);

// set the maximum number of Newton iterations to 8 so that we fail quickly (albeit
// relatively often)
SET_INT_PROP(EbosTypeTag, NewtonMaxIterations, 8);

END_PROPERTIES

namespace Ewoms {
template <class TypeTag>
class EbosProblem : public EclProblem<TypeTag>
{
    typedef EclProblem<TypeTag> ParentType;

public:
    static void registerParameters()
    {
        ParentType::registerParameters();

        Opm::BlackoilModelParametersEbos<TypeTag>::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTerminalOutput, "Do *NOT* use!");
        EWOMS_HIDE_PARAM(TypeTag, DbhpMaxRel);
        EWOMS_HIDE_PARAM(TypeTag, DwellFractionMax);
        EWOMS_HIDE_PARAM(TypeTag, MaxResidualAllowed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceMb);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnv);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnvRelaxed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWells);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWellControl);
        EWOMS_HIDE_PARAM(TypeTag, MaxWelleqIter);
        EWOMS_HIDE_PARAM(TypeTag, UseMultisegmentWell);
        EWOMS_HIDE_PARAM(TypeTag, TolerancePressureMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxPressureChangeMsWells);
        EWOMS_HIDE_PARAM(TypeTag, UseInnerIterationsMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxInnerIterMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxSinglePrecisionDays);
        EWOMS_HIDE_PARAM(TypeTag, MaxStrictIter);
        EWOMS_HIDE_PARAM(TypeTag, SolveWelleqInitially);
        EWOMS_HIDE_PARAM(TypeTag, UpdateEquationsScaling);
        EWOMS_HIDE_PARAM(TypeTag, UseUpdateStabilization);
        EWOMS_HIDE_PARAM(TypeTag, MatrixAddWellContributions);
        EWOMS_HIDE_PARAM(TypeTag, EnableTerminalOutput);
    }

    // inherit the constructors
    using ParentType::EclProblem;
};
}

#endif // EBOS_HH
