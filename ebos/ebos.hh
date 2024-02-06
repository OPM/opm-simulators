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

#include <ebos/eclproblem.hh>

#include <opm/models/utils/start.hh>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm {
template <class TypeTag>
class EbosProblem;
}

namespace Opm::Properties {

namespace TTag {
struct EbosTypeTag {
    using InheritsFrom = std::tuple<FlowModelParameters, EclBaseProblem, BlackOilModel, EclTimeSteppingParameters>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EbosTypeTag> {
    using type = EbosProblem<TypeTag>;
};

// Enable experimental features for ebos: ebos is the research simulator of the OPM
// project. If you're looking for a more stable "production quality" simulator, consider
// using `flow`
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = true;
};

// use flow's well model for now
template<class TypeTag>
struct WellModel<TypeTag, TTag::EbosTypeTag> {
    using type = BlackoilWellModel<TypeTag>;
};

// currently, ebos uses the non-multisegment well model by default to avoid
// regressions. the --use-multisegment-well=true|false command line parameter is still
// available in ebos, but hidden from view.
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = false;
};

// set some properties that are only required by the well model
template<class TypeTag>
struct MatrixAddWellContributions<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableTerminalOutput<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = false;
};

// flow's well model only works with surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = true;
};

// the values for the residual are for the whole cell instead of for a cubic meter of the cell
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = false;
};

// by default use flow's aquifer model for now
template<class TypeTag>
struct AquiferModel<TypeTag, TTag::EbosTypeTag> {
    using type = BlackoilAquiferModel<TypeTag>;
};

// use flow's linear solver backend for now
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::EbosTypeTag> {
    using type = TTag::FlowIstlSolver;
};

template<>
struct LinearSolverBackend<TTag::EbosTypeTag, TTag::FlowIstlSolverParams> {
    using type = ISTLSolver<TTag::EbosTypeTag>;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EbosTypeTag> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-1;
};

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class TypeTag>
struct EclNewtonRelaxedVolumeFraction<TypeTag, TTag::EbosTypeTag> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.05;
};

// the maximum volumetric error of a cell in the relaxed region
template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, TTag::EbosTypeTag> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e6*getPropValue<TypeTag, Properties::NewtonTolerance>();
};

// the tolerated amount of "incorrect" amount of oil per time step for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
template<class TypeTag>
struct EclNewtonSumTolerance<TypeTag, TTag::EbosTypeTag> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-5;
};

// make all Newton iterations strict, i.e., the volumetric Newton tolerance must be
// always be upheld in the majority of the spatial domain. In this context, "majority"
// means 1 - EclNewtonRelaxedVolumeFraction.
template<class TypeTag>
struct EclNewtonStrictIterations<TypeTag, TTag::EbosTypeTag> {
    static constexpr int value = 100;
};

// set the maximum number of Newton iterations to 8 so that we fail quickly (albeit
// relatively often)
template<class TypeTag>
struct NewtonMaxIterations<TypeTag, TTag::EbosTypeTag> {
    static constexpr int value = 8;
};

// if openMP is available, set the default the number of threads per process for the main
// simulation to 2 (instead of grabbing everything that is available).
#if _OPENMP
template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::EbosTypeTag> {
    static constexpr int value = 2;
};
#endif

// By default, ebos accepts the result of the time integration unconditionally if the
// smallest time step size is reached.
template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, TTag::EbosTypeTag> {
    static constexpr bool value = true;
};

} // namespace Opm::Properties

namespace Opm {
template <class TypeTag>
class EbosProblem : public EclProblem<TypeTag>
{
    typedef EclProblem<TypeTag> ParentType;

public:
    static void registerParameters()
    {
        ParentType::registerParameters();

        BlackoilModelParameters<TypeTag>::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTerminalOutput, "Do *NOT* use!");
        EWOMS_HIDE_PARAM(TypeTag, DbhpMaxRel);
        EWOMS_HIDE_PARAM(TypeTag, DwellFractionMax);
        EWOMS_HIDE_PARAM(TypeTag, MaxResidualAllowed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceMb);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceMbRelaxed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnv);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnvRelaxed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWells);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWellControl);
        EWOMS_HIDE_PARAM(TypeTag, MaxWelleqIter);
        EWOMS_HIDE_PARAM(TypeTag, UseMultisegmentWell);
        EWOMS_HIDE_PARAM(TypeTag, TolerancePressureMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxPressureChangeMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxInnerIterMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxNewtonIterationsWithInnerWellIterations);
        EWOMS_HIDE_PARAM(TypeTag, MaxInnerIterWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxSinglePrecisionDays);
        EWOMS_HIDE_PARAM(TypeTag, MinStrictCnvIter);
        EWOMS_HIDE_PARAM(TypeTag, MinStrictMbIter);
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
