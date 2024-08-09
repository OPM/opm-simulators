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
 * \brief The common settings for all flowexp variants.
 */
#ifndef FLOW_EXP_HPP
#define FLOW_EXP_HPP

#include <flowexperimental/FlowExpNewtonMethod.hpp>

#include <opm/models/discretization/common/fvbaseproblem.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/FlowProblemProperties.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm {
template <class TypeTag>
class FlowExpProblem;
}

namespace Opm::Properties {

namespace TTag {
struct FlowExpTypeTag {
    using InheritsFrom = std::tuple<FlowModelParameters, FlowBaseProblem, BlackOilModel, EclTimeSteppingParameters>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowExpTypeTag> {
    using type = FlowExpProblem<TypeTag>;
};

// Enable experimental features for flowexp: flowexp is the research simulator of the OPM
// project. If you're looking for a more stable "production quality" simulator, consider
// using `flow`
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::FlowExpTypeTag> {
    static constexpr bool value = true;
};

// use flow's well model for now
template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowExpTypeTag> {
    using type = BlackoilWellModel<TypeTag>;
};

template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::FlowExpTypeTag> {
    using type = FlowExpNewtonMethod<TypeTag>;
};

// flow's well model only works with surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::FlowExpTypeTag> {
    static constexpr bool value = true;
};

// the values for the residual are for the whole cell instead of for a cubic meter of the cell
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::FlowExpTypeTag> {
    static constexpr bool value = false;
};

// by default use flow's aquifer model for now
template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowExpTypeTag> {
    using type = BlackoilAquiferModel<TypeTag>;
};

// use flow's linear solver backend for now
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FlowExpTypeTag> {
    using type = TTag::FlowIstlSolver;
};

template<>
struct LinearSolverBackend<TTag::FlowExpTypeTag, TTag::FlowIstlSolverParams> {
    using type = ISTLSolver<TTag::FlowExpTypeTag>;
};

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowExpTypeTag> {
    using type = ISTLSolver<TypeTag>;
};

} // namespace Opm::Properties

namespace Opm::Parameters {

// if openMP is available, set the default the number of threads per process for the main
// simulation to 2 (instead of grabbing everything that is available).
#if _OPENMP
template<class TypeTag>
struct ThreadsPerProcess<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr int value = 2; };
#endif

// By default, flowexp accepts the result of the time integration unconditionally if the
// smallest time step size is reached.
template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableTerminalOutput<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr bool value = false; };

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, Properties::TTag::FlowExpTypeTag>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1e-1;
};

// set the maximum number of Newton iterations to 8 so that we fail quickly (albeit
// relatively often)
template<class TypeTag>
struct NewtonMaxIterations<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr int value = 8; };

// the maximum volumetric error of a cell in the relaxed region
template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, Properties::TTag::FlowExpTypeTag>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr auto baseValue =
        Parameters::NewtonTolerance<TypeTag,
                                    Properties::TTag::FlowExpTypeTag>::value;
    static constexpr type value = 1e6 * baseValue;
};

// currently, flowexp uses the non-multisegment well model by default to avoid
// regressions. the --use-multisegment-well=true|false command line parameter is still
// available in flowexp, but hidden from view.
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr bool value = false; };

// set some properties that are only required by the well model
template<class TypeTag>
struct MatrixAddWellContributions<TypeTag, Properties::TTag::FlowExpTypeTag>
{ static constexpr bool value = true; };

} // namespace Opm::Parameters

namespace Opm {
template <class TypeTag>
class FlowExpProblem : public FlowProblem<TypeTag> //, public FvBaseProblem<TypeTag>
{
    typedef FlowProblem<TypeTag> ParentType;
    using BaseType = ParentType; // GetPropType<TypeTag, Properties::BaseProblem>;
public:
    void writeOutput(bool verbose = true)
    {
        OPM_TIMEBLOCK(problemWriteOutput);
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        if (Parameters::get<TypeTag, Parameters::EnableWriteAllSolutions>() ||
            this->simulator().episodeWillBeOver())
        {
            // \Note: the SimulatorTimer does not carry any useful information, so PRT file (if it gets output) will contain wrong
            // timing information.
            BaseType::writeOutput(SimulatorTimer{}, verbose);
        }
    }

    static void registerParameters()
    {
        ParentType::registerParameters();

        BlackoilModelParameters<TypeTag>::registerParameters();
        Parameters::registerParam<TypeTag, Parameters::EnableTerminalOutput>("Do *NOT* use!");
        Parameters::hideParam<TypeTag, Parameters::DbhpMaxRel>();
        Parameters::hideParam<TypeTag, Parameters::DwellFractionMax>();
        Parameters::hideParam<TypeTag, Parameters::MaxResidualAllowed>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceMb>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceMbRelaxed>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceCnv>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceCnvRelaxed>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceWells>();
        Parameters::hideParam<TypeTag, Parameters::ToleranceWellControl>();
        Parameters::hideParam<TypeTag, Parameters::MaxWelleqIter>();
        Parameters::hideParam<TypeTag, Parameters::UseMultisegmentWell>();
        Parameters::hideParam<TypeTag, Parameters::TolerancePressureMsWells>();
        Parameters::hideParam<TypeTag, Parameters::MaxPressureChangeMsWells>();
        Parameters::hideParam<TypeTag, Parameters::MaxInnerIterMsWells>();
        Parameters::hideParam<TypeTag, Parameters::MaxNewtonIterationsWithInnerWellIterations>();
        Parameters::hideParam<TypeTag, Parameters::MaxInnerIterWells>();
        Parameters::hideParam<TypeTag, Parameters::MaxSinglePrecisionDays>();
        Parameters::hideParam<TypeTag, Parameters::MinStrictCnvIter>();
        Parameters::hideParam<TypeTag, Parameters::MinStrictMbIter>();
        Parameters::hideParam<TypeTag, Parameters::SolveWelleqInitially>();
        Parameters::hideParam<TypeTag, Parameters::UpdateEquationsScaling>();
        Parameters::hideParam<TypeTag, Parameters::UseUpdateStabilization>();
        Parameters::hideParam<TypeTag, Parameters::MatrixAddWellContributions>();
        Parameters::hideParam<TypeTag, Parameters::EnableTerminalOutput>();
    }

    // inherit the constructors
    using ParentType::FlowProblem;
};
}

#endif
