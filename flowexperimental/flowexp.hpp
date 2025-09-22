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

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm {
template <class TypeTag>
class FlowExpProblem;
}

namespace Opm::Properties {

namespace TTag {

struct FlowExpTypeTag
{
    using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>;
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

namespace Opm {

template <class TypeTag>
class FlowExpProblem : public FlowProblemBlackoil<TypeTag> //, public FvBaseProblem<TypeTag>
{
    using ParentType = FlowProblemBlackoil<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void writeOutput(bool verbose = true) override
    {
        OPM_TIMEBLOCK(problemWriteOutput);
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        if (Parameters::Get<Parameters::EnableWriteAllSolutions>() ||
            this->simulator().episodeWillBeOver())
        {
            // \Note: the SimulatorTimer does not carry any useful information, so PRT file (if it gets output) will contain wrong
            // timing information.
            ParentType::writeOutput(verbose);
        }
    }

    static void registerParameters()
    {
        ParentType::registerParameters();

        BlackoilModelParameters<double>::registerParameters();
        Parameters::Register<Parameters::EnableTerminalOutput>("Do *NOT* use!");
        Parameters::Hide<Parameters::DbhpMaxRel<Scalar>>();
        Parameters::Hide<Parameters::DwellFractionMax<Scalar>>();
        Parameters::Hide<Parameters::MaxResidualAllowed<Scalar>>();
        Parameters::Hide<Parameters::ToleranceMb<Scalar>>();
        Parameters::Hide<Parameters::ToleranceMbRelaxed<Scalar>>();
        Parameters::Hide<Parameters::ToleranceCnv<Scalar>>();
        Parameters::Hide<Parameters::ToleranceCnvRelaxed<Scalar>>();
        Parameters::Hide<Parameters::ToleranceWells<Scalar>>();
        Parameters::Hide<Parameters::ToleranceWellControl<Scalar>>();
        Parameters::Hide<Parameters::MaxWelleqIter>();
        Parameters::Hide<Parameters::UseMultisegmentWell>();
        Parameters::Hide<Parameters::TolerancePressureMsWells<Scalar>>();
        Parameters::Hide<Parameters::MaxPressureChangeMsWells<Scalar>>();
        Parameters::Hide<Parameters::MaxInnerIterMsWells>();
        Parameters::Hide<Parameters::MaxNewtonIterationsWithInnerWellIterations>();
        Parameters::Hide<Parameters::MaxInnerIterWells>();
        Parameters::Hide<Parameters::MaxWellStatusSwitchInInnerIterWells>();
        Parameters::Hide<Parameters::MaxWellStatusSwitchForWells>();
        Parameters::Hide<Parameters::MaxSinglePrecisionDays<Scalar>>();
        Parameters::Hide<Parameters::MinStrictCnvIter>();
        Parameters::Hide<Parameters::MinStrictMbIter>();
        Parameters::Hide<Parameters::SolveWelleqInitially>();
        Parameters::Hide<Parameters::PreSolveNetwork>();
        Parameters::Hide<Parameters::UpdateEquationsScaling>();
        Parameters::Hide<Parameters::UseUpdateStabilization>();
        Parameters::Hide<Parameters::MatrixAddWellContributions>();
        Parameters::Hide<Parameters::EnableTerminalOutput>();

        // if openMP is available, set the default the number of threads per process for the main
        // simulation to 1 (instead of grabbing everything that is available).
#if _OPENMP
        Parameters::SetDefault<Parameters::ThreadsPerProcess>(2);
#endif

        // By default, flowexp accepts the result of the time integration unconditionally if the
        // smallest time step size is reached.
        Parameters::SetDefault<Parameters::ContinueOnConvergenceError>(true);

        Parameters::SetDefault<Parameters::NewtonMaxIterations>(8);
        Parameters::SetDefault<Parameters::NewtonTolerance<Scalar>>(1e-2);
        Parameters::SetDefault<Parameters::EclNewtonRelaxedTolerance<Scalar>>(1e-1);
        Parameters::SetDefault<Parameters::EclNewtonRelaxedVolumeFraction<Scalar>>(0.0);
        Parameters::SetDefault<Parameters::EclNewtonSumTolerance<Scalar>>(1e-5);
        Parameters::SetDefault<Parameters::EnableTerminalOutput>(false);

        // currently, flowexp uses the non-multisegment well model by default to avoid
        // regressions. the --use-multisegment-well=true|false command line parameter is still
        // available in flowexp, but hidden from view.
        Parameters::SetDefault<Parameters::UseMultisegmentWell>(false);
        Parameters::SetDefault<Parameters::MatrixAddWellContributions>(false);
    }

    // inherit the constructors
    using ParentType::FlowProblemBlackoil;
};

}

#endif
