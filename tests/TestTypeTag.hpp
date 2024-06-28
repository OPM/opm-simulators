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
 * \brief Typetag used in tests
 */
#ifndef OPM_TEST_TYPETAG_HPP
#define OPM_TEST_TYPETAG_HPP

#include <opm/models/utils/start.hh>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/flow/FlowProblemProperties.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

namespace Opm::Properties {

namespace TTag {
struct TestTypeTag {
    using InheritsFrom = std::tuple<FlowModelParameters, FlowBaseProblem, BlackOilModel, EclTimeSteppingParameters>;
};
}

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::TestTypeTag> {
    using type = FlowProblem<TypeTag>;
};

// Enable experimental features for ebos: ebos is the research simulator of the OPM
// project. If you're looking for a more stable "production quality" simulator, consider
// using `flow`
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = true;
};

// use flow's well model for now
template<class TypeTag>
struct WellModel<TypeTag, TTag::TestTypeTag> {
    using type = BlackoilWellModel<TypeTag>;
};

// currently, ebos uses the non-multisegment well model by default to avoid
// regressions. the --use-multisegment-well=true|false command line parameter is still
// available in ebos, but hidden from view.
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = false;
};

// set some properties that are only required by the well model
template<class TypeTag>
struct MatrixAddWellContributions<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableTerminalOutput<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = false;
};

// flow's well model only works with surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = true;
};

// the values for the residual are for the whole cell instead of for a cubic meter of the cell
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::TestTypeTag> {
    static constexpr bool value = false;
};

// by default use flow's aquifer model for now
template<class TypeTag>
struct AquiferModel<TypeTag, TTag::TestTypeTag> {
    using type = BlackoilAquiferModel<TypeTag>;
};

// use flow's linear solver backend for now
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::TestTypeTag> {
    using type = TTag::FlowIstlSolver;
};

template<>
struct LinearSolverBackend<TTag::TestTypeTag, TTag::FlowIstlSolverParams> {
    using type = ISTLSolver<TTag::TestTypeTag>;
};

// set the maximum number of Newton iterations to 8 so that we fail quickly (albeit
// relatively often)
template<class TypeTag>
struct NewtonMaxIterations<TypeTag, TTag::TestTypeTag> {
    static constexpr int value = 8;
};

} // namespace Opm::Properties

namespace Opm::Parameters {

// if openMP is available, set the default the number of threads per process for the main
// simulation to 2 (instead of grabbing everything that is available).
#if _OPENMP
template<class TypeTag>
struct ThreadsPerProcess<TypeTag, Properties::TTag::TestTypeTag>
{ static constexpr int value = 2; };
#endif

// By default, ebos accepts the result of the time integration unconditionally if the
// smallest time step size is reached.
template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, Properties::TTag::TestTypeTag>
{ static constexpr bool value = true; };

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, Properties::TTag::TestTypeTag>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1e-1;
};

} // namespace Opm::Parameters

#endif // OPM_TEST_TYPETAG_HPP
