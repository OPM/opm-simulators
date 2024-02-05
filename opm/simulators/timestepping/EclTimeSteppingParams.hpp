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
#ifndef OPM_ECL_TIMESTEPPING_PARAMS_HPP
#define OPM_ECL_TIMESTEPPING_PARAMS_HPP

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

namespace Opm::Properties {

namespace TTag {
struct EclTimeSteppingParameters {};
}

template<class TypeTag, class MyTypeTag>
struct EnableTuning {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SolverGrowthFactor {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SolverMaxGrowth {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SolverMaxTimeStepInDays {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SolverMinTimeStep {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SolverRestartFactor {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct TimeStepAfterEventInDays {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct EnableTuning<TypeTag, TTag::EclTimeSteppingParameters> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct SolverGrowthFactor<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2.0;
};

template<class TypeTag>
struct SolverMaxGrowth<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3.0;
};

template<class TypeTag>
struct SolverMinTimeStep<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0e-12;
};

template<class TypeTag>
struct SolverMaxTimeStepInDays<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 365.0;
};

template<class TypeTag>
struct SolverRestartFactor<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.33;
};

template<class TypeTag>
struct TimeStepAfterEventInDays<TypeTag, TTag::EclTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = -1.0;
};

} // namespace Opm::Properties

namespace Opm {

template<class TypeTag>
void registerEclTimeSteppingParameters()
{
    EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTuning,
                         "Honor some aspects of the TUNING keyword.");
    EWOMS_REGISTER_PARAM(TypeTag, double, SolverGrowthFactor,
                         "The factor time steps are elongated after a successful substep");
    EWOMS_REGISTER_PARAM(TypeTag, double, SolverMaxGrowth,
                         "The maximum factor time steps are elongated after a report step");
    EWOMS_REGISTER_PARAM(TypeTag, double, SolverMaxTimeStepInDays,
                         "The maximum size of a time step in days");
    EWOMS_REGISTER_PARAM(TypeTag, double, SolverMinTimeStep,
                         "The minimum size of a time step in days for field and metric and hours for lab. If a step cannot converge without getting cut below this step size the simulator will stop");
    EWOMS_REGISTER_PARAM(TypeTag, double, SolverRestartFactor,
                         "The factor time steps are elongated after restarts");
    EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepAfterEventInDays,
                         "Time step size of the first time step after an event occurs during the simulation in days");
}

} // namespace Opm

#endif // OPM_ECL_TIME_STEPPING_PARAMS_HPP
