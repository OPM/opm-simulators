/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015, 2017 IRIS AS

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
#include "config.h"
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoilonephaseindices.hh>


BEGIN_PROPERTIES
NEW_TYPE_TAG(EclFlowProblemSimple, INHERITS_FROM(EclFlowProblem));
SET_BOOL_PROP(EclFlowProblemSimple, EnableEnergy, true);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(FluidSystem);
//! The indices required by the model
SET_PROP(EclFlowProblemSimple, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef TTAG(EclFlowProblem) BaseTypeTag;
    typedef typename GET_PROP_TYPE(BaseTypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::BlackOilOnePhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*enebledCompIdx=*/FluidSystem::waterCompIdx>
        type;
};
SET_PROP(EclFlowProblemSimple, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableBrine = GET_PROP_VALUE(TypeTag, EnableBrine) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    static const bool compositionSwitchEnabled = Indices::gasEnabled;

public:
    // typedef Opm::BlackOilFluidSystemSimple<Scalar> type;
    typedef Opm::BlackOilFluidState<Evaluation,
                                    FluidSystem,
                                    enableTemperature,
                                    enableEnergy,
                                    compositionSwitchEnabled,
                                    enableBrine,
                                    Indices::numPhases>
        type;
};

// SET_PROP(EclFlowProblemSimple, FluidSystem)
// {
// private:
//   //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
//   typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

// public:
//   typedef Opm::BlackOilFluidSystem<Scalar,Indices> type;
// };
END_PROPERTIES

int main(int argc, char** argv)
{
    using TypeTag = TTAG(EclFlowProblemSimple);
    auto mainObject = Opm::Main<TypeTag>(argc, argv);
    return mainObject.run();
}
