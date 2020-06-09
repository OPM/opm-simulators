/*
  Copyright 2013, 2014, 2015, 2019 SINTEF Digital, Mathematics and Cybernetics.
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
#include <opm/simulators/linalg/ISTLSolverEbosFlexible.hpp>
#include <ebos/eclfluxmoduleseq.hh>
BEGIN_PROPERTIES
NEW_TYPE_TAG(EclFlowProblemSimple, INHERITS_FROM(EclFlowProblem));
NEW_PROP_TAG(FluidState);
//SET_TYPE_PROP(EclBaseProblem, Problem, Opm::EclProblem<TypeTag>);
SET_PROP(EclFlowProblemSimple, FluidState)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
      enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
      enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
      enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
      enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
      static const bool compositionSwitchEnabled = Indices::gasEnabled;

    public:
//typedef Opm::BlackOilFluidSystemSimple<Scalar> type;
       typedef Opm::BlackOilFluidState<Evaluation, FluidSystem, enableTemperature, enableEnergy, compositionSwitchEnabled,  Indices::numPhases > type;
};

SET_BOOL_PROP(EclFlowProblemSimple, MatrixAddWellContributions, true);
SET_INT_PROP(EclFlowProblemSimple, LinearSolverVerbosity,0);
SET_SCALAR_PROP(EclFlowProblemSimple, LinearSolverReduction, 1e-2);
SET_INT_PROP(EclFlowProblemSimple, LinearSolverMaxIter, 100);
SET_BOOL_PROP(EclFlowProblemSimple, UseAmg, true);//probably not used
SET_BOOL_PROP(EclFlowProblemSimple, UseCpr, true);
SET_INT_PROP(EclFlowProblemSimple, CprMaxEllIter, 1);
SET_INT_PROP(EclFlowProblemSimple, CprEllSolvetype, 3);
SET_INT_PROP(EclFlowProblemSimple, CprReuseSetup, 3);
SET_INT_PROP(EclFlowProblemSimple, CprSolverVerbose, 0);
SET_STRING_PROP(EclFlowProblemSimple, LinearSolverConfiguration, "ilu0");
//SET_STRING_PROP(EclFlowProblemSimple, LinearSolverConfiguration, "file");
SET_STRING_PROP(EclFlowProblemSimple, SystemStrategy, "quasiimpes");
END_PROPERTIES

namespace Opm {
    //template <class TypeTag>
    //class EclTransExtensiveQuantitiesSeq : public EclTransExtensiveQuantitise<TypeTag>
    //{};

    // template <class TypeTag>
    // BlackOilIntensiveQuantitiesSeq    :  public BlackOilIntensiveQuantities<TypeTag>
    // {};

    // template <class TypeTag>
    // BlackOilPrimaryVariablesSeq    :  public BlackOilPrimaryVariables<TypeTag>
    // {};
    
template <class TypeTag>
struct EclTransFluxModuleSeq
{
    typedef EclTransIntensiveQuantities<TypeTag> FluxIntensiveQuantities;
    typedef EclTransExtensiveQuantitiesSeq<TypeTag> FluxExtensiveQuantities;
    typedef EclTransBaseProblem<TypeTag> FluxBaseProblem;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};
}


namespace Opm {
  namespace Properties {

    SET_PROP(EclFlowProblemSimple, FluidSystem)
    {
    private:
      //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

    public:
        typedef Opm::BlackOilFluidSystem<Scalar> type;
    };
    //NEW_TYPE_TAG(EclFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem));
     // template <class TypeTag>
    // BlackOilIntensiveQuantitiesSeq    :  public BlackOilIntensiveQuantities<TypeTag>
    // {};
    SET_TYPE_PROP(EclFlowProblemSimple, PrimaryVariables, Opm::BlackOilPrimaryVariables<TypeTag>);
    SET_TYPE_PROP(EclFlowProblemSimple, IntensiveQuantities, Opm::BlackOilIntensiveQuantities<TypeTag>);
    SET_TYPE_PROP(EclFlowProblemSimple, FluxModule, Opm::EclTransFluxModuleSeq<TypeTag>);
    //SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::ISTLSolverEbos<TypeTag>);
    //SET_TAG_PROP(EclFlowProblemSimple, LinearSolverSplice, ParallelBiCGStabLinearSolver);
    //SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::Linear::ParallelBiCGStabSolverBackend<TypeTag>);//not work
    //SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::Linear::SuperLUBackend<TypeTag>)//not work
    //SET_TAG_PROP(EclFlowProblem, FluidState, Opm::BlackOilFluidState);
    SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::ISTLSolverEbosFlexible<TypeTag>);
    SET_BOOL_PROP(EclFlowProblemSimple, EnableStorageCache, true);
    SET_BOOL_PROP(EclFlowProblemSimple, EnableIntensiveQuantityCache, true);

    //SET_INT_PROP(EclFlowProblemSimple, NumWellAdjoint, 1);
    //SET_BOOL_PROP(EclFlowProblem, EnableStorageCache, true);
    //SET_BOOL_PROP(EclFlowProblem, EnableIntensiveQuantityCache, true);
  }
}

int main(int argc, char** argv)
{
    using TypeTag = TTAG(EclFlowProblemSimple);
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
