/*
  Copyright 2017 IRIS AS

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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <opm/material/densead/Evaluation.hpp>
#include <ewoms/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/autodiff/DuneMatrix.hpp>
#include <dune/grid/CpGrid.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/autodiff/FlowMainEbos.hpp>

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclFlowTwoPhaseProblem, INHERITS_FROM(EclFlowProblem));
//! The indices required by the model
SET_TYPE_PROP(EclFlowTwoPhaseProblem, Indices,
              Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent)?1:0, GET_PROP_VALUE(TypeTag, EnablePolymer)?1:0, /*PVOffset=*/0, /*disabledCompIdx=*/2>);
}}

// ----------------- Main program -----------------
int main(int argc, char** argv)
{
    Opm::FlowMainEbos<TTAG(EclFlowTwoPhaseProblem)> mainfunc;
    return mainfunc.execute(argc, argv);
}
