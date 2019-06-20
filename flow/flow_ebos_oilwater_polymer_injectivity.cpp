/*
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

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <ewoms/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclFlowOilWaterPolymerInjectivityProblem, INHERITS_FROM(EclFlowProblem));
SET_BOOL_PROP(EclFlowOilWaterPolymerInjectivityProblem, EnablePolymer, true);
SET_BOOL_PROP(EclFlowOilWaterPolymerInjectivityProblem, EnablePolymerMW, true);
//! The indices required by the model
// For this case, there will be two primary variables introduced for the polymer
// polymer concentration and polymer molecular weight
SET_PROP(EclFlowOilWaterPolymerInjectivityProblem, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef TTAG(EclFlowProblem) BaseTypeTag;
    typedef typename GET_PROP_TYPE(BaseTypeTag, FluidSystem) FluidSystem;

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<0,
                                           2,
                                           0,
                                           /*PVOffset=*/0,
                                           /*disabledCompIdx=*/FluidSystem::gasCompIdx> type;
};
}}

namespace Opm {
/* void flowEbosOilWaterPolymerInjectivitySetDeck(Deck& deck, EclipseState& eclState)
{
    typedef TTAG(EclFlowOilWaterPolymerInjectivityProblem) TypeTag;
    typedef GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;

    Vanguard::setExternalDeck(&deck, &eclState);
} */

// ----------------- Main program -----------------
int flowEbosOilWaterPolymerInjectivityMain(int argc, char** argv)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Opm::resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    Opm::FlowMainEbos<TTAG(EclFlowOilWaterPolymerInjectivityProblem)> mainfunc;

    return mainfunc.execute(argc, argv);
}

}
