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

#include <opm/simulators/flow_ebos_polymer.hpp>

#include <opm/common/utility/ResetLocale.hpp>
#include <dune/grid/CpGrid.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/autodiff/FlowMainEbos.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclFlowPolymerProblem, INHERITS_FROM(EclFlowProblem));
SET_BOOL_PROP(EclFlowPolymerProblem, EnablePolymer, true);
}}

namespace Opm {
void flowEbosPolymerSetDeck(Deck &deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
{
    typedef TTAG(EclFlowPolymerProblem) TypeTag;
    typedef GET_PROP_TYPE(TypeTag, GridManager) GridManager;

    GridManager::setExternalDeck(&deck, &eclState, &schedule, &summaryConfig);
}

// ----------------- Main program -----------------
int flowEbosPolymerMain(int argc, char** argv)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Opm::resetLocale();

    // initialize MPI, finalize is done automatically on exit
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv).rank();
#endif

    Opm::FlowMainEbos<TTAG(EclFlowPolymerProblem)> mainfunc;
    return mainfunc.execute(argc, argv);
}

}
