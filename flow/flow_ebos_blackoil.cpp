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

#include <flow/flow_ebos_blackoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Opm {

void flowEbosBlackoilSetDeck(double setupTime, std::shared_ptr<Deck> deck,
                             std::shared_ptr<EclipseState> eclState,
                             std::shared_ptr<Schedule> schedule,
                             std::unique_ptr<UDQState> udqState,
                             std::unique_ptr<Action::State> actionState,
                             std::unique_ptr<WellTestState> wtestState,
                             std::shared_ptr<SummaryConfig> summaryConfig)
{
    using TypeTag = Properties::TTag::EclFlowProblem;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

    Vanguard::setExternalSetupTime(setupTime);
    Vanguard::setExternalDeck(std::move(deck));
    Vanguard::setExternalEclState(std::move(eclState));
    Vanguard::setExternalSchedule(std::move(schedule));
    Vanguard::setExternalUDQState(std::move(udqState));
    Vanguard::setExternalActionState(std::move(actionState));
    Vanguard::setExternalWTestState(std::move(wtestState));
    Vanguard::setExternalSummaryConfig(std::move(summaryConfig));
}

std::unique_ptr<FlowMainEbos<Properties::TTag::EclFlowProblem>>
flowEbosBlackoilMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    return std::make_unique<FlowMainEbos<Properties::TTag::EclFlowProblem>>(
        argc, argv, outputCout, outputFiles);
}

// ----------------- Main program -----------------
int flowEbosBlackoilMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    auto mainfunc = flowEbosBlackoilMainInit(argc, argv, outputCout, outputFiles);
    return mainfunc->execute();
}

}
