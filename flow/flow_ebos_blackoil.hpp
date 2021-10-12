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
#ifndef FLOW_EBOS_BLACKOIL_HPP
#define FLOW_EBOS_BLACKOIL_HPP

#include <memory>

namespace Opm {

class Deck;
class EclipseState;
template<class TypeTag> class FlowMainEbos;
class Schedule;
class SummaryConfig;
class UDQState;
namespace Action {
class State;
}
namespace Properties { namespace TTag { struct EclFlowProblem; } }

void flowEbosBlackoilSetDeck(double setupTime, std::shared_ptr<Deck> deck,
                             std::shared_ptr<EclipseState> eclState,
                             std::shared_ptr<Schedule> schedule,
                             std::unique_ptr<UDQState> udqState,
                             std::unique_ptr<Action::State> actionState,
                             std::shared_ptr<SummaryConfig> summaryConfig);

int flowEbosBlackoilMain(int argc, char** argv, bool outputCout, bool outputFiles);

std::unique_ptr<FlowMainEbos<Properties::TTag::EclFlowProblem>>
    flowEbosBlackoilMainInit(int argc, char** argv, bool outputCout, bool outputFiles);
}

#endif // FLOW_EBOS_BLACKOIL_HPP
