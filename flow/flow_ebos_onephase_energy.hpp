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

#ifndef FLOW_ONEPHASE_ENERGY_HPP
#define FLOW_ONEPHASE_ENERGY_HPP

#include <memory>

namespace Opm {

class Deck;
class EclipseState;
template<class TypeTag> class FlowMainEbos;
class Schedule;
class SummaryConfig;
class UDQState;
class WellTestState;
namespace Action {
class State;
}

void flowEbosWaterOnlyEnergySetDeck(double setupTime, std::shared_ptr<Deck> deck,
                                    std::shared_ptr<EclipseState> eclState,
                                    std::shared_ptr<Schedule> schedule,
                                    std::shared_ptr<SummaryConfig> summaryConfig);
int flowEbosWaterOnlyEnergyMain(int argc, char** argv, bool outputCout, bool outputFiles);
int flowEbosWaterOnlyEnergyMainStandalone(int argc, char** argv);
}

#endif
