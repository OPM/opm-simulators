/*
  Copyright 2020 Equinor AS.

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

#include <config.h>

#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ASTNode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/SICD.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/Valve.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WList.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WListManager.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <ebos/eclmpiserializer.hh>

#include <dune/common/parallel/mpihelper.hh>

namespace Opm {

void eclStateBroadcast(EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig,
                       UDQState& udqState,
                       Action::State& actionState)
{
    Opm::EclMpiSerializer ser(Dune::MPIHelper::getCollectiveCommunication());
    ser.broadcast(eclState);
    ser.broadcast(schedule);
    ser.broadcast(summaryConfig);
    ser.broadcast(udqState);
    ser.broadcast(actionState);
}

void eclScheduleBroadcast(Schedule& schedule)
{
    Opm::EclMpiSerializer ser(Dune::MPIHelper::getCollectiveCommunication());
    ser.broadcast(schedule);
}
}
