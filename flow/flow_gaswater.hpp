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
#ifndef FLOW_GASWATER_HPP
#define FLOW_GASWATER_HPP

#include <opm/simulators/flow/TTagFlowProblemGasWater.hpp>

#include <memory>

namespace Opm {

//! \brief Main function used in flow binary.
int flowGasWaterMain(int argc, char** argv, bool outputCout, bool outputFiles);

//! \brief Main function used in flow_gaswater binary.
int flowGasWaterMainStandalone(int argc, char** argv);

template<class TypeTag> class FlowMain;

std::unique_ptr<FlowMain<Properties::TTag::FlowGasWaterProblem>>
flowGasWaterMainInit(int argc, char** argv, bool outputCout, bool outputFiles);
}

#endif // FLOW_GASWATER_HPP
