/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.
  Copyright 2023 Inria

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
#ifndef OPM_PYMAIN_HEADER_INCLUDED
#define OPM_PYMAIN_HEADER_INCLUDED

#include <opm/simulators/flow/Main.hpp>

namespace Opm {

// ----------------- Python Main class -----------------
// Adds a python-only initialization method
class PyMain : public Main
{
public:
    using Main::Main;

    using FlowMainType = FlowMain<Properties::TTag::FlowProblemTPFA>;
    // To be called from the Python interface code. Only do the
    // initialization and then return a pointer to the FlowMain
    // object that can later be accessed directly from the Python interface
    // to e.g. advance the simulator one report step
    std::unique_ptr<FlowMainType> initFlowBlackoil(int& exitCode)
    {
        exitCode = EXIT_SUCCESS;
        if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
            // TODO: check that this deck really represents a blackoil
            // case. E.g. check that number of phases == 3
            this->setupVanguard();
            return flowBlackoilTpfaMainInit(
                argc_, argv_, outputCout_, outputFiles_);
        } else {
            //NOTE: exitCode was set by initialize_() above;
            return std::unique_ptr<FlowMainType>(); // nullptr
        }
    }
};

} // namespace Opm

#endif // OPM_PYMAIN_HEADER_INCLUDED
