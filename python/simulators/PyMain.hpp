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

#include <opm/simulators/flow/FlowMain.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

#include <flow/flow_blackoil.hpp>

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

namespace Opm {

// ----------------- Python Main class -----------------
// Adds a python-only initialization method
class PyMain : public Main
{
public:
    using FlowMainType = FlowMain<Properties::TTag::FlowProblemTPFA>;

    using Main::Main;

    void setArguments(const std::vector<std::string>& args)
    {
        if (args.empty()) {
            return;
        }

        // We have the two arguments previously setup (binary name and input
        // case name) by the main class plus whichever args are in the
        // parameter that was passed from the python side.
        this->argc_ = 2 + args.size();

        // Setup our vector of char*'s
        argv_python_.resize(2 + args.size());
        argv_python_[0] = argv_[0];
        argv_python_[1] = argv_[1];
        for (std::size_t i = 0; i < args.size(); ++i) {
            argv_python_[i+2] = const_cast<char*>(args[i].c_str());
        }

        // Finally set the main class' argv pointer to the combined
        // parameter list.
        this->argv_ = argv_python_.data();
    }

    // To be called from the Python interface code.  Only do the
    // initialization and then return a pointer to the FlowMain object that
    // can later be accessed directly from the Python interface to
    // e.g. advance the simulator one report step
    std::unique_ptr<FlowMainType> initFlowBlackoil(int& exitCode)
    {
        exitCode = EXIT_SUCCESS;

        if (this->initialize_<Properties::TTag::FlowEarlyBird>(exitCode, true)) {
            // TODO: check that this deck really represents a blackoil
            // case. E.g. check that number of phases == 3
            this->setupVanguard();

            return flowBlackoilTpfaMainInit
                (argc_, argv_, outputCout_, outputFiles_);
        }

        // NOTE: exitCode was set by initialize_() above;
        return {}; // nullptr
    }

private:
    std::vector<char*> argv_python_{};
};

} // namespace Opm

#endif // OPM_PYMAIN_HEADER_INCLUDED
