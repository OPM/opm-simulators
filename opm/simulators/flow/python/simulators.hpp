/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#ifndef OPM_SIMULATORS_HEADER_INCLUDED
#define OPM_SIMULATORS_HEADER_INCLUDED

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>

namespace Opm::Pybind {
class BlackOilSimulator
{
private:
    using FlowMainEbosType = Opm::FlowMainEbos<TTAG(EclFlowProblem)>;

public:
    BlackOilSimulator( const std::string &deckFilename);
    int run();
    int step_init();

private:
    const std::string deckFilename_;
    std::unique_ptr<FlowMainEbosType> mainEbos_;
    std::unique_ptr<Opm::Main> main_;
    bool hasRunInit_;
};

} // namespace Opm::Python
#endif // OPM_SIMULATORS_HEADER_INCLUDED
