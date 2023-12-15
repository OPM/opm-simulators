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
#ifndef FLOW_EBOS_BLACKOIL_TPFA_HPP
#define FLOW_EBOS_BLACKOIL_TPFA_HPP

#include <memory>
//#include <flow/flow_ebos_energy.hpp>
namespace Opm {
    namespace Properties {
        namespace TTag {

            struct EclFlowProblem;

            struct EclFlowProblemTPFA {
            using InheritsFrom = std::tuple<EclFlowProblem>;
            };

        }
   }
}

namespace Opm {


//! \brief Main function used in flow binary.
int flowEbosBlackoilTpfaMain(int argc, char** argv, bool outputCout, bool outputFiles);

template<class TypeTag> class FlowMainEbos;

//! \brief Initialization function used in flow binary and python simulator.
std::unique_ptr<FlowMainEbos<Properties::TTag::EclFlowProblemTPFA>>
    flowEbosBlackoilTpfaMainInit(int argc, char** argv, bool outputCout, bool outputFiles);

//! \brief Main function used in flow_brine binary.
int flowEbosBlackoilTpfaMainStandalone(int argc, char** argv);

}

#endif // FLOW_EBOS_BLACKOIL_TPFA_HPP
