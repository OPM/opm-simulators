/*
  Copyright 2024, SINTEF Digital

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

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>

#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include "flowexp_comp.hpp"

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpCompProblem<0>;
    Opm::registerEclTimeSteppingParameters<double>();

    // At the moment, this is probably as optimal as can be.
    // We only read the RUNSPEC of the Deck file to get the numComp,
    // and for this we need to first read the CLI arguments.
    Opm::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), true);

    auto inputFilename
        = Opm::FlowGenericVanguard::canonicalDeckPath(Opm::Parameters::Get<Opm::Parameters::EclDeckFileName>());

    // Only read the RUNSPEC section of the deck
    const auto deck = Opm::Parser{}.parseFile(inputFilename,
                                              Opm::ParseContext{},
                                              std::vector { Opm::Ecl::SectionType::RUNSPEC });
    const auto runspec = Opm::Runspec(deck);
    const auto numComps = runspec.numComps();

    OPM_ERROR_IF(numComps < 2 || numComps < 7,
                 "Deck has {} components, not supported. We support a maximum of 7 components, "
                 "and a minimum of 2.");

    switch (numComps) {
    case 2: return Opm::dispatchFlowExpComp<2>(argc, argv);
    case 3: return Opm::dispatchFlowExpComp<3>(argc, argv);
    case 4: return Opm::dispatchFlowExpComp<4>(argc, argv);
    case 5: return Opm::dispatchFlowExpComp<5>(argc, argv);
    case 6: return Opm::dispatchFlowExpComp<6>(argc, argv);
    case 7: return Opm::dispatchFlowExpComp<7>(argc, argv);
    }
}

