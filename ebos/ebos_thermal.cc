// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief A general-purpose simulator for ECL decks using the black-oil model.
 */
#include "config.h"

#include "ebos.hh"
#include "startEbos.hh"

namespace Opm::Properties {

namespace TTag {
struct EbosThermalTypeTag {
    using InheritsFrom = std::tuple<EbosTypeTag>;
};
}

// enable the energy extension of the black oil model
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::EbosThermalTypeTag> {
    static constexpr bool value = true;
};

} // namespace Opm::Properties

namespace Opm {

void ebosThermalSetDeck(std::unique_ptr<Opm::Deck> deck,
                        std::unique_ptr<Opm::ParseContext> parseContext,
                        std::unique_ptr<Opm::ErrorGuard> errorGuard,
                        double externalSetupTime)
{
    using ProblemTypeTag = Properties::TTag::EbosThermalTypeTag;
    using Vanguard = GetPropType<ProblemTypeTag, Properties::Vanguard>;

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(std::move(parseContext));
    Vanguard::setExternalErrorGuard(std::move(errorGuard));
    Vanguard::setExternalDeck(std::move(deck));
}

int ebosThermalMain(int argc, char **argv)
{
    using ProblemTypeTag = Properties::TTag::EbosThermalTypeTag;
    return Opm::startEbos<ProblemTypeTag>(argc, argv);
}

}
