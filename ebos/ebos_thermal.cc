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

BEGIN_PROPERTIES

NEW_TYPE_TAG(EbosThermalTypeTag, INHERITS_FROM(EbosTypeTag));

// enable the energy extension of the black oil model
SET_BOOL_PROP(EbosThermalTypeTag, EnableEnergy, true);

END_PROPERTIES

namespace Ewoms {

void ebosThermalSetDeck(Opm::Deck* deck,
                         Opm::ParseContext* parseContext,
                         Opm::ErrorGuard* errorGuard,
                         double externalSetupTime)
{
    typedef TTAG(EbosThermalTypeTag) ProblemTypeTag;
    typedef GET_PROP_TYPE(ProblemTypeTag, Vanguard) Vanguard;

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(parseContext);
    Vanguard::setExternalErrorGuard(errorGuard);
    Vanguard::setExternalDeck(deck);
}

int ebosThermalMain(int argc, char **argv)
{
    typedef TTAG(EbosThermalTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}

}
