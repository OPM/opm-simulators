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

namespace Opm {

bool ebosBlackOilDeckFileNameIsSet(int argc, char** argv)
{
    using ProblemTypeTag = Properties::TTag::EbosTypeTag;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    Parameters::reset<ProblemTypeTag>();
    setupParameters_<ProblemTypeTag>(argc,
                                     const_cast<const char**>(argv),
                                     /*doRegistration=*/true,
                                     /*allowUnused=*/true,
                                     /*handleHelp=*/false);
    bool result = EWOMS_PARAM_IS_SET(ProblemTypeTag, std::string, EclDeckFileName);
    Parameters::reset<ProblemTypeTag>();

    return result;
}

std::string ebosBlackOilGetDeckFileName(int argc, char** argv)
{
    using ProblemTypeTag = Properties::TTag::EbosTypeTag;
    using Vanguard = GetPropType<ProblemTypeTag, Properties::Vanguard>;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    Parameters::reset<ProblemTypeTag>();
    setupParameters_<ProblemTypeTag>(argc,
                                     const_cast<const char**>(argv),
                                     /*doRegistration=*/true,
                                     /*allowUnused=*/true,
                                     /*handleHelp=*/false);
    std::string rawDeckFileName = EWOMS_GET_PARAM(ProblemTypeTag, std::string, EclDeckFileName);
    std::string result = Vanguard::canonicalDeckPath(rawDeckFileName).string();
    Parameters::reset<ProblemTypeTag>();

    return result;
}

std::unique_ptr<ParseContext> ebosBlackOilCreateParseContext(int argc, char** argv)
{
    using ProblemTypeTag = Properties::TTag::EbosTypeTag;
    using Vanguard = GetPropType<ProblemTypeTag, Properties::Vanguard>;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    Parameters::reset<ProblemTypeTag>();
    setupParameters_<ProblemTypeTag>(argc,
                                     const_cast<const char**>(argv),
                                     /*doRegistration=*/true,
                                     /*allowUnused=*/true,
                                     /*handleHelp=*/false);
    std::unique_ptr<ParseContext> result = Vanguard::createParseContext();
    Parameters::reset<ProblemTypeTag>();

    return result;
}

void ebosBlackOilSetDeck(std::unique_ptr<Deck> deck,
                         std::unique_ptr<ParseContext> parseContext,
                         std::unique_ptr<ErrorGuard> errorGuard,
                         double externalSetupTime)
{
    using ProblemTypeTag = Properties::TTag::EbosTypeTag;
    using Vanguard = GetPropType<ProblemTypeTag, Properties::Vanguard>;

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(std::move(parseContext));
    Vanguard::setExternalErrorGuard(std::move(errorGuard));
    Vanguard::setExternalDeck(std::move(deck));
}

int ebosBlackOilMain(int argc, char **argv)
{
    using ProblemTypeTag = Properties::TTag::EbosTypeTag;
    return startEbos<ProblemTypeTag>(argc, argv);
}

}
