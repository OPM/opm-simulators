/*
  Copyright 2026 TNO

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

#define BOOST_TEST_MODULE TestDualPorosityKeywords

#include <opm/simulators/flow/KeywordValidation.hpp>
#include <opm/simulators/utils/UnsupportedFlowKeywords.hpp>

#include <opm/common/utility/OpmInputError.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <boost/test/unit_test.hpp>

#include <initializer_list>
#include <string>

// The dual-continuum keywords as flow validates a deck: the supported ones pass, and each one that is
// not implemented yet stops the run, unless --parsing-strictness=low, which reports it as a warning.
namespace {

enum class Section { Runspec, Grid, Props };

Opm::Deck dualContinuumDeck(const std::string& model, Section section, const std::string& keyword)
{
    const auto in = [&](Section s) { return s == section ? keyword + "\n" : std::string{}; };
    return Opm::Parser{}.parseString("RUNSPEC\nDIMENS\n 2 1 2 /\nTABDIMS\n/\n" + model + "\n" + in(Section::Runspec)
                                     + "GRID\n" + in(Section::Grid)
                                     + "PROPS\n" + in(Section::Props));
}

// Validates the deck as flow does with its list of unsupported keywords; lowStrictness is
// --parsing-strictness=low. A critical problem at normal strictness throws.
bool hasErrors(const Opm::Deck& deck, bool lowStrictness)
{
    const Opm::KeywordValidation::KeywordValidator validator(Opm::FlowKeywordValidation::unsupportedKeywords(),
                                                             {{}, {}, {}}, {{}, {}, {}}, {});
    Opm::ErrorGuard errors;
    validator.validateDeck(deck, Opm::ParseContext{}, lowStrictness, errors);
    const bool found = static_cast<bool>(errors);
    errors.clear();
    return found;
}

bool throwsInputError(const Opm::Deck& deck)
{
    try {
        hasErrors(deck, /*lowStrictness=*/false);
    }
    catch (const Opm::OpmInputError&) {
        return true;
    }
    return false;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(SupportedDualPorosityKeywordsPassValidation)
{
    const std::string grid = "DPGRID\nSIGMA\n 0.1 /\nSIGMAV\n 4*0.1 /";
    for (const auto* model : {"DUALPORO\nNODPPM", "DUALPERM\nNODPPM"}) {
        BOOST_CHECK_MESSAGE(!throwsInputError(dualContinuumDeck(model, Section::Grid, grid)),
                            std::string{model} + " with DPGRID, SIGMA and SIGMAV must pass the validation");
    }
}

BOOST_AUTO_TEST_CASE(DeferredDualPorosityKeywordsStopTheRun)
{
    struct Deferred { const char* name; Section section; const char* keyword; };
    for (const Deferred& kw : std::initializer_list<Deferred>{
        {"DIFFMMF", Section::Grid, "DIFFMMF\n 4*1.0 /"},
        {"DPKRMOD", Section::Props, "DPKRMOD\n /"},
        {"DPNUM", Section::Grid, "DPNUM\n 4*1 /"},
        {"DZMTRX", Section::Grid, "DZMTRX\n 1.0 /"},
        {"DZMTRXV", Section::Grid, "DZMTRXV\n 4*1.0 /"},
        {"GRAVDR", Section::Runspec, "GRAVDR"},
        {"GRAVDRB", Section::Runspec, "GRAVDRB"},
        {"GRAVDRM", Section::Runspec, "GRAVDRM\n /"},
        {"LTOSIGMA", Section::Grid, "LTOSIGMA\n /"},
        {"LX", Section::Grid, "LX\n 4*1.0 /"},
        {"LY", Section::Grid, "LY\n 4*1.0 /"},
        {"LZ", Section::Grid, "LZ\n 4*1.0 /"},
        {"NMATOPTS", Section::Grid, "NMATOPTS\n /"},
        {"NMATRIX", Section::Runspec, "NMATRIX\n 2 /"},
        {"SIGMAGDV", Section::Grid, "SIGMAGDV\n 4*0.1 /"},
        {"SIGMATH", Section::Grid, "SIGMATH\n 4*0.1 /"},
    }) {
        const auto deck = dualContinuumDeck("DUALPORO", kw.section, kw.keyword);
        BOOST_CHECK_MESSAGE(throwsInputError(deck), std::string{kw.name} + " must stop the run");
        BOOST_CHECK_MESSAGE(!hasErrors(deck, /*lowStrictness=*/true),
                            std::string{kw.name} + " must only warn with --parsing-strictness=low");
    }
}
