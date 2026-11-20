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

#include <opm/simulators/utils/UnsupportedFlowKeywords.hpp>

#include <boost/test/unit_test.hpp>

// The supported dual-porosity subset must not be blocked...
BOOST_AUTO_TEST_CASE(DualPorositySubsetNotBlocked)
{
    const auto& blocked = Opm::FlowKeywordValidation::unsupportedKeywords();
    for (const auto* kw : {"DUALPORO", "DUALPERM", "NODPPM", "SIGMA", "SIGMAV", "DPGRID"}) {
        BOOST_CHECK_MESSAGE(blocked.find(kw) == blocked.end(),
                            std::string{kw} + " must not be in the unsupported-keyword list");
    }
}

// ...while every deferred dual-continuum keyword must remain a hard error.
