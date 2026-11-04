/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE InjectedPhaseTest
#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/InjectedPhase.hpp>

#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

namespace {

using Indices = Opm::BlackOilDefaultFluidSystemIndices;

Opm::PhaseUsageInfo<Indices> phaseUsage(const bool oil, const bool gas, const bool water)
{
    Opm::PhaseUsageInfo<Indices> pu;
    pu.initFromPhases(Opm::Phases{oil, gas, water});
    return pu;
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(CanonicalPhaseOfInjectorType)
{
    using Opm::InjectorType;
    BOOST_CHECK_EQUAL(*Opm::injectedCanonicalPhaseIdx<Indices>(InjectorType::WATER),
                      Indices::waterPhaseIdx);
    BOOST_CHECK_EQUAL(*Opm::injectedCanonicalPhaseIdx<Indices>(InjectorType::GAS),
                      Indices::gasPhaseIdx);
    BOOST_CHECK_EQUAL(*Opm::injectedCanonicalPhaseIdx<Indices>(InjectorType::OIL),
                      Indices::oilPhaseIdx);
    BOOST_CHECK(!Opm::injectedCanonicalPhaseIdx<Indices>(InjectorType::MULTI).has_value());
}

BOOST_AUTO_TEST_CASE(ActivePhaseThreePhase)
{
    using Opm::InjectorType;
    const auto pu = phaseUsage(true, true, true);
    // Active order is water, oil, gas.
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::WATER), 0);
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::OIL), 1);
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::GAS), 2);
    BOOST_CHECK(!Opm::injectedActivePhaseIdx(pu, InjectorType::MULTI).has_value());
}

BOOST_AUTO_TEST_CASE(ActivePhaseGasWater)
{
    using Opm::InjectorType;
    // CO2STORE-like gas-water run: gas sits at active index 1, not at its
    // canonical index 2.
    const auto pu = phaseUsage(false, true, true);
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::WATER), 0);
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::GAS), 1);
    BOOST_CHECK(!Opm::injectedActivePhaseIdx(pu, InjectorType::OIL).has_value());
}

BOOST_AUTO_TEST_CASE(ActivePhaseOilGas)
{
    using Opm::InjectorType;
    const auto pu = phaseUsage(true, true, false);
    BOOST_CHECK(!Opm::injectedActivePhaseIdx(pu, InjectorType::WATER).has_value());
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::OIL), 0);
    BOOST_CHECK_EQUAL(*Opm::injectedActivePhaseIdx(pu, InjectorType::GAS), 1);
}
