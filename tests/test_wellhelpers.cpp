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

#include <config.h>

#define BOOST_TEST_MODULE WellHelpersTest
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/UDAValue.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Well/WellInjectionControls.hpp>
#include <opm/input/eclipse/Schedule/Well/WellProductionControls.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/wells/WellHelpers.hpp>

namespace {

// The SI rate a rate UDA evaluating to zero or a negative value gives, see
// UDA::eval_well_uda() and UDA::eval_group_uda().
double udaEpsilonRate(const double unit_rate)
{
    return Opm::UDAValue{}.epsilonLimit() * unit_rate;
}

const double sm3_per_day = Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
const double mscf_per_day = 1000.0 * Opm::unit::cubic(Opm::unit::feet) / Opm::unit::day;

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(ZeroRateTarget)
{
    using Opm::wellhelpers::isZeroRateTarget;

    BOOST_CHECK(isZeroRateTarget(0.0));
    BOOST_CHECK(isZeroRateTarget(-0.0));
    BOOST_CHECK(isZeroRateTarget(udaEpsilonRate(sm3_per_day)));
    BOOST_CHECK(isZeroRateTarget(udaEpsilonRate(mscf_per_day)));
    BOOST_CHECK(isZeroRateTarget(0.0f));
    BOOST_CHECK(isZeroRateTarget(static_cast<float>(udaEpsilonRate(sm3_per_day))));

    BOOST_CHECK(!isZeroRateTarget(1.0e-6 * sm3_per_day));
    BOOST_CHECK(!isZeroRateTarget(1.0 * sm3_per_day));
    BOOST_CHECK(!isZeroRateTarget(1.0e-6f));
}

BOOST_AUTO_TEST_CASE(ZeroProductionRateTarget)
{
    using Opm::WellProducerCMode;
    using Opm::wellhelpers::rateControlWithZeroProdTarget;

    Opm::WellProductionControls controls(0);
    controls.oil_rate = 100.0 * sm3_per_day;
    controls.water_rate = 0.0;
    controls.gas_rate = udaEpsilonRate(sm3_per_day);
    controls.liquid_rate = 1.0e-3 * sm3_per_day;
    controls.resv_rate = udaEpsilonRate(sm3_per_day);
    controls.prediction_mode = true;

    BOOST_CHECK(!rateControlWithZeroProdTarget(controls, WellProducerCMode::ORAT));
    BOOST_CHECK(rateControlWithZeroProdTarget(controls, WellProducerCMode::WRAT));
    BOOST_CHECK(rateControlWithZeroProdTarget(controls, WellProducerCMode::GRAT));
    BOOST_CHECK(!rateControlWithZeroProdTarget(controls, WellProducerCMode::LRAT));
    BOOST_CHECK(rateControlWithZeroProdTarget(controls, WellProducerCMode::RESV));
    BOOST_CHECK(!rateControlWithZeroProdTarget(controls, WellProducerCMode::BHP));
    BOOST_CHECK(!rateControlWithZeroProdTarget(controls, WellProducerCMode::THP));
}

BOOST_AUTO_TEST_CASE(ZeroInjectionRateTarget)
{
    using Opm::WellInjectorCMode;
    using Opm::wellhelpers::rateControlWithZeroInjTarget;

    Opm::WellInjectionControls controls(0);
    controls.surface_rate = udaEpsilonRate(sm3_per_day);
    controls.reservoir_rate = 10.0 * sm3_per_day;

    BOOST_CHECK(rateControlWithZeroInjTarget(controls, WellInjectorCMode::RATE));
    BOOST_CHECK(!rateControlWithZeroInjTarget(controls, WellInjectorCMode::RESV));
    BOOST_CHECK(!rateControlWithZeroInjTarget(controls, WellInjectorCMode::BHP));

    controls.reservoir_rate = 0.0;
    BOOST_CHECK(rateControlWithZeroInjTarget(controls, WellInjectorCMode::RESV));
}
