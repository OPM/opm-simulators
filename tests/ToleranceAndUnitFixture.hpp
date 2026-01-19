/*
  Copyright 2023 Equinor ASA

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

#pragma once

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <boost/test/unit_test.hpp>

//! \brief Shared test fixture with unit conversion helpers and tolerance checks
//!
//! This fixture provides standardized tolerance levels and unit conversion utilities
//! for OPM reservoir simulation tests. It's designed to be used as a base class
//! or included directly in test suites that need consistent numerical comparisons.
struct ToleranceAndUnitFixture
{
    //! Relative tolerance for exact mathematical relationships
    static constexpr double tight_tol = 1e-12;

    //! Relative tolerance for rate values with minor numerical differences
    static constexpr double rate_tol = 1e-8;

    //! Relative tolerance for algorithm results (e.g., scale factors)
    static constexpr double algo_tol = 1e-3;

    //! Convert SI rate (m³/s) to metric rate (SM3/day)
    static double metric_rate(double si_rate)
    {
        using namespace Opm::unit;
        return convert::to(si_rate, cubic(meter) / day);
    }

    //! Convert metric rate (SM3/day) to SI rate (m³/s)
    static double si_rate(double metric_rate)
    {
        using namespace Opm::unit;
        return convert::from(metric_rate, cubic(meter) / day);
    }

    //! Check that two values are close within tight tolerance (for exact relationships)
    static void checkClose(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, tight_tol);
    }

    //! Check that two rate values are close
    static void checkRate(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, rate_tol);
    }

    //! Check that two algorithm results are close (looser tolerance)
    static void checkAlgo(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, algo_tol);
    }
};
