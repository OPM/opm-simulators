/*
  Copyright 2024 Equinor AS

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

#define BOOST_TEST_MODULE TestWaterPhaseConsistencyChecks

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/WaterPhaseConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <limits>
#include <string>
#include <vector>

// ###########################################################################

namespace Checks = Opm::Satfunc::PhaseChecks::Water;

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Sw_min)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SWmin<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SWL");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.125f; // >= 0 && < 1

        check.test(endPoints);
    }

    {
        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 0.125f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SWmin<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Swl value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();

        auto check = Checks::SWmin<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Swl value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SWmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swl = -0.01; // < 0

    check.test(endPoints);

    {
        auto value = 1.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, -0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Is_One)
{
    auto check = Checks::SWmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swl = 1.0; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.0, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Exceeds_One)
{
    auto check = Checks::SWmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swl = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sw_min

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sw_max)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SWmax<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SWU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swu = 0.125f; // >= 0 && < 1

        check.test(endPoints);
    }

    {
        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 0.125f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(All_Good_Swu_Is_One)
{
    auto check = Checks::SWmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swu = 1.0; // <= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.0, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SWmax<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Swu value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swu = std::numeric_limits<float>::infinity();

        auto check = Checks::SWmax<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Swu value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SWmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swu = -0.01; // < 0

    check.test(endPoints);

    {
        auto value = 1.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, -0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Zero)
{
    auto check = Checks::SWmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swu = 0.0; // <= 0

    check.test(endPoints);

    {
        auto value = 1.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 0.0, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Exceeds_One)
{
    auto check = Checks::SWmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Swu = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sw_max

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sw_cr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SWcr<double>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{3});

    {
        auto columns = std::vector<std::string>(3);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SWCR");
        BOOST_CHECK_EQUAL(columns[2], "SWU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.09;
        endPoints.Swcr = 0.12;
        endPoints.Swu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.09, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.12, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.90, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(All_Good_Swcr_Same_As_Swl)
{
    auto check = Checks::SWcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.01;
        endPoints.Swcr = 0.01;
        endPoints.Swu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.90, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.125f;
        endPoints.Swu = 0.75f;

        auto check = Checks::SWcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.75;
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = 0.125f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.75f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.125f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.125f;
        endPoints.Swu = 0.75f;

        auto check = Checks::SWcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.75;
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = 0.125f;
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.75f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.125f;
        endPoints.Swu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.01f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Swcr_TooSmall)
{
    auto check = Checks::SWcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.15;
        endPoints.Swcr = 0.125;
        endPoints.Swu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.125, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.9, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Swcr_TooLarge)
{
    auto check = Checks::SWcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.15;
        endPoints.Swcr = 0.65;
        endPoints.Swu = 0.6;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.6, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Swcr_Same_As_Swu)
{
    auto check = Checks::SWcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.15;
        endPoints.Swcr = 0.65;
        endPoints.Swu = 0.65;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.65, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Swu_TooSmall)
{
    auto check = Checks::SWcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.15;
        endPoints.Swcr = 0.15;
        endPoints.Swu = 0.10;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.10, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // So_min
