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

#define BOOST_TEST_MODULE TestGasPhaseConsistencyChecks

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/GasPhaseConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <limits>
#include <string>
#include <vector>

// ###########################################################################

namespace Checks = Opm::Satfunc::PhaseChecks::Gas;

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Sg_min)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SGmin<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SGL");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.125f; // >= 0 && < 1

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
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SGmin<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Sgl value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::infinity();

        auto check = Checks::SGmin<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Sgl value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SGmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgl = -0.01; // < 0

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
    auto check = Checks::SGmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgl = 1.0; // >= 1

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
    auto check = Checks::SGmin<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgl = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sg_min

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sg_max)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SGmax<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SGU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgu = 0.125f; // >= 0 && < 1

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
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SGmax<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Sgu value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        auto check = Checks::SGmax<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Sgu value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SGmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgu = -0.01; // < 0

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
    auto check = Checks::SGmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgu = 1.0; // >= 1

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
    auto check = Checks::SGmax<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sgu = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sg_max

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sg_cr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SGcr<double>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{3});

    {
        auto columns = std::vector<std::string>(3);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SGL");
        BOOST_CHECK_EQUAL(columns[1], "SGCR");
        BOOST_CHECK_EQUAL(columns[2], "SGU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.09;
        endPoints.Sgcr = 0.12;
        endPoints.Sgu = 0.9;

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

BOOST_AUTO_TEST_CASE(All_Good_Sgcr_Same_As_Sgl)
{
    auto check = Checks::SGcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.01;
        endPoints.Sgcr = 0.01;
        endPoints.Sgu = 0.9;

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
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = 0.75f;

        auto check = Checks::SGcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.75;
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.75f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = 0.75f;

        auto check = Checks::SGcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.75;
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.75f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.75f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.125f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.125f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.01f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();;

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Sgcr_TooSmall)
{
    auto check = Checks::SGcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.15;
        endPoints.Sgcr = 0.125;
        endPoints.Sgu = 0.9;

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

BOOST_AUTO_TEST_CASE(Sgcr_TooLarge)
{
    auto check = Checks::SGcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.15;
        endPoints.Sgcr = 0.65;
        endPoints.Sgu = 0.6;

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

BOOST_AUTO_TEST_CASE(Sgcr_Same_As_Sgu)
{
    auto check = Checks::SGcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.15;
        endPoints.Sgcr = 0.65;
        endPoints.Sgu = 0.65;

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

BOOST_AUTO_TEST_CASE(Sgu_TooSmall)
{
    auto check = Checks::SGcr<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.15;
        endPoints.Sgcr = 0.15;
        endPoints.Sgu = 0.10;

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
