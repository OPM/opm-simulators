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

#define BOOST_TEST_MODULE TestOilPhaseConsistencyChecks

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <limits>
#include <string>
#include <vector>

// ###########################################################################

namespace Checks = Opm::Satfunc::PhaseChecks::Oil;

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Gas_Oil_System)

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sogcr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SOcr_GO<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SOGCR");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sogcr = 0.3f; // >= 0 && < 1

        check.test(endPoints);
    }

    {
        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 0.3f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SOcr_GO<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Sogcr value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sogcr = std::numeric_limits<float>::infinity();

        auto check = Checks::SOcr_GO<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Sogcr value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SOcr_GO<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sogcr = -0.01; // < 0

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
    auto check = Checks::SOcr_GO<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sogcr = 1.0; // >= 1

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
    auto check = Checks::SOcr_GO<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sogcr = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sogcr

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(So_min)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SOmin_GO<double>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{3});

    {
        auto columns = std::vector<std::string>(3);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SGU");
        BOOST_CHECK_EQUAL(columns[2], "SWL + SGU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.09;
        endPoints.Sgu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.09, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.9 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.99, 1.0e-8);
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
        endPoints.Sgu = 0.75f;

        auto check = Checks::SOmin_GO<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swl + Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.75f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swl + Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swl + Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.75f;

        auto check = Checks::SOmin_GO<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swl + Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.75f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgu value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swl + Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgu value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swl + Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Swl_TooLarge)
{
    auto check = Checks::SOmin_GO<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.15;
        endPoints.Sgu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.9 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 1.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sgu_TooLarge)
{
    auto check = Checks::SOmin_GO<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl = 0.1;
        endPoints.Sgu = 0.95;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.1 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.95, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 1.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // So_min

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Mobile_Oil_Sgmin)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::MobileOil_GO_SGmin<float>{};

    const auto expectNumExported = std::size_t{4};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExported);

    {
        auto columns = std::vector<std::string>(expectNumExported);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SGL");
        BOOST_CHECK_EQUAL(columns[2], "SOGCR");
        BOOST_CHECK_EQUAL(columns[3], "1 - SWL - SGL");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl   = 0.15f;
        endPoints.Sgl   = 0.01f;
        endPoints.Sogcr = 0.15f;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.01f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 1.0f - 0.15f - 0.01, 1.0e-5f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = 0.32f;

        auto check = Checks::MobileOil_GO_SGmin<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = 0.32f;

        auto check = Checks::MobileOil_GO_SGmin<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Swl_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.70;
        endPoints.Sgl   = 0.01;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.70, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.70 - 0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sgl_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgl   = 0.05;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.05, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sogcr_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgl   = 0.04;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.04, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.04, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Mobile_Oil_Sgmin

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Mobile_Oil_Sgcr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::MobileOil_GO_SGcr<float>{};

    const auto expectNumExported = std::size_t{4};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExported);

    {
        auto columns = std::vector<std::string>(expectNumExported);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SGCR");
        BOOST_CHECK_EQUAL(columns[2], "SOGCR");
        BOOST_CHECK_EQUAL(columns[3], "1 - SWL - SGCR");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl   = 0.15f;
        endPoints.Sgcr  = 0.01f;
        endPoints.Sogcr = 0.15f;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.01f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 1.0f - 0.15f - 0.01, 1.0e-5f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = 0.32f;

        auto check = Checks::MobileOil_GO_SGcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = 0.32f;

        auto check = Checks::MobileOil_GO_SGcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.01f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Swl_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.70;
        endPoints.Sgcr  = 0.01;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.70, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.70 - 0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sgcr_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgcr  = 0.05;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.05, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sogcr_TooLarge)
{
    auto check = Checks::MobileOil_GO_SGcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgcr  = 0.04;
        endPoints.Sogcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.04, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.04, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Mobile_Oil_Sgcr

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END() // Gas_Oil_System

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Oil_Water_System)

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Sowcr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SOcr_OW<float>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{1});

    {
        auto column = std::string{};
        check.columnNames(&column);

        BOOST_CHECK_EQUAL(column, "SOWCR");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sowcr = 0.3f; // >= 0 && < 1

        check.test(endPoints);
    }

    {
        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 0.3f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();

        auto check = Checks::SOcr_OW<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isnan(value), "Sowcr value must be NaN");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sowcr = std::numeric_limits<float>::infinity();

        auto check = Checks::SOcr_OW<float>{};
        check.test(endPoints);

        auto value = -0.1f;
        check.exportCheckValues(&value);

        BOOST_CHECK_MESSAGE(std::isinf(value), "Sowcr value must be Inf");
        BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
        BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
    }
}

BOOST_AUTO_TEST_CASE(Negative)
{
    auto check = Checks::SOcr_OW<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sowcr = -0.01; // < 0

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
    auto check = Checks::SOcr_OW<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sowcr = 1.0; // >= 1

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
    auto check = Checks::SOcr_OW<double>{};
    auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
    endPoints.Sowcr = 1.2; // >= 1

    check.test(endPoints);

    {
        auto value = 0.0;
        check.exportCheckValues(&value);

        BOOST_CHECK_CLOSE(value, 1.2, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Sowcr

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(So_min)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::SOmin_OW<double>{};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), std::size_t{3});

    {
        auto columns = std::vector<std::string>(3);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SGL");
        BOOST_CHECK_EQUAL(columns[1], "SWU");
        BOOST_CHECK_EQUAL(columns[2], "SGL + SWU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.09;
        endPoints.Swu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.09, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.9 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.99, 1.0e-8);
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
        endPoints.Swu = 0.75f;

        auto check = Checks::SOmin_OW<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgl + Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.75f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgl + Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgl + Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.75f;

        auto check = Checks::SOmin_OW<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgl + Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.75f;
        endPoints.Swu = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.75f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swu value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgl + Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(3);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swu value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgl + Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Sgl_TooLarge)
{
    auto check = Checks::SOmin_OW<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.15;
        endPoints.Swu = 0.9;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.9 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 1.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Swu_TooLarge)
{
    auto check = Checks::SOmin_OW<double>{};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl = 0.1;
        endPoints.Swu = 0.95;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(3);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.1 , 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.95, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 1.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Mobile_Oil_Swmin)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::MobileOil_OW_SWmin<float>{};

    const auto expectNumExported = std::size_t{4};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExported);

    {
        auto columns = std::vector<std::string>(expectNumExported);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SGL");
        BOOST_CHECK_EQUAL(columns[2], "SOWCR");
        BOOST_CHECK_EQUAL(columns[3], "1 - SWL - SGL");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl   = 0.15f;
        endPoints.Sgl   = 0.01f;
        endPoints.Sowcr = 0.15f;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.01f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 1.0f - 0.15f - 0.01, 1.0e-5f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = 0.32f;

        auto check = Checks::MobileOil_OW_SWmin<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Swl - Sgl value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = 0.32f;

        auto check = Checks::MobileOil_OW_SWmin<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.6f;
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Swl - Sgl value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Swl_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.70;
        endPoints.Sgl   = 0.01;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.70, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.70 - 0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sgl_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgl   = 0.05;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.05, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sowcr_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWmin<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Swl   = 0.65;
        endPoints.Sgl   = 0.04;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.04, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.04, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Mobile_Oil_Swmin

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Mobile_Oil_Swcr)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::MobileOil_OW_SWcr<float>{};

    const auto expectNumExported = std::size_t{4};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExported);

    {
        auto columns = std::vector<std::string>(expectNumExported);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SGL");
        BOOST_CHECK_EQUAL(columns[1], "SWCR");
        BOOST_CHECK_EQUAL(columns[2], "SOWCR");
        BOOST_CHECK_EQUAL(columns[3], "1 - SWCR - SGL");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl   = 0.15f;
        endPoints.Swcr  = 0.01f;
        endPoints.Sowcr = 0.15f;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.01f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 1.0f - 0.15f - 0.01, 1.0e-5f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = 0.32f;

        auto check = Checks::MobileOil_OW_SWcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Sgl - Swcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.6f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Sgl - Swcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.6f;
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Sgl - Swcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Sgl - Swcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgu value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1 - Sgl - Swcr value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        const auto expectNumExported = std::size_t{4};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = 0.32f;

        auto check = Checks::MobileOil_OW_SWcr<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Sgl - Swcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.6f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Sgl - Swcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.6f;
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE  (values[0], 0.6f, 1.0e-6f);
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[3], 1.0f - 0.6f - 0.01f, 1.0e-6);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.32f;
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE  (values[2], 0.32f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Sgl - Swcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.01f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE  (values[1], 0.01f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Sgl - Swcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExported);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1 - Sgl - Swcr value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Sgl_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl   = 0.70;
        endPoints.Swcr  = 0.01;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.70, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.70 - 0.01, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Swcr_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl   = 0.65;
        endPoints.Swcr  = 0.05;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.05, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.05, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sowcr_TooLarge)
{
    auto check = Checks::MobileOil_OW_SWcr<double>{};

    const auto expectNumExported = std::size_t{4};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        endPoints.Sgl   = 0.65;
        endPoints.Swcr  = 0.04;
        endPoints.Sowcr = 0.32;

        check.test(endPoints);
    }

    {
        auto values = std::vector<double>(expectNumExported);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(values[1], 0.04, 1.0e-8);
        BOOST_CHECK_CLOSE(values[2], 0.32, 1.0e-8);
        BOOST_CHECK_CLOSE(values[3], 1.0 - 0.65 - 0.04, 1.0e-8);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Mobile_Oil_Swcr

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END() // Oil_Water_System
