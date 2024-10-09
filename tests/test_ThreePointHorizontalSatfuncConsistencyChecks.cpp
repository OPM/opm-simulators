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

#define BOOST_TEST_MODULE TestThreePointHorizontalConsistencyChecks

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/ThreePointHorizontalConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <limits>
#include <string>
#include <vector>

// ###########################################################################

namespace Checks = Opm::Satfunc::PhaseChecks::ThreePointHorizontal;

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Displacing_Oil_in_Gas_Oil)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExportedCheckValues);

    {
        auto columns = std::vector<std::string>(expectNumExportedCheckValues);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SWL");
        BOOST_CHECK_EQUAL(columns[1], "SOGCR");
        BOOST_CHECK_EQUAL(columns[2], "SGCR");
        BOOST_CHECK_EQUAL(columns[3], "1-SOGCR-SWL");
        BOOST_CHECK_EQUAL(columns[4], "SGU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(All_Good_2)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.4f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.2f;
        endPoints.Sgu = 0.6;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.4f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.2f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.4f, 1.0e-5f);
        BOOST_CHECK_CLOSE(values[4], 0.6f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        constexpr auto expectNumExportedCheckValues = std::size_t{5};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        auto check = Checks::DisplacingOil_GO<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = 0.7f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sogcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sgu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sogcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sgcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sogcr-Swl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Sgu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        constexpr auto expectNumExportedCheckValues = std::size_t{5};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        auto check = Checks::DisplacingOil_GO<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = 0.7f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = 0.3f;
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = 0.25f;
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Swl = std::numeric_limits<float>::infinity();
        endPoints.Sgcr = std::numeric_limits<float>::infinity();
        endPoints.Sogcr = std::numeric_limits<float>::infinity();
        endPoints.Sgu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sogcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sgcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sogcr-Swl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Sgu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Sr_TooSmall)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.55f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.35f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.55f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.35f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.1f, 3.0e-5f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_Is_Lower_Bound)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.50f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.35f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.50f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.35f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.15f, 2.0e-5f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_TooLarge)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.0f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.25f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_Is_Upper_Bound)
{
    auto check = Checks::DisplacingOil_GO<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Swl = 0.05f;
        endPoints.Sgcr = 0.15f;
        endPoints.Sogcr = 0.25f;
        endPoints.Sgu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.05f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.7f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Displacing_Oil_in_Gas_Oil

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Displacing_Oil_in_Oil_Water)

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto check = Checks::DisplacingOil_OW<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    BOOST_CHECK_EQUAL(check.numExportedCheckValues(), expectNumExportedCheckValues);

    {
        auto columns = std::vector<std::string>(expectNumExportedCheckValues);
        check.columnNames(columns.data());

        BOOST_CHECK_EQUAL(columns[0], "SGL");
        BOOST_CHECK_EQUAL(columns[1], "SOWCR");
        BOOST_CHECK_EQUAL(columns[2], "SWCR");
        BOOST_CHECK_EQUAL(columns[3], "1-SOWCR-SGL");
        BOOST_CHECK_EQUAL(columns[4], "SWU");
    }

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(! check.isViolated(), "Test must not be violated");
    BOOST_CHECK_MESSAGE(! check.isCritical(), "Test must not be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Non_Finite)
{
    // NaN
    if constexpr (std::numeric_limits<float>::has_quiet_NaN) {
        constexpr auto expectNumExportedCheckValues = std::size_t{5};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        auto check = Checks::DisplacingOil_OW<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = 0.7f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Sowcr = std::numeric_limits<float>::quiet_NaN();
        endPoints.Swu = std::numeric_limits<float>::quiet_NaN();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isnan(values[0]), "Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[1]), "Sowcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[2]), "Swcr value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[3]), "1-Sowcr-Sgl value must be NaN");
            BOOST_CHECK_MESSAGE(std::isnan(values[4]), "Swu value must be NaN");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }

    // Inf
    if constexpr (std::numeric_limits<float>::has_infinity) {
        constexpr auto expectNumExportedCheckValues = std::size_t{5};

        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        auto check = Checks::DisplacingOil_OW<float>{};
        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.7;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_CLOSE(values[3], 0.45f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = 0.7f;

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = 0.3f;
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_CLOSE(values[1], 0.3f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = 0.25f;
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_CLOSE(values[0], 0.25f, 1.0e-6f);
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }

        endPoints.Sgl = std::numeric_limits<float>::infinity();
        endPoints.Swcr = std::numeric_limits<float>::infinity();
        endPoints.Sowcr = std::numeric_limits<float>::infinity();
        endPoints.Swu = std::numeric_limits<float>::infinity();

        check.test(endPoints);

        {
            auto values = std::vector<float>(expectNumExportedCheckValues);
            check.exportCheckValues(values.data());

            BOOST_CHECK_MESSAGE(std::isinf(values[0]), "Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[1]), "Sowcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[2]), "Swcr value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[3]), "1-Sowcr-Sgl value must be Inf");
            BOOST_CHECK_MESSAGE(std::isinf(values[4]), "Swu value must be Inf");

            BOOST_CHECK_MESSAGE(check.isViolated(), "Check must be violated");
            BOOST_CHECK_MESSAGE(check.isCritical(), "Check must be violated at critical level");
        }
    }
}

BOOST_AUTO_TEST_CASE(Sr_TooSmall)
{
    auto check = Checks::DisplacingOil_OW<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.55f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.35f;
        endPoints.Swu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.55f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.35f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.1f, 3.0e-5f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_Is_Lower_Bound)
{
    auto check = Checks::DisplacingOil_OW<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.50f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.35f;
        endPoints.Swu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.50f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.35f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.15f, 2.0e-5f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_TooLarge)
{
    auto check = Checks::DisplacingOil_OW<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.0f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.25f;
        endPoints.Swu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_CASE(Sr_Is_Upper_Bound)
{
    auto check = Checks::DisplacingOil_OW<float>{};

    constexpr auto expectNumExportedCheckValues = std::size_t{5};

    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        endPoints.Sgl = 0.05f;
        endPoints.Swcr = 0.15f;
        endPoints.Sowcr = 0.25f;
        endPoints.Swu = 0.7;

        check.test(endPoints);
    }

    {
        auto values = std::vector<float>(expectNumExportedCheckValues);
        check.exportCheckValues(values.data());

        BOOST_CHECK_CLOSE(values[0], 0.05f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[1], 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[2], 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[3], 0.7f, 1.0e-6f);
        BOOST_CHECK_CLOSE(values[4], 0.7f, 1.0e-6f);
    }

    BOOST_CHECK_MESSAGE(check.isViolated(), "Test must be violated");
    BOOST_CHECK_MESSAGE(check.isCritical(), "Test must be violated at critical level");
}

BOOST_AUTO_TEST_SUITE_END() // Displacing_Oil_in_Oil_Water
