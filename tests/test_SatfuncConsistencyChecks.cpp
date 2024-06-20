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

#define BOOST_TEST_MODULE TestSatfuncConsistencyChecks

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

#include <fmt/format.h>

BOOST_AUTO_TEST_SUITE(NoFailures)

namespace {
    class Standard : public Opm::SatfuncConsistencyChecks<double>::Check
    {
    public:
        void test(const Opm::EclEpsScalingPointsInfo<double>&) override {}
        bool isViolated() const override { return false; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 1; }

        void exportCheckValues(double* exportedCheckValues) const override
        {
            *exportedCheckValues = 17.29;
        }

        std::string description() const override
        {
            return "Water Phase End-Point";
        }

        std::string condition() const override
        {
            return "0 <= SWL < 1";
        }

        void columnNames(std::string* headers) const override
        {
            *headers = "SWL";
        }
    };

    Opm::EclEpsScalingPointsInfo<double> makePoints() { return {}; }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(All_Good)
{
    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Standard>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(! checker.anyFailedChecks(),
                        "There must be no failed checks");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_MESSAGE(rpt.empty(), "There must be no output from reportFailures()");
}

BOOST_AUTO_TEST_SUITE_END()

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Single_Exported_Value)

namespace {
    class StandardViolation : public Opm::SatfuncConsistencyChecks<double>::Check
    {
    public:
        void test(const Opm::EclEpsScalingPointsInfo<double>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 1; }

        void exportCheckValues(double* exportedCheckValues) const override
        {
            *exportedCheckValues = 17.29;
        }

        std::string description() const override
        {
            return "Water Phase End-Point";
        }

        std::string condition() const override
        {
            return "0 <= SWL < 1";
        }

        void columnNames(std::string* headers) const override
        {
            *headers = "SWL";
        }
    };

    class CriticalViolation : public Opm::SatfuncConsistencyChecks<double>::Check
    {
    public:
        void test(const Opm::EclEpsScalingPointsInfo<double>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return true; }
        std::size_t numExportedCheckValues() const override { return 1; }

        void exportCheckValues(double* exportedCheckValues) const override
        {
            *exportedCheckValues = 314.15926;
        }

        std::string description() const override
        {
            return "Minimum Pressure";
        }

        std::string condition() const override
        {
            return "PRESS > 350";
        }

        void columnNames(std::string* headers) const override
        {
            *headers = "PRESS";
        }
    };

    Opm::EclEpsScalingPointsInfo<double> makePoints() { return {}; }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Standard_Violation)
{
    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<StandardViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point
  0 <= SWL < 1
  Total Violations: 1

List of Violations
+------+---------------+
| Cell | SWL           |
+------+---------------+
| 1234 |  1.729000e+01 |
+------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Standard_Violation_ReportIJK)
{
    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<StandardViolation>());
    checker.finaliseCheckSet();

    checker.setPointIDFormatCallback([](const std::size_t)
    {
        return std::string { "(11, 22, 33)" };
    });

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point
  0 <= SWL < 1
  Total Violations: 1

List of Violations
+--------------+---------------+
| Cell         | SWL           |
+--------------+---------------+
| (11, 22, 33) |  1.729000e+01 |
+--------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Critical_Violation)
{
    auto checker = Opm::SatfuncConsistencyChecks<double>{"PVTNUM", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<CriticalViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(42, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedCriticalChecks(),
                        "There must be at least one failed Critical check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Critical,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Minimum Pressure
  PRESS > 350
  Total Violations: 1

List of Violations
+--------+---------------+
| PVTNUM | PRESS         |
+--------+---------------+
| 42     |  3.141593e+02 |
+--------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END()     // Single_Exported_Value

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Two_Exported_Values)

namespace {
    class Violation : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 2; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] = 1.6f;
            exportedCheckValues[1] = 1.6f + 0.5f;
        }

        std::string description() const override
        {
            return "Sum";
        }

        std::string condition() const override
        {
            return "a + 1/2 < 2";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "a";
            headers[1] = "a + 1/2";
        }
    };

    Opm::EclEpsScalingPointsInfo<float> makePoints() { return {}; }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Standard)
{
    auto checker = Opm::SatfuncConsistencyChecks<float>{"Bucket", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Sum
  a + 1/2 < 2
  Total Violations: 1

List of Violations
+--------+---------------+---------------+
| Bucket | a             | a + 1/2       |
+--------+---------------+---------------+
| 1234   |  1.600000e+00 |  2.100000e+00 |
+--------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END()     // Two_Exported_Values

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Five_Exported_Values)

namespace {
    class Violation : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 5; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] = 0.1f;
            exportedCheckValues[1] = 0.7f;
            exportedCheckValues[2] = 0.3f;
            exportedCheckValues[3] = 0.2f;
            exportedCheckValues[4] = 0.88f;
        }

        std::string description() const override
        {
            return "Water Phase End-Point Displacing Saturation";
        }

        std::string condition() const override
        {
            return "SWCR < 1-SOWCR-SGL < SWU";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SOWCR";
            headers[2] = "SWCR";
            headers[3] = "1-SOWCR-SGL";
            headers[4] = "SWU";
        }
    };

    Opm::EclEpsScalingPointsInfo<float> makePoints() { return {}; }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Standard)
{
    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.setPointIDFormatCallback([](const std::size_t)
    {
        return std::string { "(121, 323, 42)" };
    });

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point Displacing Saturation
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 1

List of Violations
+----------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block     | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+----------------+---------------+---------------+---------------+---------------+---------------+
| (121, 323, 42) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+----------------+---------------+---------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 5};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints( 1234, makePoints());
    checker.checkEndpoints( 1729, makePoints());
    checker.checkEndpoints( 1618, makePoints());
    checker.checkEndpoints(31415, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    // Note that grid blocks are reported in sorted order rather in the
    // order of insertion.
    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point Displacing Saturation
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| 1234       |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 1618       |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 1729       |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 31415      |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END()     // Five_Exported_Values

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Multiple_Failing_Tests)

namespace {
    class MonotoneWater : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 3; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] = 0.1f;
            exportedCheckValues[1] = 0.3f;
            exportedCheckValues[2] = 0.3f;
        }

        std::string description() const override
        {
            return "Water Phase End-Point Monotonicity";
        }

        std::string condition() const override
        {
            return "SWL <= SWCR < SWU";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SWCR";
            headers[2] = "SWU";
        }
    };

    class MonotoneGas : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 3; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] =  0.0f;
            exportedCheckValues[1] = -0.1f;
            exportedCheckValues[2] =  0.8f;
        }

        std::string description() const override
        {
            return "Gas Phase End-Point Monotonicity";
        }

        std::string condition() const override
        {
            return "SGL <= SGCR < SGU";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SGCR";
            headers[2] = "SGU";
        }
    };

    class NonNegOSAT_OW : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 3; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] = 0.1f;
            exportedCheckValues[1] = 1.0f;
            exportedCheckValues[2] = 1.1f;
        }

        std::string description() const override
        {
            return "Oil Phase Non-Negative Saturation (O/W)";
        }

        std::string condition() const override
        {
            return "SGL + SWU <= 1";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SWU";
            headers[2] = "SGL + SWU";
        }
    };

    class NonNegOSAT_GO : public Opm::SatfuncConsistencyChecks<float>::Check
    {
        void test(const Opm::EclEpsScalingPointsInfo<float>&) override {}
        bool isViolated() const override { return true; }
        bool isCritical() const override { return false; }
        std::size_t numExportedCheckValues() const override { return 3; }

        void exportCheckValues(float* exportedCheckValues) const override
        {
            exportedCheckValues[0] = 0.25f;
            exportedCheckValues[1] = 0.8f;
            exportedCheckValues[2] = 1.05f;
        }

        std::string description() const override
        {
            return "Oil Phase Non-Negative Saturation (G/O)";
        }

        std::string condition() const override
        {
            return "SWL + SGU <= 1";
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGU";
            headers[2] = "SWL + SGU";
        }
    };

    Opm::EclEpsScalingPointsInfo<float> makePoints() { return {}; }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Standard)
{
    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    checker.setPointIDFormatCallback([](const std::size_t)
    {
        return std::string { "(121, 323, 42)" };
    });

    checker.checkEndpoints(1234, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Gas Phase End-Point Monotonicity
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+----------------+---------------+---------------+---------------+
| Grid Block     | SGL           | SGCR          | SGU           |
+----------------+---------------+---------------+---------------+
| (121, 323, 42) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+----------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 1

List of Violations
+----------------+---------------+---------------+---------------+
| Grid Block     | SWL           | SGU           | SWL + SGU     |
+----------------+---------------+---------------+---------------+
| (121, 323, 42) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+----------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+----------------+---------------+---------------+---------------+
| Grid Block     | SWL           | SWCR          | SWU           |
+----------------+---------------+---------------+---------------+
| (121, 323, 42) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+----------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 1

List of Violations
+----------------+---------------+---------------+---------------+
| Grid Block     | SGL           | SWU           | SGL + SWU     |
+----------------+---------------+---------------+---------------+
| (121, 323, 42) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+----------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 5};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    checker.checkEndpoints( 1234, makePoints());
    checker.checkEndpoints( 1729, makePoints());
    checker.checkEndpoints( 1618, makePoints());
    checker.checkEndpoints(31415, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    auto rpt = std::string{};
    checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                           [&rpt](std::string_view record)
                           {
                               rpt += fmt::format("{}\n", record);
                           });

    // Note that grid blocks are reported in sorted order rather in the
    // order of insertion.
    BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Gas Phase End-Point Monotonicity
  SGL <= SGCR < SGU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| 1234       |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 1618       |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 1729       |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 31415      |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| 1234       |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 1618       |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 1729       |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 31415      |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| 1234       |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 1618       |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 1729       |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 31415      |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| 1234       |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 1618       |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 1729       |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 31415      |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END()     // Multiple_Failing_Tests
