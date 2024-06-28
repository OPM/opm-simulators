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

#define BOOST_TEST_MODULE TestSatfuncConsistencyChecks_Parallel

#define BOOST_TEST_NO_MAIN

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

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

#include <fmt/format.h>

namespace {

#if HAVE_MPI
    struct MPIError
    {
        MPIError(std::string_view errstr, const int ec)
            : errorstring { errstr }
            , errorcode   { ec }
        {}

        std::string errorstring;
        int errorcode;
    };

    void MPI_err_handler(MPI_Comm*, int* err_code, ...)
    {
        std::array<char, MPI_MAX_ERROR_STRING> err_string_vec{'\0'};
        auto err_length = 0;

        MPI_Error_string(*err_code, err_string_vec.data(), &err_length);

        auto err_string = std::string_view {
            err_string_vec.data(),
            static_cast<std::string_view::size_type>(err_length)
        };

        std::cerr << "An MPI Error ocurred:\n  -> " << err_string << '\n';

        throw MPIError { err_string, *err_code };
    }

    // Register a throwing error handler to allow for debugging with
    //
    //   catch throw
    //
    // in GDB.
    void register_error_handler()
    {
        MPI_Errhandler handler{};

        MPI_Comm_create_errhandler(MPI_err_handler, &handler);
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
    }

#else // !HAVE_MPI

    void register_error_handler()
    {}

#endif // HAVE_MPI

    class NProc_Is
    {
    public:
        explicit NProc_Is(const int expectNP)
            : expectNP_ { expectNP }
        {}

        boost::test_tools::assertion_result
        operator()(boost::unit_test::test_unit_id) const
        {
            auto comm = Opm::Parallel::Communication {
                Dune::MPIHelper::getCommunicator()
            };

            if (comm.size() == this->expectNP_) {
                return true;
            }

            boost::test_tools::assertion_result response(false);
            response.message() << "Number of MPI processes ("
                               << comm.size()
                               << ") differs from expected "
                               << this->expectNP_;

            return response;
        }

    private:
        int expectNP_{};
    };

    bool init_unit_test_func()
    {
        return true;
    }

} // Anonymous namespace

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
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 1};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Standard>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank(), makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(! checker.anyFailedChecks(),
                        "There must be no failed checks");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    auto rpt = std::string{};
    if (comm.rank() == 0) {
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });
    }

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

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Standard_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<StandardViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point
  0 <= SWL < 1
  Total Violations: 2

List of Violations
+------+---------------+
| Cell | SWL           |
+------+---------------+
| 1    |  1.729000e+01 |
| 2    |  1.729000e+01 |
+------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Critical_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"PVTNUM", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<CriticalViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(42 * (comm.rank() + 1), makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(! checker.anyFailedChecks(),
                        "There must be no failed standard level checks");

    BOOST_CHECK_MESSAGE(checker.anyFailedCriticalChecks(),
                        "There must be at least one failed Critical check");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Critical,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Minimum Pressure
  PRESS > 350
  Total Violations: 2

List of Violations
+--------+---------------+
| PVTNUM | PRESS         |
+--------+---------------+
| 42     |  3.141593e+02 |
| 84     |  3.141593e+02 |
+--------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Standard_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<StandardViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point
  0 <= SWL < 1
  Total Violations: 3

List of Violations
+------+---------------+
| Cell | SWL           |
+------+---------------+
| 1    |  1.729000e+01 |
| 2    |  1.729000e+01 |
| 3    |  1.729000e+01 |
+------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Critical_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"PVTNUM", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<CriticalViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(42 * (comm.rank() + 1), makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(! checker.anyFailedChecks(),
                        "There must be no failed standard level checks");

    BOOST_CHECK_MESSAGE(checker.anyFailedCriticalChecks(),
                        "There must be at least one failed Critical check");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Critical,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Minimum Pressure
  PRESS > 350
  Total Violations: 3

List of Violations
+--------+---------------+
| PVTNUM | PRESS         |
+--------+---------------+
| 42     |  3.141593e+02 |
| 84     |  3.141593e+02 |
| 126    |  3.141593e+02 |
+--------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Standard_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"Cell", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<StandardViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point
  0 <= SWL < 1
  Total Violations: 4

List of Violations
+------+---------------+
| Cell | SWL           |
+------+---------------+
| 1    |  1.729000e+01 |
| 2    |  1.729000e+01 |
| 3    |  1.729000e+01 |
| 4    |  1.729000e+01 |
+------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Critical_Violation)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<double>{"PVTNUM", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<CriticalViolation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(42 * (comm.rank() + 1), makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(! checker.anyFailedChecks(),
                        "There must be no failed standard level checks");

    BOOST_CHECK_MESSAGE(checker.anyFailedCriticalChecks(),
                        "There must be at least one failed Critical check");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<double>::ViolationLevel::Critical,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Minimum Pressure
  PRESS > 350
  Total Violations: 4

List of Violations
+--------+---------------+
| PVTNUM | PRESS         |
+--------+---------------+
| 42     |  3.141593e+02 |
| 84     |  3.141593e+02 |
| 126    |  3.141593e+02 |
| 168    |  3.141593e+02 |
+--------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

// ---------------------------------------------------------------------------

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

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Bucket", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Sum
  a + 1/2 < 2
  Total Violations: 2

List of Violations
+--------+---------------+---------------+
| Bucket | a             | a + 1/2       |
+--------+---------------+---------------+
| 1      |  1.600000e+00 |  2.100000e+00 |
| 2      |  1.600000e+00 |  2.100000e+00 |
+--------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Bucket", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Sum
  a + 1/2 < 2
  Total Violations: 3

List of Violations
+--------+---------------+---------------+
| Bucket | a             | a + 1/2       |
+--------+---------------+---------------+
| 1      |  1.600000e+00 |  2.100000e+00 |
| 2      |  1.600000e+00 |  2.100000e+00 |
| 3      |  1.600000e+00 |  2.100000e+00 |
+--------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Bucket", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Sum
  a + 1/2 < 2
  Total Violations: 4

List of Violations
+--------+---------------+---------------+
| Bucket | a             | a + 1/2       |
+--------+---------------+---------------+
| 1      |  1.600000e+00 |  2.100000e+00 |
| 2      |  1.600000e+00 |  2.100000e+00 |
| 3      |  1.600000e+00 |  2.100000e+00 |
| 4      |  1.600000e+00 |  2.100000e+00 |
+--------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

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

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point Displacing Saturation
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 2

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 2          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    checker.collectFailures(0, comm);

    if (comm.rank() == 0) {
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
  Total Violations: 8

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point Displacing Saturation
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 3

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 2          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 3          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    checker.collectFailures(0, comm);

    if (comm.rank() == 0) {
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
  Total Violations: 12

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Water Phase End-Point Displacing Saturation
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 2          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 3          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| 4          |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<Violation>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    checker.collectFailures(0, comm);

    if (comm.rank() == 0) {
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
  Total Violations: 16

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (1, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (2, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (3, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (4, 1234)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (4, 1618)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (4, 1729)  |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
| (4, 31415) |  1.000000e-01 |  7.000000e-01 |  3.000000e-01 |  2.000000e-01 |  8.800000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

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

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Gas Phase End-Point Monotonicity
  SGL <= SGCR < SGU
  Total Violations: 2

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| 1          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 2          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 2

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| 1          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 2          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 2

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 2          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 2

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 2          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
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
  Total Violations: 8

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 8

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 8

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 8

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Gas Phase End-Point Monotonicity
  SGL <= SGCR < SGU
  Total Violations: 3

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| 1          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 2          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 3          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 3

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| 1          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 2          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 3          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 3

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 2          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 3          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 3

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 2          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 3          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
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
  Total Violations: 12

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 12

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 12

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 12

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Standard)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 4};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    checker.checkEndpoints(comm.rank() + 1, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
        auto rpt = std::string{};
        checker.reportFailures(Opm::SatfuncConsistencyChecks<float>::ViolationLevel::Standard,
                               [&rpt](std::string_view record)
                               {
                                   rpt += fmt::format("{}\n", record);
                               });

        BOOST_CHECK_EQUAL(rpt, R"(Consistency Problem:
  Gas Phase End-Point Monotonicity
  SGL <= SGCR < SGU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| 1          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 2          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 3          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| 4          |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| 1          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 2          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 3          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| 4          |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 2          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 3          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| 4          |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 4

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| 1          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 2          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 3          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| 4          |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_CASE(Standard_Multiple_Failing_Points)
{
    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto checker = Opm::SatfuncConsistencyChecks<float>{"Grid Block", 16};

    checker.resetCheckSet();
    checker.addCheck(std::make_unique<MonotoneGas>());
    checker.addCheck(std::make_unique<NonNegOSAT_GO>());
    checker.addCheck(std::make_unique<MonotoneWater>());
    checker.addCheck(std::make_unique<NonNegOSAT_OW>());
    checker.finaliseCheckSet();

    const auto rankMultiplier = 100'000;

    checker.setPointIDFormatCallback([rankMultiplier](const std::size_t pointID)
    {
        const auto rank = pointID / rankMultiplier;
        const auto pt   = pointID % rankMultiplier;

        return fmt::format("({}, {})", rank, pt);
    });

    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1234, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1729, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) +  1618, makePoints());
    checker.checkEndpoints(rankMultiplier*(comm.rank() + 1) + 31415, makePoints());

    checker.collectFailures(0, comm);

    BOOST_CHECK_MESSAGE(checker.anyFailedChecks(),
                        "There must be at least one failed check");

    BOOST_CHECK_MESSAGE(! checker.anyFailedCriticalChecks(),
                        "There must be no failed critical checks");

    if (comm.rank() == 0) {
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
  Total Violations: 16

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (1, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (2, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (3, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (4, 1234)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (4, 1618)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (4, 1729)  |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
| (4, 31415) |  0.000000e+00 | -1.000000e-01 |  8.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (G/O)
  SWL + SGU <= 1
  Total Violations: 16

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (1, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (2, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (3, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (4, 1234)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (4, 1618)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (4, 1729)  |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
| (4, 31415) |  2.500000e-01 |  8.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Water Phase End-Point Monotonicity
  SWL <= SWCR < SWU
  Total Violations: 16

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (1, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (2, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (3, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (4, 1234)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (4, 1618)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (4, 1729)  |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
| (4, 31415) |  1.000000e-01 |  3.000000e-01 |  3.000000e-01 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Oil Phase Non-Negative Saturation (O/W)
  SGL + SWU <= 1
  Total Violations: 16

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (1, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (2, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (3, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (4, 1234)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (4, 1618)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (4, 1729)  |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
| (4, 31415) |  1.000000e-01 |  1.000000e+00 |  1.100000e+00 |
+------------+---------------+---------------+---------------+


)");
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()     // Multiple_Failing_Tests

// ===========================================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    register_error_handler();

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
