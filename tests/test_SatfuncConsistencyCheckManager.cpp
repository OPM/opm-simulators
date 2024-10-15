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

#define BOOST_TEST_MODULE TestSatfuncConsistencyCheckManager

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

#include <opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include <fmt/format.h>

#if HAVE_MPI
#include <array>
#include <iostream>
#endif  // HAVE_MPI

// ###########################################################################

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

    bool init_unit_test_func()
    {
        return true;
    }

    // -----------------------------------------------------------------------

    template <typename Scalar>
    using CheckMgr = Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>;

    template <typename Scalar>
    using ViolationLevel = typename CheckMgr<Scalar>::ViolationLevel;

    constexpr auto root = 0;
    constexpr auto numSamplePoints = std::size_t{1};
    const auto localToGlobal = [](const int) { return std::size_t{0}; };

    Opm::Deck makeEpsDeck(std::string_view epsSpec)
    {
        return Opm::Parser{}.parseString(fmt::format(R"(RUNSPEC
DIMENS
 1 1 1 /

OIL
GAS
WATER

TABDIMS
/

ENDSCALE
/

-- =================================================================
GRID

DXV
 100 /

DYV
 100 /

DZV
 5 /

DEPTHZ
 4*2000 /

PERMX
 100 /

PERMY
 100 /

PERMZ
 10 /

PORO
 0.3 /

-- =================================================================
PROPS

SGOF
 0.00  0.0  1.0 0.0
 0.80   1*  0.0 0.0
 0.85  1.0  0.0 0.0
 /

SWOF
 0.15  0.0  1.0 0.0
 0.80  0.8  0.0 0.0
 1.0   1.0  0.0 0.0
/
{}

-- =================================================================
REGIONS

SATNUM
 1 /

END
)", epsSpec));
    }

    std::pair<Dune::CpGrid, Opm::EclipseState>
    setup(std::string_view epsSpec)
    {
        auto ret = std::pair<Dune::CpGrid, Opm::EclipseState> {
            std::piecewise_construct,
            std::forward_as_tuple(),
            std::forward_as_tuple(makeEpsDeck(epsSpec))
        };

        auto& [cpgrid, es] = ret;

        cpgrid.processEclipseFormat(&es.getInputGrid(), &es,
                                    /* periodic_extension = */ false,
                                    /* turn_normals = */       false,
                                    /* clip_z = */             false,
                                    /* pinchActive = */        false);

        return ret;
    }

} // Anonymous namespace

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Water_Phase)

BOOST_AUTO_TEST_CASE(SWL_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWL
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum water saturation
  0 <= SWL < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SWL           |
+------------+---------------+
| (1, 1, 1)  | -5.000000e-02 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SWL_Too_High)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWL
 1.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SWL + SGU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  8.500000e-01 |  1.850000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in G/O system at minimum gas saturation
  SOGCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOGCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  0.000000e+00 |  5.000000e-02 |  0.000000e+00 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in G/O system at critical gas saturation
  SOGCR < 1 - SWL - SGCR
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SWL           | SGCR          | SOGCR         | 1 - SWL - SGCR |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  1.000000e+00 |  0.000000e+00 |  5.000000e-02 |   0.000000e+00 |
+------------+---------------+---------------+---------------+----------------+


Consistency Problem:
  Mobile oil saturation in O/W system at minimum water saturation
  SOWCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOWCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  0.000000e+00 |  2.000000e-01 |  0.000000e+00 |
+------------+---------------+---------------+---------------+---------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum water saturation
  0 <= SWL < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SWL           |
+------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |
+------------+---------------+


Consistency Problem:
  Mobile water saturation
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  1.500000e-01 |  1.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SWU_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWU
 0.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Positive maximum water saturation
  0 < SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SWU           |
+------------+---------------+
| (1, 1, 1)  |  0.000000e+00 |
+------------+---------------+


Consistency Problem:
  Mobile water saturation
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  1.500000e-01 |  0.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SWU_Too_High)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWU
 1.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SGL + SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  0.000000e+00 |  1.050000e+00 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Positive maximum water saturation
  0 < SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SWU           |
+------------+---------------+
| (1, 1, 1)  |  1.050000e+00 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SWCR_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWCR
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile water saturation
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 | -5.000000e-02 |  1.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SWCR_Too_High)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWCR
 1.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in O/W system at critical water saturation
  SOWCR < 1 - SWCR - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SGL           | SWCR          | SOWCR         | 1 - SWCR - SGL |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  0.000000e+00 |  1.050000e+00 |  2.000000e-01 |  -5.000000e-02 |
+------------+---------------+---------------+---------------+----------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile water saturation
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  1.050000e+00 |  1.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END() // Water_Phase

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Gas_Phase)

BOOST_AUTO_TEST_CASE(SGL_Too_Low)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGL
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum gas saturation
  0 <= SGL < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SGL           |
+------------+---------------+
| (1, 1, 1)  | -5.000000e-02 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SGL_Too_High)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGL
 1.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in G/O system at minimum gas saturation
  SOGCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOGCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  1.000000e+00 |  5.000000e-02 | -1.500000e-01 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SGL + SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  1.000000e+00 |  2.000000e+00 |
+------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in O/W system at minimum water saturation
  SOWCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOWCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  1.000000e+00 |  2.000000e-01 | -1.500000e-01 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in O/W system at critical water saturation
  SOWCR < 1 - SWCR - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SGL           | SWCR          | SOWCR         | 1 - SWCR - SGL |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  1.000000e+00 |  1.500000e-01 |  2.000000e-01 |  -1.500000e-01 |
+------------+---------------+---------------+---------------+----------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum gas saturation
  0 <= SGL < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SGL           |
+------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |
+------------+---------------+


Consistency Problem:
  Mobile gas saturation
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |  0.000000e+00 |  8.500000e-01 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SGU_Too_Low)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGU
 0.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Positive maximum gas saturation must not exceed one
  0 < SGU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SGU           |
+------------+---------------+
| (1, 1, 1)  |  0.000000e+00 |
+------------+---------------+


Consistency Problem:
  Mobile gas saturation
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  0.000000e+00 |  0.000000e+00 |  0.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SGU_Too_High)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGU
 1.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SWL + SGU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  1.050000e+00 |  1.200000e+00 |
+------------+---------------+---------------+---------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Positive maximum gas saturation must not exceed one
  0 < SGU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SGU           |
+------------+---------------+
| (1, 1, 1)  |  1.050000e+00 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SGCR_Too_Low)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGCR
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile gas saturation
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  0.000000e+00 | -5.000000e-02 |  8.500000e-01 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SGCR_Too_High)
{
    using Scalar = float;

    const auto& [grid, es] = setup(R"(
SGCR
 1.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in G/O system at critical gas saturation
  SOGCR < 1 - SWL - SGCR
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SWL           | SGCR          | SOGCR         | 1 - SWL - SGCR |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  1.500000e-01 |  1.050000e+00 |  5.000000e-02 |  -1.999999e-01 |
+------------+---------------+---------------+---------------+----------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile gas saturation
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  0.000000e+00 |  1.050000e+00 |  8.500000e-01 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END() // Gas_Phase

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Oil_Phase)

BOOST_AUTO_TEST_CASE(So_At_SGMax_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SGU
  0.9 /
SWL
  0.15 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedCriticalChecks(),
                        "There must not be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SWL + SGU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGU           | SWL + SGU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  9.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(So_At_SWMax_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWU
  0.9 /
SGL
  0.15 /
SGCR
  0.15 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedCriticalChecks(),
                        "There must not be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SGL + SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  9.000000e-01 |  1.050000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SOWCR_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SOWCR
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative critical oil saturation in O/W system
  0 <= SOWCR < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SOWCR         |
+------------+---------------+
| (1, 1, 1)  | -5.000000e-02 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SOWCR_Too_High)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SOWCR
 1.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in O/W system at minimum water saturation
  SOWCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOWCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  0.000000e+00 |  1.000000e+00 |  8.500000e-01 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in O/W system at critical water saturation
  SOWCR < 1 - SWCR - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SGL           | SWCR          | SOWCR         | 1 - SWCR - SGL |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  0.000000e+00 |  1.500000e-01 |  1.000000e+00 |   8.500000e-01 |
+------------+---------------+---------------+---------------+----------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative critical oil saturation in O/W system
  0 <= SOWCR < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SOWCR         |
+------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SOGCR_Too_Low)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SOGCR
 -0.05 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative critical oil saturation in G/O system
  0 <= SOGCR < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SOGCR         |
+------------+---------------+
| (1, 1, 1)  | -5.000000e-02 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(SOGCR_Too_High)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SOGCR
 1.0 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in G/O system at minimum gas saturation
  SOGCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOGCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  0.000000e+00 |  1.000000e+00 |  8.500000e-01 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Mobile oil saturation in G/O system at critical gas saturation
  SOGCR < 1 - SWL - SGCR
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SWL           | SGCR          | SOGCR         | 1 - SWL - SGCR |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  1.500000e-01 |  0.000000e+00 |  1.000000e+00 |   8.500000e-01 |
+------------+---------------+---------------+---------------+----------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Non-negative critical oil saturation in G/O system
  0 <= SOGCR < 1
  Total Violations: 1

List of Violations
+------------+---------------+
| Grid Block | SOGCR         |
+------------+---------------+
| (1, 1, 1)  |  1.000000e+00 |
+------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Mobile_Oil_OW_System)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWCR
  0.42 /
SOWCR
  0.63 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedCriticalChecks(),
                        "There must not be any failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in O/W system at critical water saturation
  SOWCR < 1 - SWCR - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SGL           | SWCR          | SOWCR         | 1 - SWCR - SGL |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  0.000000e+00 |  4.200000e-01 |  6.300000e-01 |   5.800000e-01 |
+------------+---------------+---------------+---------------+----------------+


)");
}

BOOST_AUTO_TEST_CASE(Mobile_Oil_OW_System_At_SWL)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SWL
  0.42 /
SGU
  0.57 /
SOWCR
  0.63 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in O/W system at minimum water saturation
  SOWCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOWCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  4.200000e-01 |  0.000000e+00 |  6.300000e-01 |  5.800000e-01 |
+------------+---------------+---------------+---------------+---------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile water saturation
  SWL <= SWCR < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SWL           | SWCR          | SWU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  4.200000e-01 |  1.500000e-01 |  1.000000e+00 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Mobile_Oil_GO_System)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SGCR
  0.27 /
SOGCR
  0.63 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedCriticalChecks(),
                        "There must not be any failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in G/O system at critical gas saturation
  SOGCR < 1 - SWL - SGCR
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+----------------+
| Grid Block | SWL           | SGCR          | SOGCR         | 1 - SWL - SGCR |
+------------+---------------+---------------+---------------+----------------+
| (1, 1, 1)  |  1.500000e-01 |  2.700000e-01 |  6.300000e-01 |   5.800000e-01 |
+------------+---------------+---------------+---------------+----------------+


)");
}

BOOST_AUTO_TEST_CASE(Mobile_Oil_GO_System_At_SGL)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SGL
  0.27 /
SOGCR
  0.63 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(checkMgr.anyFailedStandardChecks(),
                        "There must be failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Standard,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile oil saturation in G/O system at minimum gas saturation
  SOGCR < 1 - SWL - SGL
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SGL           | SOGCR         | 1 - SWL - SGL |
+------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  2.700000e-01 |  6.300000e-01 |  5.800000e-01 |
+------------+---------------+---------------+---------------+---------------+


Consistency Problem:
  Non-negative minimum oil saturation in G/O system
  SGL + SWU <= 1
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SWU           | SGL + SWU     |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  2.700000e-01 |  1.000000e+00 |  1.270000e+00 |
+------------+---------------+---------------+---------------+


)");

    msg.clear();
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile gas saturation
  SGL <= SGCR < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+
| Grid Block | SGL           | SGCR          | SGU           |
+------------+---------------+---------------+---------------+
| (1, 1, 1)  |  2.700000e-01 |  0.000000e+00 |  8.500000e-01 |
+------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Displacing_Oil_OW_System)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SCALECRS
 'YES' /
SWCR
  0.4 /
SWU
  0.6 /
SGL
  0.15 /
SGCR
  0.15 /
SOWCR
  0.2 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile displacing oil in three point horizontally scaled oil/water system
  SWCR < 1-SOWCR-SGL < SWU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SGL           | SOWCR         | SWCR          | 1-SOWCR-SGL   | SWU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.500000e-01 |  2.000000e-01 |  4.000000e-01 |  6.500000e-01 |  6.000000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_CASE(Displacing_Oil_GO_System)
{
    using Scalar = double;

    const auto& [grid, es] = setup(R"(
SCALECRS
 'YES' /
SWL
  0.1 /
SWCR
  0.2 /
SGCR
  0.2 /
SGU
  0.6 /
SOGCR
  0.2 /
)");

    auto checkMgr = CheckMgr<Scalar>(numSamplePoints, es, localToGlobal);

    checkMgr.collectFailuresTo(root)
        .run(grid.leafGridView(), [](const auto&) { return 0; });

    BOOST_CHECK_MESSAGE(! checkMgr.anyFailedStandardChecks(),
                        "There must not be any failed checks at the standard level");
    BOOST_CHECK_MESSAGE(checkMgr.anyFailedCriticalChecks(),
                        "There must be failed checks at the critical level");

    auto msg = std::string{};
    checkMgr.reportFailures(ViolationLevel<Scalar>::Critical,
                            [&msg](std::string_view record)
                            { msg += fmt::format("{}\n", record); });

    BOOST_CHECK_EQUAL(msg, R"(Consistency Problem:
  Mobile displacing oil in three point horizontally scaled gas/oil system
  SGCR < 1-SOGCR-SWL < SGU
  Total Violations: 1

List of Violations
+------------+---------------+---------------+---------------+---------------+---------------+
| Grid Block | SWL           | SOGCR         | SGCR          | 1-SOGCR-SWL   | SGU           |
+------------+---------------+---------------+---------------+---------------+---------------+
| (1, 1, 1)  |  1.000000e-01 |  2.000000e-01 |  2.000000e-01 |  7.000000e-01 |  6.000000e-01 |
+------------+---------------+---------------+---------------+---------------+---------------+


)");
}

BOOST_AUTO_TEST_SUITE_END() // Oil_Phase

// ===========================================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    register_error_handler();

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
