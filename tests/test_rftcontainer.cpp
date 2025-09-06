/*
  Copyright 2025 Equinor AS

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

#define BOOST_TEST_MODULE Test_RFTContainer_Parallel

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

#include <opm/simulators/flow/RFTContainer.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/output/data/Wells.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

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

    using TestFluidSystem = Opm::BlackOilFluidSystem
        <double, Opm::BlackOilDefaultFluidSystemIndices>;

    using RFTContainer = Opm::RFTContainer<TestFluidSystem>;

    struct Case
    {
        explicit Case(const std::string& spec)
            : Case { Opm::Parser{}.parseString(spec) }
        {}

        explicit Case(const Opm::Deck& deck)
            : es    { deck }
            , sched { deck, es, std::make_shared<Opm::Python>() }
        {}

        Opm::EclipseState es{};
        Opm::Schedule     sched{};
    };

    Case singleCentreWell()
    {
        return Case { R"(RUNSPEC
DIMENS
 3 3 3 /
START
  3 SEP 2025 /
TABDIMS
/
OIL
GAS
WATER
GRID
DXV
 3*100 /
DYV
 3*100 /
DZV
 3*5 /
DEPTHZ
 16*2000 /
EQUALS
  PERMX 100 /
  PERMY 100 /
  PERMZ  10 /
  PORO   0.3 /
/
PROPS
DENSITY
  800 1000 0.5 /
SCHEDULE
WELSPECS
  R G 2 2 1* OIL /
/
COMPDAT
  R 1* 1* 1 3 OPEN /
/
WCONPROD
  R SHUT ORAT 0.0 /
/
WRFTPLT
  R 'YES' /
/
DATES
  10 SEP 2025 /
  20 SEP 2025 /
  30 SEP 2025 /
/
END
)" };
    }

    Case quarterFiveSpot()
    {
        return Case { R"(RUNSPEC
DIMENS
 3 3 3 /
START
  3 SEP 2025 /
TABDIMS
/
OIL
GAS
WATER
GRID
DXV
 3*100 /
DYV
 3*100 /
DZV
 3*5 /
DEPTHZ
 16*2000 /
EQUALS
  PERMX 100 /
  PERMY 100 /
  PERMZ  10 /
  PORO   0.3 /
/
PROPS
DENSITY
  800 1000 0.5 /
SCHEDULE
WELSPECS
  R1 G 1 1 1* OIL /
  R2 G 3 3 1* OIL /
/
COMPDAT
  R* 1* 1* 1 3 OPEN /
/
WCONPROD
  R* SHUT ORAT 0.0 /
/
WRFTPLT
  R* 'YES' /
/
DATES
  10 SEP 2025 /
  20 SEP 2025 /
  30 SEP 2025 /
/
END
)" };
    }
} // Anonymous namespace

BOOST_AUTO_TEST_SUITE(Single_Rank)

BOOST_AUTO_TEST_CASE(Single_Well)
{
    const auto cse = singleCentreWell();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string&) { return r == 0; },
        [r = comm.rank()](const std::string&) { return r == 0; }
    };

    rftc.allocate(0);

    // IJK = (2,2,1)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 0),
                /* Po = */ []() { return 100;  },
                /* Sw = */ []() { return 0.01; },
                /* Sg = */ []() { return 0.8;  });

    // IJK = (2,2,2)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 1),
                /* Po = */ []() { return 150; },
                /* Sw = */ []() { return 0.1; },
                /* Sg = */ []() { return 0.5; });

    // IJK = (2,2,3)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 2),
                /* Po = */ []() { return 215; },
                /* Sw = */ []() { return 0.4; },
                /* Sg = */ []() { return 0.1; });

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    if (comm.rank() != 0) {
        return;
    }

    const auto rPos = xw.find("R");

    BOOST_REQUIRE_MESSAGE(rPos != xw.end(),
                          R"(There must be dynamic results for well "R")");

    BOOST_CHECK_EQUAL(rPos->second.connections.size(), std::size_t{3});

    // (2,2,1)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
    }

    // (2,2,2)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
    }

    // (2,2,3)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string&) { return r == 0; },
        [r = comm.rank()](const std::string&) { return r == 0; }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    // IJK = (1,1,1)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                /* Po = */ []() { return 100;  },
                /* Sw = */ []() { return 0.01; },
                /* Sg = */ []() { return 0.8;  });

    // IJK = (1,1,2)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                /* Po = */ []() { return 150; },
                /* Sw = */ []() { return 0.1; },
                /* Sg = */ []() { return 0.5; });

    // IJK = (1,1,3)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                /* Po = */ []() { return 215; },
                /* Sw = */ []() { return 0.4; },
                /* Sg = */ []() { return 0.1; });

    // --------------------------------------------------------------------

    // IJK = (3,3,1)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                /* Po = */ []() { return 123;  },
                /* Sw = */ []() { return 0.02; },
                /* Sg = */ []() { return 0.63; });

    // IJK = (3,3,2)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                /* Po = */ []() { return 152;  },
                /* Sw = */ []() { return 0.12; },
                /* Sg = */ []() { return 0.48; });

    // IJK = (3,3,3)
    rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                /* Po = */ []() { return 199;  },
                /* Sw = */ []() { return 0.42; },
                /* Sg = */ []() { return 0.07; });

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    if (comm.rank() != 0) {
        return;
    }

    // --------------------------------------------------------------------

    const auto r1Pos = xw.find("R1");

    BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                          R"(There must be dynamic results for well "R1")");

    BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

    // (1,1,1)
    {
        const auto* xc = r1Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R1" must have dynamic results in connection (1,1,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
    }

    // (1,1,2)
    {
        const auto* xc = r1Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R1" must have dynamic results in connection (1,1,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
    }

    // (1,1,3)
    {
        const auto* xc = r1Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R1" must have dynamic results in connection (1,1,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
    }

    // --------------------------------------------------------------------

    const auto r2Pos = xw.find("R2");

    BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                          R"(There must be dynamic results for well "R2")");

    BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

    // (3,3,1)
    {
        const auto* xc = r2Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R2" must have dynamic results in connection (3,3,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
    }

    // (3,3,2)
    {
        const auto* xc = r2Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R2" must have dynamic results in connection (3,3,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
    }

    // (3,3,3)
    {
        const auto* xc = r2Pos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R2" must have dynamic results in connection (3,3,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END()

// ===========================================================================

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Single_Split_Well)
{
    const auto cse = singleCentreWell();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string&) { return r == 0; },
        [](const std::string&) { return true; }
    };

    rftc.allocate(0);

    if (comm.rank() == 0) {
        // IJK = (2,2,1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });
    }
    else {
        // IJK = (2,2,2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (2,2,3)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    if (comm.rank() != 0) {
        return;
    }

    const auto rPos = xw.find("R");

    BOOST_REQUIRE_MESSAGE(rPos != xw.end(),
                          R"(There must be dynamic results for well "R")");

    BOOST_CHECK_EQUAL(rPos->second.connections.size(), std::size_t{3});

    // (2,2,1)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
    }

    // (2,2,2)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
    }

    // (2,2,3)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Non_Split_Wells)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 0) && (wname == "R2"))
                || ((r == 1) && (wname == "R1"));
        },
        [r = comm.rank()](const std::string& wname) {
            return ((r == 0) && (wname == "R2"))
                || ((r == 1) && (wname == "R1"));
        }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (1,1,1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (1,1,2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (1,1,3)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        // IJK = (3,3,1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });

        // IJK = (3,3,2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });

        // IJK = (3,3,3)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Split_Wells)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // Rank 0 owns R2, rank 1 ownes R1.  Both wells intersected on both
    // ranks.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 0) && (wname == "R2"))
                || ((r == 1) && (wname == "R1"));
        },
        [](const std::string&) { return true; }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });

        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });

        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ===========================================================================

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Single_Split_Well)
{
    const auto cse = singleCentreWell();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // Well owned by rank 2, but intersected on all ranks.
    const auto owner = 2;
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank(), owner](const std::string&) { return r == owner; },
        [](const std::string&) { return true; }
    };

    rftc.allocate(0);

    if (comm.rank() == 0) {
        // IJK = (2,2,1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });
    }
    else if (comm.rank() == 1) {
        // IJK = (2,2,2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });
    }
    else {
        // IJK = (2,2,3)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    if (comm.rank() != owner) {
        return;
    }

    const auto rPos = xw.find("R");

    BOOST_REQUIRE_MESSAGE(rPos != xw.end(),
                          R"(There must be dynamic results for well "R")");

    BOOST_CHECK_EQUAL(rPos->second.connections.size(), std::size_t{3});

    // (2,2,1)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
    }

    // (2,2,2)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
    }

    // (2,2,3)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Owned_And_Split)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // R1 is owned by rank 2 and intersected on rank 2.  R2 is owned by rank
    // 1 and intersected on ranks 0 and 1.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 2) && (wname == "R1"))
                || ((r == 1) && (wname == "R2"));
        },
        [r = comm.rank()](const std::string& wname) {
            return (((r == 0) || (r == 1)) && (wname == "R2"))
                || ((r == 2) && (wname == "R1"));
        }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });

        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });
    }

    if (comm.rank() == 0) {
        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Split_Wells)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // R1 is owned by rank 2.  R2 is owned by rank 1.  Both wells
    // intersected on all ranks.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 2) && (wname == "R1"))
                || ((r == 1) && (wname == "R2"));
        },
        [](const std::string&) { return true; }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });

        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    if (comm.rank() == 0) {
        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ===========================================================================

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Single_Split_Well)
{
    const auto cse = singleCentreWell();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // Well owned by rank 2, but intersected on ranks 0, 1, and 3.  Does
    // typically not happen in a real simulator run...
    const auto owner = 2;
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank(), owner](const std::string&) { return r == owner; },
        [r = comm.rank(), owner](const std::string&) { return r != owner; }
    };

    rftc.allocate(0);

    if (comm.rank() == 0) {
        // IJK = (2,2,1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });
    }
    else if (comm.rank() == 1) {
        // IJK = (2,2,2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });
    }
    else if (comm.rank() == 3) {
        // IJK = (2,2,3)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(1, 1, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    if (comm.rank() != owner) {
        return;
    }

    const auto rPos = xw.find("R");

    BOOST_REQUIRE_MESSAGE(rPos != xw.end(),
                          R"(There must be dynamic results for well "R")");

    BOOST_CHECK_EQUAL(rPos->second.connections.size(), std::size_t{3});

    // (2,2,1)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 0));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,1))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
    }

    // (2,2,2)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 1));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,2))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
    }

    // (2,2,3)
    {
        const auto* xc = rPos->second.find_connection
            (cse.es.getInputGrid().getGlobalIndex(1, 1, 2));

        BOOST_REQUIRE_MESSAGE(xc != nullptr,
                              R"(Well "R" must have dynamic results in connection (2,2,3))");

        BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
        BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Owned_And_Split)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // R1 is owned by rank 2 and intersected on rank 2.  R2 is owned by rank
    // 3 and intersected on ranks 0, 1, and 3.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 2) && (wname == "R1"))
                || ((r == 3) && (wname == "R2"));
        },
        [r = comm.rank()](const std::string& wname) {
            return (((r == 0) || (r == 1) || (r == 3)) && (wname == "R2"))
                || ((r == 2) && (wname == "R1"));
        }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });

        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });
    }

    if (comm.rank() == 3) {
        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });
    }

    if (comm.rank() == 0) {
        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 3) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Split_Wells_Even_Odd)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // R1 is owned by rank 2 and intersected on ranks 0 and 2.  R2 is owned
    // by rank 1 and intersected on ranks 1 and 3.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 2) && (wname == "R1"))
                || ((r == 1) && (wname == "R2"));
        },
        [r = comm.rank()](const std::string& wname) {
            return ((r % 2 == 0) && (wname == "R1"))
                || ((r % 2 == 1) && (wname == "R2"));
        }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    if (comm.rank() == 2) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    if (comm.rank() == 3) {
        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });

        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(Quarter_Five_Spot_Split_Wells_Cross)
{
    const auto cse = quarterFiveSpot();

    const auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    // R1 is owned by rank 2 and intersected on ranks 1 and 2.  R2 is owned
    // by rank 0 and intersected on ranks 0, 2 and 3.
    auto rftc = RFTContainer {
        cse.es, cse.sched,
        [r = comm.rank()](const std::string& wname) {
            return ((r == 2) && (wname == "R1"))
                || ((r == 0) && (wname == "R2"));
        },
        [r = comm.rank()](const std::string& wname) {
            return (((r == 1) || (r == 2)) && (wname == "R1"))
                || ((r != 1) && (wname == "R2"));
        }
    };

    rftc.allocate(0);

    // --------------------------------------------------------------------


    if (comm.rank() == 0) {
        // IJK = (3,3,3) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 2),
                    /* Po = */ []() { return 199;  },
                    /* Sw = */ []() { return 0.42; },
                    /* Sg = */ []() { return 0.07; });
    }

    if (comm.rank() == 2) {
        // IJK = (1,1,1) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 0),
                    /* Po = */ []() { return 100;  },
                    /* Sw = */ []() { return 0.01; },
                    /* Sg = */ []() { return 0.8;  });

        // IJK = (3,3,2) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 1),
                    /* Po = */ []() { return 152;  },
                    /* Sw = */ []() { return 0.12; },
                    /* Sg = */ []() { return 0.48; });

        // IJK = (1,1,2) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 1),
                    /* Po = */ []() { return 150; },
                    /* Sw = */ []() { return 0.1; },
                    /* Sg = */ []() { return 0.5; });
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 1) {
        // IJK = (1,1,3) (R1)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(0, 0, 2),
                    /* Po = */ []() { return 215; },
                    /* Sw = */ []() { return 0.4; },
                    /* Sg = */ []() { return 0.1; });
    }

    if (comm.rank() == 3) {
        // IJK = (3,3,1) (R2)
        rftc.assign(cse.es.getInputGrid().getGlobalIndex(2, 2, 0),
                    /* Po = */ []() { return 123;  },
                    /* Sw = */ []() { return 0.02; },
                    /* Sg = */ []() { return 0.63; });
    }

    // --------------------------------------------------------------------

    auto xw = Opm::data::Wells{};
    rftc.addToWells(xw, 0, comm);

    // --------------------------------------------------------------------

    if (comm.rank() == 2) {
        const auto r1Pos = xw.find("R1");

        BOOST_REQUIRE_MESSAGE(r1Pos != xw.end(),
                              R"(There must be dynamic results for well "R1")");

        BOOST_CHECK_EQUAL(r1Pos->second.connections.size(), std::size_t{3});

        // (1,1,1)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 100.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.01, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.8 , 1.0e-8);
        }

        // (1,1,2)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 150.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.1, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.5, 1.0e-8);
        }

        // (1,1,3)
        {
            const auto* xc = r1Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(0, 0, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R1" must have dynamic results in connection (1,1,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 215.0, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.4, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.1, 1.0e-8);
        }
    }

    // --------------------------------------------------------------------

    if (comm.rank() == 0) {
        const auto r2Pos = xw.find("R2");

        BOOST_REQUIRE_MESSAGE(r2Pos != xw.end(),
                              R"(There must be dynamic results for well "R2")");

        BOOST_CHECK_EQUAL(r2Pos->second.connections.size(), std::size_t{3});

        // (3,3,1)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 0));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,1))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 123.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.02, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.63, 1.0e-8);
        }

        // (3,3,2)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 1));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,2))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 152.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.12, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.48, 1.0e-8);
        }

        // (3,3,3)
        {
            const auto* xc = r2Pos->second.find_connection
                (cse.es.getInputGrid().getGlobalIndex(2, 2, 2));

            BOOST_REQUIRE_MESSAGE(xc != nullptr,
                                  R"(Well "R2" must have dynamic results in connection (3,3,3))");

            BOOST_CHECK_CLOSE(xc->cell_pressure        , 199.0 , 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_water,   0.42, 1.0e-8);
            BOOST_CHECK_CLOSE(xc->cell_saturation_gas  ,   0.07, 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

// ===========================================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    register_error_handler();

    // Just enough initialisation for a working TestFluidSystem::phaseIsActive().
    //
    // This is a **massive** hack.
    TestFluidSystem::initBegin(1);

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
