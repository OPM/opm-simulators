/*
  Copyright 2023 Equinor.

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

#define BOOST_TEST_MODULE Parallel_WBPn_Calculation

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/ParallelWBPCalculation.hpp>

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvg.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgDynamicSourceData.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

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
#endif // HAVE_MPI

    double standardGravity()
    {
        return Opm::unit::gravity;
    }

    std::size_t globIndex(const std::array<int,3>& ijk,
                          const std::array<int,3>& dims)
    {
        return ijk[0] + dims[0]*(ijk[1] + static_cast<std::size_t>(dims[1])*ijk[2]);
    }

    std::array<int,3> cellIJK(int                      cell,
                              const std::array<int,3>& dims)
    {
        auto ijk = std::array<int,3>{};

        ijk[0] = cell % dims[0];  cell /= dims[0];
        ijk[1] = cell % dims[1];
        ijk[2] = cell / dims[1];

        return ijk;
    }

    namespace Rank {
        namespace Top {
            int globalToLocal(const std::size_t global)
            {
                return (global >= 5 * 5 * 5)
                    ? -1
                    : static_cast<int>(global);
            }

            bool isInRange(const std::array<int,3>& ijk)
            {
                // Well block column in top half covers zero-based index
                // range
                //
                //   [ 1..3, 1..3, 2..4 ]
                //
                // of index range
                //
                //   [ 0..4, 0..4 , 0..4 ]

                return (ijk[0] >= 1) && (ijk[0] <= 3)
                    && (ijk[1] >= 1) && (ijk[1] <= 3)
                    && (ijk[2] >= 2) && (ijk[2] <= 4);
            }

            std::size_t fieldIx(const std::array<int,3>& ijk)
            {
                return globIndex({ ijk[0] - 1, ijk[1] - 1, ijk[2] - 2 }, {3, 3, 3});
            }

            double fieldValue(const int cell, std::initializer_list<double> field)
            {
                const auto ijk = cellIJK(cell, { 5, 5, 5 });

                return isInRange(ijk) ? *(std::data(field) + fieldIx(ijk)) : 0.0;
            }

            // Octave: 1234 + fix(100 * rand([3, 3, 6])) -- top half
            double pressure(const int cell)
            {
                return fieldValue(cell, {
                        // K=2
                        1.302000e+03, 1.308000e+03, 1.279000e+03,
                        1.242000e+03, 1.256000e+03, 1.325000e+03,
                        1.249000e+03, 1.316000e+03, 1.287000e+03,

                        // K=3
                        1.333000e+03, 1.241000e+03, 1.278000e+03,
                        1.244000e+03, 1.330000e+03, 1.234000e+03,
                        1.311000e+03, 1.315000e+03, 1.320000e+03,

                        // K=4
                        1.242000e+03, 1.273000e+03, 1.259000e+03,
                        1.314000e+03, 1.277000e+03, 1.325000e+03,
                        1.252000e+03, 1.260000e+03, 1.248000e+03,
                    });
            }

            // Octave: fix(1e6 * (123.4 + 56.7*rand([3, 3, 6]))) / 1e6 -- top half
            double porevol(const int cell)
            {
                return fieldValue(cell, {
                        // K=2
                        1.301471680e+02, 1.516572410e+02, 1.778174820e+02,
                        1.426998700e+02, 1.565846810e+02, 1.360901360e+02,
                        1.659968420e+02, 1.378638930e+02, 1.520877640e+02,

                        // K=3
                        1.630376500e+02, 1.739142140e+02, 1.777918230e+02,
                        1.544271200e+02, 1.312600050e+02, 1.318649700e+02,
                        1.380007180e+02, 1.710686680e+02, 1.378177990e+02,

                        // K=4
                        1.695699490e+02, 1.372078650e+02, 1.760892470e+02,
                        1.432440790e+02, 1.345469500e+02, 1.376364540e+02,
                        1.583297330e+02, 1.502354770e+02, 1.433390940e+02,
                    });
            }

            // Octave: 0.1 + round(0.1 * rand([3, 3, 6]), 2) -- top half
            double density(const int cell)
            {
                return fieldValue(cell, {
                        // K=2
                        0.120, 0.120, 0.120,
                        0.140, 0.130, 0.190,
                        0.140, 0.120, 0.190,

                        // K=3
                        0.200, 0.140, 0.110,
                        0.130, 0.140, 0.160,
                        0.130, 0.160, 0.170,

                        // K=4
                        0.120, 0.110, 0.130,
                        0.130, 0.140, 0.150,
                        0.110, 0.130, 0.180,
                    });
            }

            void cellSource(const int                                          cell,
                            Opm::PAvgDynamicSourceData::SourceDataSpan<double> src)
            {
                using Item = Opm::PAvgDynamicSourceData::SourceDataSpan<double>::Item;

                src .set(Item::Pressure      , pressure(cell))
                    .set(Item::PoreVol       , porevol (cell))
                    .set(Item::MixtureDensity, density (cell));
            }

            std::vector<int> localConnIdx()
            {
                auto localIdx = std::vector<int>(6, -1);
                for (auto perf = 0; perf < 3; ++perf) {
                    localIdx[perf] = perf;
                }

                return localIdx;
            }

            Opm::ParallelWBPCalculation::EvaluatorFactory connSource()
            {
                return []() {
                    auto rho = std::vector { 0.1, 0.12, 0.14, };

                    return [rho = std::move(rho)]
                        (const int                                          connIx,
                         Opm::PAvgDynamicSourceData::SourceDataSpan<double> src)
                    {
                        using Item = Opm::PAvgDynamicSourceData::SourceDataSpan<double>::Item;

                        src .set(Item::Pressure      , 1222.0)
                            .set(Item::PoreVol       ,    1.25)
                            .set(Item::MixtureDensity, rho[connIx]);
                    };
                };
            }
        } // namespace Top

        namespace Bottom {
            int globalToLocal(const std::size_t global)
            {
                constexpr auto middle = 5 * 5 * 5;

                return (global < middle)
                    ? -1
                    : static_cast<int>(global - middle);
            }

            bool isInRange(const std::array<int,3>& ijk)
            {
                // Well block column in bottom half covers zero-based index
                // range
                //
                //   [ 1..3, 1..3, 0..2 ]
                //
                // of index range
                //
                //   [ 0..4, 0..4 , 0..4 ]

                return (ijk[0] >= 1) && (ijk[0] <= 3)
                    && (ijk[1] >= 1) && (ijk[1] <= 3)
                    && (ijk[2] >= 0) && (ijk[2] <= 2);
            }

            std::size_t fieldIx(const std::array<int,3>& ijk)
            {
                return globIndex({ ijk[0] - 1, ijk[1] - 1, ijk[2] }, {3, 3, 3});
            }

            double fieldValue(const int cell, std::initializer_list<double> field)
            {
                const auto ijk = cellIJK(cell, { 5, 5, 5 });

                return isInRange(ijk) ? *(std::data(field) + fieldIx(ijk)) : 0.0;
            }

            // Octave: 1234 + fix(100 * rand([3, 3, 6])) -- bottom half
            double pressure(const int cell)
            {
                return fieldValue(cell, {
                    // K=5
                    1.247000e+03, 1.320000e+03, 1.291000e+03,
                    1.288000e+03, 1.248000e+03, 1.319000e+03,
                    1.296000e+03, 1.269000e+03, 1.285000e+03,

                    // K=6
                    1.274000e+03, 1.241000e+03, 1.257000e+03,
                    1.246000e+03, 1.252000e+03, 1.257000e+03,
                    1.275000e+03, 1.238000e+03, 1.324000e+03,

                    // K=7
                    1.328000e+03, 1.283000e+03, 1.282000e+03,
                    1.267000e+03, 1.324000e+03, 1.270000e+03,
                    1.245000e+03, 1.312000e+03, 1.272000e+03,
                    });
            }

            // Octave: fix(1e6 * (123.4 + 56.7*rand([3, 3, 6]))) / 1e6 -- bottom half
            double porevol(const int cell)
            {
                return fieldValue(cell, {
                        // K=5
                        1.705079830e+02, 1.565844730e+02, 1.545693280e+02,
                        1.754048800e+02, 1.396070720e+02, 1.663332520e+02,
                        1.661364390e+02, 1.449712790e+02, 1.555954870e+02,

                        // K=6
                        1.277009380e+02, 1.264589710e+02, 1.534962210e+02,
                        1.675787810e+02, 1.763584050e+02, 1.307656820e+02,
                        1.556523010e+02, 1.500144490e+02, 1.240748470e+02,

                        // K=7
                        1.425148530e+02, 1.325957360e+02, 1.684359330e+02,
                        1.410458920e+02, 1.533678280e+02, 1.327922820e+02,
                        1.575323760e+02, 1.383104710e+02, 1.604862840e+02,
                    });
            }

            // Octave: 0.1 + round(0.1 * rand([3, 3, 6]), 2) -- bottom half
            double density(const int cell)
            {
                return fieldValue(cell, {
                        // K=5
                        0.100, 0.190, 0.170,
                        0.150, 0.160, 0.120,
                        0.150, 0.200, 0.150,

                        // K=6
                        0.150, 0.120, 0.150,
                        0.160, 0.170, 0.140,
                        0.140, 0.200, 0.100,

                        // K=7
                        0.190, 0.190, 0.180,
                        0.110, 0.130, 0.130,
                        0.170, 0.110, 0.170,
                    });
            }

            void cellSource(const int                                          cell,
                            Opm::PAvgDynamicSourceData::SourceDataSpan<double> src)
            {
                using Item = Opm::PAvgDynamicSourceData::SourceDataSpan<double>::Item;

                src .set(Item::Pressure      , pressure(cell))
                    .set(Item::PoreVol       , porevol (cell))
                    .set(Item::MixtureDensity, density (cell));
            }

            std::vector<int> localConnIdx()
            {
                auto localIdx = std::vector<int>(6, -1);
                for (auto perf = 0; perf < 3; ++perf) {
                    localIdx[3 + perf] = perf;
                }

                return localIdx;
            }

            Opm::ParallelWBPCalculation::EvaluatorFactory connSource()
            {
                return []() {
                    auto rho = std::vector { 0.16, 0.18, 0.2, };

                    return [rho = std::move(rho)]
                        (const int                                          connIx,
                         Opm::PAvgDynamicSourceData::SourceDataSpan<double> src)
                    {
                        using Item = Opm::PAvgDynamicSourceData::SourceDataSpan<double>::Item;

                        src .set(Item::Pressure      , 1222.0)
                            .set(Item::PoreVol       ,    1.25)
                            .set(Item::MixtureDensity, rho[connIx]);
                    };
                };
            }
        } // namespace Bottom
    } // namespace Rank

    std::shared_ptr<Opm::WellConnections>
    centreConnections(const int topConn, const int numConns)
    {
        auto conns = std::vector<Opm::Connection>{};

        const auto dims = std::array { 5, 5, 10 };

        const auto i = 2;
        const auto j = 2;
        const auto kMax = std::min(dims[2] - 1, topConn + numConns);

        const auto state = std::array {
            Opm::Connection::State::OPEN,
            Opm::Connection::State::SHUT,
            Opm::Connection::State::OPEN,
        };

        for (auto k = topConn; k < kMax; ++k) {
            const auto depth = 2000 + (2*k + 1) / static_cast<double>(2);

            auto ctf_props = Opm::Connection::CTFProperties{};

            // 0.03, 0.0, 0.01, 0.02, 0.03, ...
            ctf_props.CF = ((k + 3 - topConn) % 4) / 100.0;

            ctf_props.Kh = 1.0;
            ctf_props.Ke = 1.0;
            ctf_props.rw = 1.0;
            ctf_props.r0 = 0.5;
            ctf_props.re = 0.5;
            ctf_props.connection_length = 1.0;

            conns.emplace_back(i, j, k, globIndex({i, j, k}, dims), k - topConn,

                               // Open, Shut, Open, Open, Shut, ...
                               state[(k - topConn) % state.size()],

                               Opm::Connection::Direction::Z,
                               Opm::Connection::CTFKind::DeckValue,
                               0, depth, ctf_props, k - topConn, false);
        }

        return std::make_shared<Opm::WellConnections>
            (Opm::Connection::Order::INPUT, i, j, conns);
    }

    Opm::Well producerWell()
    {
        auto w = Opm::Well {
            "P", "G", 0, 0, 2, 2, 2000.5,
            Opm::WellType { true, Opm::Phase::OIL }, // Oil producer
            Opm::Well::ProducerCMode::ORAT,
            Opm::Connection::Order::INPUT,
            Opm::UnitSystem::newMETRIC(),
            -3.0e+20,           // UDQ undefined
            0.0, true, true, 0,
            Opm::Well::GasInflowEquation::STD
        };

        w.updateConnections(centreConnections(2, 6), true);

        return w;
    }

    Opm::ParallelWellInfo
    parallelWellInfo(const Opm::Parallel::Communication& comm)
    {
        auto pwi = Opm::ParallelWellInfo {
            std::pair { std::string{ "P" }, true }, comm
        };

        pwi.beginReset();

        const auto numLocalPerf = 3;
        const auto perfOffset = comm.rank() * numLocalPerf;

        auto prev = Opm::ParallelWellInfo::INVALID_ECL_INDEX;
        for (auto perf = 0; perf < numLocalPerf; ++perf) {
            const auto curr = perfOffset + perf;
            pwi.pushBackEclIndex(prev, curr);
            prev = curr;
        }

        pwi.endReset();

        pwi.communicateFirstPerforation(comm.rank() == 0);

        return pwi;
    }

    void setCallbacksTop(Opm::ParallelWBPCalculation& wbpCalcService)
    {
        wbpCalcService
            .localCellIndex(&Rank::Top::globalToLocal)
            .evalCellSource(&Rank::Top::cellSource);
    }

    void setCallbacksBottom(Opm::ParallelWBPCalculation& wbpCalcService)
    {
        wbpCalcService
            .localCellIndex(&Rank::Bottom::globalToLocal)
            .evalCellSource(&Rank::Bottom::cellSource);
    }

    void setCallbacks(const int                    rank,
                      Opm::ParallelWBPCalculation& wbpCalcService)
    {
        if (rank == 0) {
            setCallbacksTop(wbpCalcService);
        }
        else {
            setCallbacksBottom(wbpCalcService);
        }
    }

    Opm::ParallelWBPCalculation::EvaluatorFactory connSource(const int rank)
    {
        if (rank == 0) {
            return Rank::Top::connSource();
        }
        else {
            return Rank::Bottom::connSource();
        }
    }

    std::vector<int> localConnIdx(const int rank)
    {
        if (rank == 0) {
            return Rank::Top::localConnIdx();
        }
        else {
            return Rank::Bottom::localConnIdx();
        }
    }

    bool init_unit_test_func()
    {
        return true;
    }

    struct Setup
    {
        Setup()
            : comm           { Dune::MPIHelper::getCommunicator() }
            , cellIndexMap   { 5, 5, 10 }
            , wbpCalcService { cellIndexMap, comm }
            , pwi            { parallelWellInfo(comm) }
        {
            setCallbacks(this->comm.rank(), this->wbpCalcService);

            this->wbpCalcService
                .createCalculator(producerWell(), this->pwi,
                                  localConnIdx(this->comm.rank()),
                                  connSource(this->comm.rank()));
        }

        Opm::Parallel::Communication comm;
        Opm::GridDims cellIndexMap;
        Opm::ParallelWBPCalculation wbpCalcService;
        Opm::ParallelWellInfo pwi;
    };

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Create)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto wbpCalcService = Opm::ParallelWBPCalculation {
        Opm::GridDims { 5, 5, 10 },
        comm
    };

    setCallbacks(comm.rank(), wbpCalcService);
    const auto pwi = parallelWellInfo(comm);

    const auto calcIdx = wbpCalcService
        .createCalculator(producerWell(), pwi,
                          localConnIdx(comm.rank()),
                          connSource(comm.rank()));

    BOOST_CHECK_EQUAL(calcIdx, std::size_t{0});
}

BOOST_AUTO_TEST_CASE(TopOfFormation_Well_OpenConns)
{
    // Producer connected in Z direction in cells (3,3,3), (3,3,4), (3,3,5),
    // (3,3,6), (3,3,7), and (3,3,8).  Connections (3,3,4) and (3,3,7) are
    // shut.

    Setup cse {};

    cse.wbpCalcService.defineCommunication();
    cse.wbpCalcService.collectDynamicValues();

    const auto calcIndex = std::size_t{0};
    const auto controls  = Opm::PAvg{};
    const auto gravity   = standardGravity();
    const auto refDepth  = 2002.5; // BHP reference depth.  Depth correction in layers 4..8.
    cse.wbpCalcService.inferBlockAveragePressures(calcIndex, controls, gravity, refDepth);

    const auto avgPress = cse.wbpCalcService.averagePressures(calcIndex);
    using WBPMode = Opm::PAvgCalculator::Result::WBPMode;

    BOOST_CHECK_CLOSE(avgPress.value(WBPMode::WBP) , 1254.806625666667, 1.0e-8);
    BOOST_CHECK_CLOSE(avgPress.value(WBPMode::WBP4), 1295.348292333333, 1.0e-8);
    BOOST_CHECK_CLOSE(avgPress.value(WBPMode::WBP5), 1275.077459000000, 1.0e-8);
    BOOST_CHECK_CLOSE(avgPress.value(WBPMode::WBP9), 1269.379542333333, 1.0e-8);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_MPI
    // Register a throwing error handler to allow for debugging with
    //
    //   catch throw
    //
    // in GDB.
    MPI_Errhandler handler{};
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif // HAVE_MPI

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
