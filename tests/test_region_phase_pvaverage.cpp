// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor.

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

#define BOOST_TEST_MODULE Region_Phase_PVAverage

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/RegionPhasePVAverage.hpp>

#include <opm/simulators/flow/partitionCells.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

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

    template <typename GlobalVector>
    auto selectSubset(const std::vector<bool>& keep,
                      const GlobalVector&      global)
    {
        using VT = std::remove_cv_t<std::remove_reference_t<
            typename std::iterator_traits<decltype(global.begin())>::value_type>>;

        auto local = std::vector<VT>{};

        const auto n = keep.size();
        auto elm = global.begin();
        for (auto i = 0*n; i < n; ++i, ++elm) {
            if (keep[i]) {
                local.push_back(*elm);
            }
        }

        return local;
    }

    std::vector<bool> cellSubset(const int rank, const std::vector<int>& p)
    {
        auto keep = std::vector<bool>(p.size(), false);

        const auto n = p.size();
        for (auto i = 0*n; i < n; ++i) {
            keep[i] = p[i] == rank;
        }

        return keep;
    }

    std::vector<int> fipnum()
    {
        return {
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        };
    }

    std::vector<int> fiplayer()
    {
        return {
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            1, 1, 1,
            1, 1, 1,
            1, 1, 1,

            2, 2, 2,
            2, 2, 2,
            2, 2, 2,
        };
    }

    class RegionSets
    {
    public:
        RegionSets();

        RegionSets& selectSubset(const std::vector<bool>& keep);

        std::vector<std::string> names() const;
        const std::vector<int>& operator()(const std::string& rname) const;

    private:
        std::unordered_map<std::string, std::vector<int>> regions_;
    };

    RegionSets::RegionSets()
        : regions_ {
                { "FIPNUM", fipnum() },
                { "FIPLRS", fiplayer() },
            }
    {}

    RegionSets& RegionSets::selectSubset(const std::vector<bool>& keep)
    {
        for (auto& [name, region] : this->regions_) {
            if (region.size() != keep.size()) {
                continue;
            }

            auto local = ::selectSubset(keep, region);
            region.swap(local);
        }

        return *this;
    }

    std::vector<std::string> RegionSets::names() const
    {
        auto rnames = std::vector<std::string> {};
        rnames.reserve(this->regions_.size());

        for (const auto& region : this->regions_) {
            rnames.push_back(region.first);
        }

        return rnames;
    }

    const std::vector<int>& RegionSets::operator()(const std::string& rname) const
    {
        auto rpos = this->regions_.find(rname);
        if (rpos == this->regions_.end()) {
            throw std::invalid_argument {
                fmt::format("Unknown region set '{}'", rname)
            };
        }

        return rpos->second;
    }

    bool init_unit_test_func()
    {
        return true;
    }

} // Anonymous namespace

// ===========================================================================

using PVAvg = Opm::RegionPhasePoreVolAverage;

BOOST_AUTO_TEST_SUITE(Sequential)

BOOST_AUTO_TEST_CASE(Single_Phase_Single_Val_Per_Layer)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto rset = RegionSets{};
    const auto numPhases = std::size_t{1};
    const auto p = PVAvg::Phase {0};

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    BOOST_CHECK_MESSAGE(std::isnan(avgCalc.fieldValue(p)),
                        "Field value must be NaN prior to accumulation");

    avgCalc.prepareAccumulation();

    {
        const auto activeCell = 0*9 + std::size_t{0};
        const auto cv = PVAvg::CellValue {
            1.0, 0.25, 200.0
        };

        avgCalc.addCell(activeCell, p, cv);
    }

    {
        const auto activeCell = 1*9 + std::size_t{0};
        const auto cv = PVAvg::CellValue {
            2.0, 0.0, 200.0
        };

        avgCalc.addCell(activeCell, p, cv);
    }

    {
        const auto activeCell = 2*9 + std::size_t{0};
        const auto cv = PVAvg::CellValue {
            3.0, 0.5, 200.0
        };

        avgCalc.addCell(activeCell, p, cv);
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(p), 7.0/3.0, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", p, PVAvg::Region{0}), 7.0/3.0, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{0}), 1.0, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{1}), 2.0, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{2}), 3.0, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Single_Phase_Full)
{
    const auto x = std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
    };

    const auto s = std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto rset = RegionSets{};
    const auto numPhases = std::size_t{1};
    const auto p = PVAvg::Phase {0};

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
        const auto cv = PVAvg::CellValue {
            x[c], s[c], 200.0
        };

        avgCalc.addCell(c, p, cv);
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(p), 1.61, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", p, PVAvg::Region{0}), 1.61, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{0}), 1.11, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{1}), 2.11, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{2}), 3.11, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Single_Phase_Full_Varying_PV)
{
    const auto x = std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
    };

    const auto s = std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    };

    const auto pv = std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
    };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto rset = RegionSets{};
    const auto numPhases = std::size_t{1};
    const auto p = PVAvg::Phase {0};

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
        const auto cv = PVAvg::CellValue {
            x[c], s[c], pv[c]
        };

        avgCalc.addCell(c, p, cv);
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(p), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", p, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Three_Phase_All_Present_Full_Varying_PV)
{
    const auto x = std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
    };

    const auto s = std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    };

    const auto pv = std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
    };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto transform = std::array {
        std::function { [](const double s0) { return s0 / 2; } },
        std::function { [](const double s0) { return s0 / 4; } },
        std::function { [](const double s0) { return 1 - (3*s0/4); } },
    };

    const auto rset = RegionSets{};
    const auto numPhases = transform.size();

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    {
        auto p = PVAvg::Phase {0};

        for (const auto& t : transform) {
            for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
                const auto cv = PVAvg::CellValue {
                    x[c], t(s[c]), pv[c]
                };

                avgCalc.addCell(c, p, cv);
            }

            ++p.ix;
        }
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{1}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{2}), 1.471079741628658, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{0}, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{1}, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{2}, PVAvg::Region{0}), 1.471079741628658, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{0}), 1.111326195193250, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{1}), 2.156805281711179, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Three_Phase_Two_Present_Full_Varying_PV)
{
    const auto x = std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
    };

    const auto s = std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.25, 0.30, 0.25,
        0.30, 0.35, 0.30,
        0.25, 0.30, 0.25,
    };

    const auto pv = std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
    };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto transform = std::array {
        std::function { [](const double s0) { return s0; } },
        std::function { [](const double s0) { return 0 * s0; } },
        std::function { [](const double s0) { return 1 - s0; } },
    };

    const auto rset = RegionSets{};
    const auto numPhases = transform.size();

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    {
        auto p = PVAvg::Phase {0};

        for (const auto& t : transform) {
            for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
                const auto cv = PVAvg::CellValue {
                    x[c], t(s[c]), pv[c]
                };

                avgCalc.addCell(c, p, cv);
            }

            ++p.ix;
        }
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{0}), 2.203870000146282, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{1}), 1.460875637211809, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{2}), 1.388846263879987, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{0}, PVAvg::Region{0}), 2.203870000146282, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{1}, PVAvg::Region{0}), 1.460875637211809, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{2}, PVAvg::Region{0}), 1.388846263879987, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{2}), 3.081133011812399, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{0}), 1.111895506705607, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{1}), 2.156118568232662, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{0}), 1.111126811726796, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{1}), 2.157055514795793, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{2}), 3.079851244928804, 1.0e-8);
}

BOOST_AUTO_TEST_SUITE_END()     // Sequential

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Parallel)

BOOST_AUTO_TEST_CASE(Single_Phase_Full)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto keep = cellSubset(comm.rank(), Opm::partitionCellsSimple(3 * 3 * 3, comm.size()).first);

    const auto x = selectSubset(keep, std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
        });

    const auto s = selectSubset(keep, std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        });

    const auto rset = RegionSets{}.selectSubset(keep);
    const auto numPhases = std::size_t{1};
    const auto p = PVAvg::Phase {0};

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
        const auto cv = PVAvg::CellValue {
            x[c], s[c], 200.0
        };

        avgCalc.addCell(c, p, cv);
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(p), 1.61, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", p, PVAvg::Region{0}), 1.61, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{0}), 1.11, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{1}), 2.11, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{2}), 3.11, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Single_Phase_Full_Varying_PV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto keep = cellSubset(comm.rank(), Opm::partitionCellsSimple(3 * 3 * 3, comm.size()).first);

    const auto x = selectSubset(keep, std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
        });

    const auto s = selectSubset(keep, std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        });

    const auto pv = selectSubset(keep, std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
        });

    const auto rset = RegionSets{}.selectSubset(keep);
    const auto numPhases = std::size_t{1};
    const auto p = PVAvg::Phase {0};

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
        const auto cv = PVAvg::CellValue {
            x[c], s[c], pv[c]
        };

        avgCalc.addCell(c, p, cv);
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(p), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", p, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", p, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Three_Phase_All_Present_Full_Varying_PV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto keep = cellSubset(comm.rank(), Opm::partitionCellsSimple(3 * 3 * 3, comm.size()).first);

    const auto x = selectSubset(keep, std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
        });

    const auto s = selectSubset(keep, std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.0, 0.0, 0.0,          // s=0 in layer 3 => use PV average instead
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        });

    const auto pv = selectSubset(keep, std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
        });

    const auto transform = std::array {
        std::function { [](const double s0) { return s0 / 2; } },
        std::function { [](const double s0) { return s0 / 4; } },
        std::function { [](const double s0) { return 1 - (3*s0/4); } },
    };

    const auto rset = RegionSets{}.selectSubset(keep);
    const auto numPhases = transform.size();

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    {
        auto p = PVAvg::Phase {0};

        for (const auto& t : transform) {
            for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
                const auto cv = PVAvg::CellValue {
                    x[c], t(s[c]), pv[c]
                };

                avgCalc.addCell(c, p, cv);
            }

            ++p.ix;
        }
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{1}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{2}), 1.471079741628658, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{0}, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{1}, PVAvg::Region{0}), 1.128401081038469, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{2}, PVAvg::Region{0}), 1.471079741628658, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{0}), 1.111326195193250, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{1}), 2.156805281711179, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Three_Phase_Two_Present_Full_Varying_PV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto keep = cellSubset(comm.rank(), Opm::partitionCellsSimple(3 * 3 * 3, comm.size()).first);

    const auto x = selectSubset(keep, std::array {
        // K=1
        1.0, 1.01, 1.02,
        1.1, 1.11, 1.12,
        1.2, 1.21, 1.22,

        // K=2
        2.0, 2.01, 2.02,
        2.1, 2.11, 2.12,
        2.2, 2.21, 2.22,

        // K=3
        3.0, 3.01, 3.02,
        3.1, 3.11, 3.12,
        3.2, 3.21, 3.22,
        });

    const auto s = selectSubset(keep, std::array {
        // K=1
        0.0, 0.1, 0.0,
        0.1, 0.4, 0.1,
        0.0, 0.1, 0.0,

        // K=2
        0.1, 0.0, 0.1,
        0.0, 0.4, 0.0,
        0.1, 0.0, 0.1,

        // K=3
        0.25, 0.30, 0.25,
        0.30, 0.35, 0.30,
        0.25, 0.30, 0.25,
        });

    const auto pv = selectSubset(keep, std::array {
        // K=1
        1729.0,  271.82, 3141.5 ,
           0.1, 1000.0,   321.09,
        4321.0, 1234.5,    67.89,

        // K=2
        0.12, 0.34, 0.56,
        0.78, 0.90, 1.23,
        1.45, 1.67, 1.89,

        // K=3
        500.0, 450.0, 400.0,
        350.0,   1.0, 300.0,
        250.0, 200.0, 150.0,
        });

    const auto transform = std::array {
        std::function { [](const double s0) { return s0; } },
        std::function { [](const double s0) { return 0 * s0; } },
        std::function { [](const double s0) { return 1 - s0; } },
    };

    const auto rset = RegionSets{}.selectSubset(keep);
    const auto numPhases = transform.size();

    auto avgCalc = PVAvg {
        comm, numPhases, rset.names(), rset
    };

    avgCalc.prepareAccumulation();

    {
        auto p = PVAvg::Phase {0};

        for (const auto& t : transform) {
            for (auto nc = x.size(), c = 0*nc; c < nc; ++c) {
                const auto cv = PVAvg::CellValue {
                    x[c], t(s[c]), pv[c]
                };

                avgCalc.addCell(c, p, cv);
            }

            ++p.ix;
        }
    }

    avgCalc.accumulateParallel();

    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{0}), 2.203870000146282, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{1}), 1.460875637211809, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.fieldValue(PVAvg::Phase{2}), 1.388846263879987, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{0}, PVAvg::Region{0}), 2.203870000146282, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{1}, PVAvg::Region{0}), 1.460875637211809, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPNUM", PVAvg::Phase{2}, PVAvg::Region{0}), 1.388846263879987, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{0}), 1.127070395417597, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{1}), 2.146062992125984, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{0}, PVAvg::Region{2}), 3.081133011812399, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{0}), 1.111895506705607, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{1}), 2.156118568232662, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{1}, PVAvg::Region{2}), 3.080203767781622, 1.0e-8);

    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{0}), 1.111126811726796, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{1}), 2.157055514795793, 1.0e-8);
    BOOST_CHECK_CLOSE(avgCalc.value("FIPLRS", PVAvg::Phase{2}, PVAvg::Region{2}), 3.079851244928804, 1.0e-8);
}

BOOST_AUTO_TEST_SUITE_END()     // Parallel

// ===========================================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    register_error_handler();

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
