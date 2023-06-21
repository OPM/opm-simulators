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

#define BOOST_TEST_MODULE Parallel_PAvg_Dynamic_Source_Data

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/input/eclipse/Schedule/Well/PAvgDynamicSourceData.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <cstddef>
#include <functional>
#include <iostream>
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
    std::vector<char> err_string_vec(MPI_MAX_ERROR_STRING);
    auto err_length = 0;

    MPI_Error_string(*err_code, err_string_vec.data(), &err_length);

    auto err_string = std::string_view {
        err_string_vec.data(), static_cast<std::string_view::size_type>(err_length)
    };

    std::cerr << "An MPI Error ocurred:\n  -> " << err_string << '\n';

    throw MPIError { err_string, *err_code };
}
#endif // HAVE_MPI

bool init_unit_test_func()
{
    return true;
}

std::vector<std::size_t> sourceLocations(const std::size_t numLoc)
{
    auto srcLoc = std::vector<std::size_t>(numLoc);
    std::iota(srcLoc.begin(), srcLoc.end(), std::size_t{0});
    return srcLoc;
}

class LocalCellIndex
{
public:
    explicit LocalCellIndex(const std::size_t rank,
                            const std::size_t size)
        : rank_ { rank }
        , size_ { size }
    {}

    int operator()(const std::size_t i) const
    {
        return ((i % this->size_) == this->rank_)
            ? static_cast<int>(i / this->size_)
            : -1;
    }

private:
    std::size_t rank_{};
    std::size_t size_{};
};

class CalculateSourceTerm
{
public:
    using SrcTerm = Opm::PAvgDynamicSourceData::SourceDataSpan<double>;

    explicit CalculateSourceTerm(const std::size_t rank)
        : rank_ { rank }
    {}

    void operator()(const int i, SrcTerm source_term) const
    {
        using Item = typename SrcTerm::Item;

        source_term
            .set(Item::Pressure      , this->rank_*314.15 -      i)
            .set(Item::PoreVol       , this->rank_*172.9  + 10.0*i)
            .set(Item::MixtureDensity, this->rank_*852.96 +      i);
    }

private:
    std::size_t rank_{};
};

std::size_t
sourceTermsAreCorrect(const std::size_t                 comm_size,
                      const std::size_t                 num_src,
                      const Opm::PAvgDynamicSourceData& source_data)
{
    using Item = Opm::PAvgDynamicSourceData::SourceDataSpan<const double>::Item;

    auto num_correct = 0*num_src;

    for (auto srcID = 0*num_src; srcID < num_src; ++srcID) {
        const auto rank = srcID % comm_size;
        const auto locI = srcID / comm_size;

        const auto src = source_data[srcID];
        const auto ok =
            (src[Item::Pressure]       == rank*314.15 -    locI) &&
            (src[Item::PoreVol]        == rank*172.9  + 10*locI) &&
            (src[Item::MixtureDensity] == rank*852.96 +    locI);

        num_correct += ok;
    }

    return num_correct == num_src;
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Eval_and_collect)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    const auto comm_rank = static_cast<std::size_t>(comm.rank());
    const auto comm_size = static_cast<std::size_t>(comm.size());

    const auto num_src_loc = std::size_t{50};

    auto source_data = Opm::ParallelPAvgDynamicSourceData {
        comm, sourceLocations(num_src_loc),
        LocalCellIndex { comm_rank, comm_size }
    };

    source_data.collectLocalSources(CalculateSourceTerm { comm_rank });
    source_data.synchroniseSources();

    const auto num_rank_correct = comm.sum
        (sourceTermsAreCorrect(comm_size, num_src_loc, source_data));

    if (comm_rank == 0) {
        BOOST_CHECK_EQUAL(num_rank_correct, comm_size);
    }
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
