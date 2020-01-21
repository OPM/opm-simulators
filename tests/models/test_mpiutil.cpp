/*
  Copyright 2020 Equinor ASA.

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

#include <opm/models/parallel/mpiutil.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <cassert>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

bool noStrings(int, int)
{
    std::string empty;
    auto res = Opm::gatherStrings(empty);
    assert(res.empty());
    return true;
}

bool oddRankStrings(int size, int rank)
{
    std::string what = (rank % 2 == 1) ? "An error on rank " + std::to_string(rank) : std::string();
    auto res = Opm::gatherStrings(what);
    assert(int(res.size()) == size/2);
    for (int i = 0; i < size/2; ++i) {
        assert(res[i] == "An error on rank " + std::to_string(2*i + 1));
    }
    return true;
}

bool allRankStrings(int size, int rank)
{
    std::string what = "An error on rank " + std::to_string(rank);
    auto res = Opm::gatherStrings(what);
    assert(int(res.size()) == size);
    for (int i = 0; i < size; ++i) {
        assert(res[i] == "An error on rank " + std::to_string(i));
    }
    return true;
}


int testMain(int size, int rank)
{
    bool ok = noStrings(size, rank);
    ok = ok && oddRankStrings(size, rank);
    ok = ok && allRankStrings(size, rank);
    if (ok) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}


int main(int argc, char** argv)
{
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
    int mpiSize = mpiHelper.size();
    int mpiRank = mpiHelper.rank();
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return testMain(mpiSize, mpiRank);
}
