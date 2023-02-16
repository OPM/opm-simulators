/*
  Copyright 2021 Equinor.

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

#include <opm/common/utility/FileSystem.hpp>

#include <opm/simulators/utils/HDF5File.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#define BOOST_TEST_MODULE HDF5FileParallelTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>

using namespace Opm;

BOOST_AUTO_TEST_CASE(ReadWrite)
{
    std::string path;
    Parallel::Communication comm;
    if (comm.rank() == 0) {
        path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    }
    std::size_t size = path.size();
    comm.broadcast(&size, 1, 0);
    if (comm.rank() != 0) {
        path.resize(size);
    }
    comm.broadcast(path.data(), size, 0);
    std::filesystem::create_directory(path);
    auto rwpath = (std::filesystem::path(path) / "rw.hdf5").string();
    std::vector<char> test_data(10);
    std::iota(test_data.begin(), test_data.end(), 10*comm.rank());
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE, comm);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
    }
    {
        Opm::HDF5File in_file(rwpath, Opm::HDF5File::OpenMode::READ, comm);
        std::vector<char> data;
        BOOST_CHECK_NO_THROW(in_file.read("/test_data", "d1", data));
        BOOST_CHECK_EQUAL_COLLECTIONS(data.begin(), data.end(),
                                      test_data.begin(), test_data.end());
    }
    std::filesystem::remove(rwpath);
    std::filesystem::remove(path);
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
