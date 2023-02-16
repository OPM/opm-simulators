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

#define BOOST_TEST_MODULE HDF5FileTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <stdexcept>

BOOST_AUTO_TEST_CASE(ReadWrite)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "rw.hdf5").string();
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
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

BOOST_AUTO_TEST_CASE(ThrowOpenNonexistent)
{
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif
    BOOST_CHECK_THROW(Opm::HDF5File out_file("no_such_file.hdf5", Opm::HDF5File::OpenMode::READ, comm), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ReadNonExistentDset)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "existent_dset.hdf5").string();
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
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
        BOOST_CHECK_THROW(in_file.read("/test_data", "d2", data), std::runtime_error);
    }
    std::filesystem::remove(rwpath);
    std::filesystem::remove(path);
}

BOOST_AUTO_TEST_CASE(WriteExistentDset)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "existent_dset.hdf5").string();
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE, comm);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
        BOOST_CHECK_THROW(out_file.write("/test_data", "d1", test_data), std::runtime_error);
    }
    std::filesystem::remove(rwpath);
    std::filesystem::remove(path);
}

BOOST_AUTO_TEST_CASE(List)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "existent_dset.hdf5").string();
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE, comm);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d2", test_data));
        BOOST_CHECK_NO_THROW(out_file.write("/test_data/test", "d2", test_data));
    }
    {
        Opm::HDF5File in_file(rwpath, Opm::HDF5File::OpenMode::READ, comm);

        auto res1 = in_file.list("/");
        BOOST_CHECK_EQUAL(res1.size(), 1u);
        BOOST_CHECK_EQUAL(res1[0], "test_data");

        auto res2 = in_file.list("/test_data");
        BOOST_CHECK_EQUAL(res2.size(), 3u);
        BOOST_CHECK_EQUAL(res2[0], "d1");
        BOOST_CHECK_EQUAL(res2[1], "d2");
        BOOST_CHECK_EQUAL(res2[2], "test");

        BOOST_CHECK_THROW(in_file.list("/not_there"), std::runtime_error);
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
