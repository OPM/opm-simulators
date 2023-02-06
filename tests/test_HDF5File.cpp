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
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <stdexcept>

using namespace Opm;

BOOST_AUTO_TEST_CASE(ReadWrite)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "rw.hdf5").string();
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
    }
    {
        Opm::HDF5File in_file(rwpath, Opm::HDF5File::OpenMode::READ);
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
    BOOST_CHECK_THROW(Opm::HDF5File out_file("no_such_file.hdf5", Opm::HDF5File::OpenMode::READ), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ReadNonExistentDset)
{
    auto path = std::filesystem::temp_directory_path() / Opm::unique_path("hdf5test%%%%%");
    std::filesystem::create_directory(path);
    auto rwpath = (path / "existent_dset.hdf5").string();
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
    }
    {
        Opm::HDF5File in_file(rwpath, Opm::HDF5File::OpenMode::READ);
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
    const std::vector<char> test_data{1,2,3,4,5,6,8,9};
    {
        Opm::HDF5File out_file(rwpath, Opm::HDF5File::OpenMode::OVERWRITE);
        BOOST_CHECK_NO_THROW(out_file.write("/test_data", "d1", test_data));
        BOOST_CHECK_THROW(out_file.write("/test_data", "d1", test_data), std::runtime_error);
    }
    std::filesystem::remove(rwpath);
    std::filesystem::remove(path);
}

