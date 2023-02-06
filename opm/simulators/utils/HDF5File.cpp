/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/utils/HDF5File.hpp>

#include <filesystem>
#include <stdexcept>

namespace {

bool groupExists(hid_t parent, const std::string& path)
{
  // turn off errors to avoid cout spew
  H5E_BEGIN_TRY {
#if H5_VERS_MINOR > 8
      return H5Lexists(parent, path.c_str(), H5P_DEFAULT) == 1;
#else
      return H5Gget_objinfo(static_cast<hid_t>(parent), path.c_str(), 0, nullptr) == 0;
#endif
  } H5E_END_TRY;
  return false;
}

}

namespace Opm {

HDF5File::HDF5File(const std::string& fileName, OpenMode mode)
{
    bool exists = std::filesystem::exists(fileName);
    if (mode == OpenMode::OVERWRITE ||
        (mode == OpenMode::APPEND && !exists)) {
        m_file = H5Fcreate(fileName.c_str(),
                           H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    } else {
        m_file = H5Fopen(fileName.c_str(),
                         mode == OpenMode::READ ? H5F_ACC_RDONLY : H5F_ACC_RDWR,
                         H5P_DEFAULT);
    }
    if (m_file == H5I_INVALID_HID) {
        throw std::runtime_error(std::string("HDF5File: Failed to ") +
                                 ( mode == OpenMode::OVERWRITE ||
                                  (mode == OpenMode::APPEND && !exists) ? "create" : "open") +
                                  fileName);
    }
}

HDF5File::~HDF5File()
{
    if (m_file != H5I_INVALID_HID) {
        H5Fclose(m_file);
    }
}

void HDF5File::write(const std::string& group,
                     const std::string& dset,
                     const std::vector<char>& buffer)
{
    hid_t grp;
    if (groupExists(m_file, group)) {
        grp = H5Gopen2(m_file, group.c_str(), H5P_DEFAULT);
    } else {
        grp = H5Gcreate2(m_file, group.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);
    }

    hsize_t size = buffer.size();
    hsize_t start = 0;

    hid_t space = H5Screate_simple(1, &size, nullptr);
    hid_t dataset_id = H5Dcreate2(grp, dset.c_str(), H5T_NATIVE_CHAR, space,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        throw std::runtime_error("HDF5File: Trying to write already existing dataset " + group + '/' + dset);
    }

    if (size > 0) {
        hid_t filespace = H5Dget_space(dataset_id);
        hsize_t stride = 1;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, &stride, &size, nullptr);
        hid_t memspace = H5Screate_simple(1, &size, nullptr);
        H5Dwrite(dataset_id, H5T_NATIVE_CHAR, memspace, filespace, H5P_DEFAULT, buffer.data());
        H5Sclose(memspace);
        H5Sclose(filespace);
    }
    H5Dclose(dataset_id);
    H5Sclose(space);
    H5Gclose(grp);
}

void HDF5File::read(const std::string& group,
                    const std::string& dset,
                    std::vector<char>& buffer)
{
    hid_t dataset_id = H5Dopen2(m_file, (group + "/"+ dset).c_str(), H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        throw std::runtime_error("HDF5File: Trying to read non-existing dataset " + group + '/' + dset);
    }

    hid_t space = H5Dget_space(dataset_id);
    hsize_t size = H5Sget_simple_extent_npoints(space);
    buffer.resize(size);
    H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
    H5Dclose(dataset_id);
}

}
