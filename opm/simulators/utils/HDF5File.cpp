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

#include <opm/common/utility/String.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <cassert>
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

HDF5File::HDF5File(const std::string& fileName,
                   OpenMode mode,
                   Parallel::Communication comm)
    : comm_(comm)
{
    bool exists = std::filesystem::exists(fileName);
    hid_t acc_tpl = H5P_DEFAULT;
    if (comm.size() > 1) {
#if HAVE_MPI
        MPI_Info info = MPI_INFO_NULL;
        acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(acc_tpl, comm_, info);
#else
        assert(false); // should be unreachable
#endif
    }

    if (mode == OpenMode::OVERWRITE ||
        (mode == OpenMode::APPEND && !exists)) {
        m_file = H5Fcreate(fileName.c_str(),
                           H5F_ACC_TRUNC,
                           H5P_DEFAULT, acc_tpl);
    } else {
        m_file = H5Fopen(fileName.c_str(),
                         mode == OpenMode::READ ? H5F_ACC_RDONLY : H5F_ACC_RDWR,
                         acc_tpl);
    }
    if (m_file == H5I_INVALID_HID) {
        throw std::runtime_error(std::string("HDF5File: Failed to ") +
                                 ( mode == OpenMode::OVERWRITE ||
                                  (mode == OpenMode::APPEND && !exists) ? "create" : "open") +
                                  fileName);
    }

    if (comm_.size() > 1) {
        H5Pclose(acc_tpl);
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
                     const std::vector<char>& buffer,
                     DataSetMode mode) const
{
    hid_t grp = H5I_INVALID_HID;
    std::string realGroup = group;
    if (mode == DataSetMode::PROCESS_SPLIT) {
        if (group != "/")
            realGroup += '/';
        realGroup += dset;
    }

    OPM_BEGIN_PARALLEL_TRY_CATCH();

    if (groupExists(m_file, realGroup)) {
        grp = H5Gopen2(m_file, realGroup.c_str(), H5P_DEFAULT);
    } else {
        auto grps = split_string(realGroup, '/');
        std::string curr;
        for (size_t i = 0; i < grps.size(); ++i) {
            if (grps[i].empty())
                continue;
            curr += '/';
            curr += grps[i];
            if (!groupExists(m_file, curr)) {
                hid_t subgrp = H5Gcreate2(m_file, curr.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);
                if (subgrp == H5I_INVALID_HID) {
                    throw std::runtime_error("Failed to create group '" + curr + "'");
                }
                if (i == grps.size() - 1) {
                    grp = subgrp;
                } else {
                    H5Gclose(subgrp);
                }
            } else if (i == grps.size() - 1) {
                grp = H5Gopen2(m_file, realGroup.c_str(), H5P_DEFAULT);
            }
        }
    }

    if (grp == H5I_INVALID_HID) {
        throw std::runtime_error("Failed to create group '" + realGroup + "'");
    }

    if (mode == DataSetMode::PROCESS_SPLIT) {
        writeSplit(grp, buffer, realGroup);
    } else if (mode == DataSetMode::ROOT_ONLY) {
        writeRootOnly(grp, buffer, group, dset);
    }
    H5Gclose(grp);

    OPM_END_PARALLEL_TRY_CATCH("HDF5File: Error writing data: ", comm_);
}

void HDF5File::read(const std::string& group,
                    const std::string& dset,
                    std::vector<char>& buffer,
                    DataSetMode mode) const
{
    std::string realSet = group + '/' + dset;
    if (mode == DataSetMode::PROCESS_SPLIT) {
        realSet += '/' + std::to_string(comm_.rank());
    }
    hid_t dataset_id = H5Dopen2(m_file, realSet.c_str(), H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        throw std::runtime_error("Trying to read non-existing dataset " + group + '/' + dset);
    }

    hid_t space = H5Dget_space(dataset_id);
    hsize_t size = H5Sget_simple_extent_npoints(space);
    buffer.resize(size);
    H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
    H5Dclose(dataset_id);
}

std::vector<std::string> HDF5File::list(const std::string& group) const
{
    // Lambda function pushing the group entries to a vector
    auto&& list_group = [] (hid_t, const char* name, const H5L_info_t*, void* data) -> herr_t
    {
        auto& list = *static_cast<std::vector<std::string>*>(data);
        list.push_back(name);
        return 0;
    };

    hsize_t idx = 0;
    std::vector<std::string> result;
    if (H5Literate_by_name(m_file, group.c_str(),
                           H5_INDEX_NAME, H5_ITER_INC,
                           &idx, list_group, &result, H5P_DEFAULT) < 0) {
        throw std::runtime_error("Failure while listing group '" + group + "'");
    }

    return result;
}

void HDF5File::writeSplit(hid_t grp,
                          const std::vector<char>& buffer,
                          const std::string& dset) const
{
    std::vector<hsize_t> proc_sizes(comm_.size());
    if (comm_.size() > 1) {
        hsize_t lsize = buffer.size();
        comm_.allgather(&lsize, 1, proc_sizes.data());
    } else {
        proc_sizes[0] = buffer.size();
    }

    for (int i = 0; i < comm_.size(); ++i) {
        hid_t space = H5Screate_simple(1, &proc_sizes[i], nullptr);
        hid_t dataset_id = H5Dcreate2(grp,
                                      std::to_string(i).c_str(),
                                      H5T_NATIVE_CHAR, space,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id == H5I_INVALID_HID) {
            throw std::runtime_error("Trying to write already existing dataset '" +
                                     dset + '/' + std::to_string(i) + "'");
        }

        if (i == comm_.rank()) {
            hid_t filespace = H5Dget_space(dataset_id);
            hsize_t stride = 1;
            hsize_t start = 0;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, &stride, &proc_sizes[i], nullptr);
            hid_t memspace = H5Screate_simple(1, &proc_sizes[i], nullptr);
            H5Dwrite(dataset_id, H5T_NATIVE_CHAR, memspace, filespace, H5P_DEFAULT, buffer.data());
            H5Sclose(memspace);
            H5Sclose(filespace);
        }
        H5Dclose(dataset_id);
        H5Sclose(space);
    }
}

void HDF5File::writeRootOnly(hid_t grp,
                             const std::vector<char>& buffer,
                             const std::string& group,
                             const std::string& dset) const
{
    hsize_t size = buffer.size();
    comm_.broadcast(&size, 1, 0);
    hid_t space = H5Screate_simple(1, &size, nullptr);
    hid_t dataset_id = H5Dcreate2(grp,
                                  dset.c_str(),
                                  H5T_NATIVE_CHAR, space,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        throw std::runtime_error("Trying to write already existing dataset '" +
                                 group + '/' + dset + "'");
    }

    if (comm_.rank() == 0) {
        hid_t filespace = H5Dget_space(dataset_id);
        hsize_t stride = 1;
        hsize_t start = 0;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, &stride, &size, nullptr);
        hid_t memspace = H5Screate_simple(1, &size, nullptr);
        H5Dwrite(dataset_id, H5T_NATIVE_CHAR, memspace, filespace, H5P_DEFAULT, buffer.data());
        H5Sclose(memspace);
        H5Sclose(filespace);
    }
    H5Dclose(dataset_id);
    H5Sclose(space);
}

}
