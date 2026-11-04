/*
  Copyright 2026, SINTEF Digital

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
/*!
 * \file
 * \brief A temporary directory for tests that is removed when they end
 */
#ifndef OPM_TEMPORARY_DIRECTORY_HPP
#define OPM_TEMPORARY_DIRECTORY_HPP

#include <opm/common/utility/FileSystem.hpp>

#include <filesystem>
#include <string>
#include <system_error>

namespace Opm {

/*!
 * \brief A uniquely named directory under the system's temporary directory
 *
 * The directory and its contents are removed when the object goes out of
 * scope, also when a test case ends early on a failed requirement or an
 * exception.
 */
class TemporaryDirectory
{
public:
    explicit TemporaryDirectory(const std::string& prefix)
        : path_(std::filesystem::temp_directory_path() / unique_path(prefix + "-%%%%%%%%"))
    {
        std::filesystem::create_directories(path_);
    }

    ~TemporaryDirectory()
    {
        // A destructor must not throw; a directory left behind is harmless.
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    const std::filesystem::path& path() const
    { return path_; }

private:
    std::filesystem::path path_;
};

} // namespace Opm

#endif // OPM_TEMPORARY_DIRECTORY_HPP
