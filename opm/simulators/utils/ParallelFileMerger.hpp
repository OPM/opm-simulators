/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 STATOIL AS.

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

#ifndef OPM_PARALLELFILEMERGER_HEADER_INCLUDED
#define OPM_PARALLELFILEMERGER_HEADER_INCLUDED

#include <filesystem>
#include <iosfwd>
#include <memory>
#include <regex>
#include <string>

namespace Opm {
namespace detail {

namespace fs = ::std::filesystem;

/// \brief A functor that merges multiple files of a parallel run to one file.
///
/// Without care multiple processes might log messages in a parallel run.
/// Non-root processes will do that to seperate files
/// <basename>.<rank>.<extension. This functor will append those file
/// to usual ones and delete the other files.
class ParallelFileMerger
{
public:
    /// \brief Constructor
    /// \param output_dir The output directory to use for reading/Writing.
    /// \param deckname The name of the deck.
    ParallelFileMerger(const fs::path& output_dir,
                       const std::string& deckname,
                       bool show_fallout = false);

    void operator()(const fs::path& file);

private:
    /// \brief Append contents of a file to a stream
    /// \brief of The output stream to use.
    /// \brief file The file whose content to append.
    /// \brief rank The rank that wrote the file.
    void appendFile(std::ofstream& of, const fs::path& file, const std::string& rank);

    /// \brief Regex to capture *.DBG
    std::regex debugFileRegex_;
    /// \brief Regex to capture  *.PRT
    std::regex logFileRegex_;
    /// \brief Regex to capture  CASENAME.[0-9]+.[A-Z]+
    std::regex fileWarningRegex_;
    /// \brief Stream to *.DBG file
    std::unique_ptr<std::ofstream> debugStream_;
    /// \brief Stream to *.PRT file
    std::unique_ptr<std::ofstream> logStream_;
    /// \brief Whether to show any logging fallout
    bool show_fallout_;
};

} // end namespace detail
} // end namespace Opm

#endif // end header guard
