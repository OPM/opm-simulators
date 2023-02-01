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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/utils/ParallelFileMerger.hpp>
#include <fstream>
#include <iostream>

namespace Opm {
namespace detail {

ParallelFileMerger::ParallelFileMerger(const fs::path& output_dir,
                                       const std::string& deckname,
                                       bool show_fallout)
    : debugFileRegex_(deckname+"\\.\\d+\\.DBG")
    , logFileRegex_(deckname+"\\.\\d+\\.PRT")
    , fileWarningRegex_(deckname+"\\.(\\d+)\\.[^.]+")
    , show_fallout_(show_fallout)
{
    if (show_fallout_) {
        auto debugPath = output_dir;
        debugPath /= (deckname + ".DBG");
        debugStream_ = std::make_unique<std::ofstream>(debugPath,
                                                       std::ofstream::app);
        auto logPath = output_dir;
        logPath /= ( deckname + ".PRT");
        logStream_ = std::make_unique<std::ofstream>(logPath,
                                                     std::ofstream::app);
    }
}

void ParallelFileMerger::operator()(const fs::path& file)
{
    std::smatch matches;
    std::string filename = file.filename().native();

    if (std::regex_match(filename, matches, fileWarningRegex_)) {
        std::string rank = std::regex_replace(filename, fileWarningRegex_, "$1");

        if (std::regex_match(filename, logFileRegex_)) {
            if (show_fallout_) {
                appendFile(*logStream_, file, rank);
            } else {
                fs::remove(file);
            }
        } else {
            if (std::regex_match(filename, debugFileRegex_)) {
                if (show_fallout_) {
                    appendFile(*debugStream_, file, rank);
                } else {
                    fs::remove(file);
                }
            } else {
                if (show_fallout_) {
                    std::cerr << "WARNING: Unrecognized file with name "
                              << filename
                              << " that might stem from a parallel run."
                              << std::endl;
                }
            }
        }
    }
}

void ParallelFileMerger::appendFile(std::ofstream& of, const fs::path& file, const std::string& rank)
{
    if (fs::file_size(file)) {
        std::cerr << "WARNING: There has been logging to file "
                      << file.string() <<" by process "
                      << rank << std::endl;

        std::ifstream in(file);
        of<<std::endl<< std::endl;
        of<<"=======================================================";
        of<<std::endl<<std::endl;
        of << " Output written by rank " << rank << " to file " << file.string();
        of << ":" << std::endl << std::endl;
        of << in.rdbuf() << std::endl << std::endl;
        of << "======================== end output =====================";
        of << std::endl;
        in.close();
    }
    fs::remove(file);
}

} // end namespace detail
} // end namespace Opm
