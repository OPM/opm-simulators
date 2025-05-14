// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#include <opm/models/io/restart.hpp>

namespace Opm {

std::string Restart::restartFileName_(int rank,
                                      const std::string& outputDir,
                                      const std::string& simName,
                                      double t)
{
    std::string dir = outputDir;
    if (dir == ".") {
        dir = "";
    }
    else if (!dir.empty() && dir.back() != '/') {
        dir += "/";
    }

    std::ostringstream oss;
    oss << dir << simName << "_time=" << t << "_rank=" << rank << ".ers";
    return oss.str();
}

void Restart::serializeSectionBegin(const std::string& cookie)
{
    outStream_ << cookie << "\n";
}

void Restart::serializeSectionEnd()
{
    outStream_ << "\n";
}

void Restart::deserializeSectionBegin(const std::string& cookie)
{
    if (!inStream_.good()) {
        throw std::runtime_error("Encountered unexpected EOF in restart file.");
    }
    std::string buf;
    std::getline(inStream_, buf);
    if (buf != cookie) {
        throw std::runtime_error("Could not start section '"+cookie+"'");
    }
}

void Restart::deserializeSectionEnd()
{
    std::string dummy;
    std::getline(inStream_, dummy);
    for (unsigned i = 0; i < dummy.length(); ++i) {
        if (!std::isspace(dummy[i])) {
            throw std::logic_error("Encountered unread values while deserializing");
        }
    }
}

void Restart::deserializeEnd()
{
    inStream_.close();
}

void Restart::serializeEnd()
{
    outStream_.close();
}

void Restart::openInputStream(const std::string& cookie)
{
    // open input file and read magic cookie
    inStream_.open(fileName_.c_str());
    if (!inStream_.good()) {
        throw std::runtime_error("Restart file '" + fileName_ +
                                 "' could not be opened properly");
    }

    // make sure that we don't open an empty file
    inStream_.seekg(0, std::ios::end);
    auto pos = inStream_.tellg();
    if (pos == 0) {
        throw std::runtime_error("Restart file '" + fileName_ + "' is empty");
    }
    inStream_.seekg(0, std::ios::beg);

    deserializeSectionBegin(cookie);
    deserializeSectionEnd();
}

void Restart::openOutputStream(const std::string& cookie)
{
    outStream_.open(fileName_);
    outStream_.precision(20);

    serializeSectionBegin(cookie);
    serializeSectionEnd();
}

} // namespace Opm
