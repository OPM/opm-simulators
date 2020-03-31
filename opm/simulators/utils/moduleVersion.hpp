/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_MODULEVERSION_HEADER_INCLUDED
#define OPM_MODULEVERSION_HEADER_INCLUDED

#include <string>

namespace Opm
{

    /// Return the version name of the module, for example "2015.10"
    /// (for a release branch) or "2016.04-pre" (for a master branch).
    std::string moduleVersionName();

    /// Return a (short) git hash for the current version of the
    /// module if this is a Release build (as defined by CMake), or
    /// "debug" for Debug builds.
    std::string moduleVersionHash();

    /// Return a string containing both the name and hash, if N is the
    /// name and H is the hash it will be "N (H)". For example
    /// "2016.04-pre (f15be17)" or "2016.04-pre (debug)".
    std::string moduleVersion();

    /// Return a string "dd-mm-yyyy at HH::MM::SS hrs" which is the time
    /// the binary was compiled.
    std::string compileTimestamp();

} // namespace Opm

#endif // OPM_MODULEVERSION_HEADER_INCLUDED
