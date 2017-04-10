/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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

#ifndef OPM_ENSUREDIRECTORYEXISTS_HEADER_INCLUDED
#define OPM_ENSUREDIRECTORYEXISTS_HEADER_INCLUDED

#include <boost/filesystem.hpp>

namespace Opm
{

    /// The directory pointed to by 'dirpath' will be created if it
    /// does not already exist. Will throw an exception if this cannot
    /// be done.
    ///
    /// Note that std::string can be passed to this functions, as they
    /// can be implicitly converted to boost::filesystem::path objects.
    void ensureDirectoryExists(const boost::filesystem::path& dirpath);


} // namespace Opm

#endif // OPM_ENSUREDIRECTORYEXISTS_HEADER_INCLUDED
