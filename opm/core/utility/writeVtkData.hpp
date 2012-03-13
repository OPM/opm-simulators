/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WRITEVTKDATA_HEADER_INCLUDED
#define OPM_WRITEVTKDATA_HEADER_INCLUDED


#include <string>
#include <map>
#include <vector>
#include <tr1/array>
#include <iosfwd>

struct UnstructuredGrid;

namespace Opm
{

    /// Intended to map strings (giving the output field names) to data.
    typedef std::map<std::string, const std::vector<double>*> DataMap;

    /// Vtk output for cartesian grids.
    void writeVtkData(const std::tr1::array<int, 3>& dims,
		      const std::tr1::array<double, 3>& cell_size,
		      const DataMap& data,
		      std::ostream& os);

    /// Vtk output for general grids.
    void writeVtkData(const UnstructuredGrid& grid,
		      const DataMap& data,
		      std::ostream& os);
} // namespace Opm

#endif // OPM_WRITEVTKDATA_HEADER_INCLUDED
