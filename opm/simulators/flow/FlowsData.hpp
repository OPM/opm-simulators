/*
  Copyright 2023 SINTEF Digital

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

#ifndef OPM_FLOWS_DATA_HEADER_INCLUDED
#define OPM_FLOWS_DATA_HEADER_INCLUDED

#include <cstddef>
#include <string>
#include <vector>

namespace Opm {

//! \brief Simple container for FLOWS data.
template<class Scalar>
struct FlowsData
{
    //! \brief Resize data vectors.
    void resize(const std::size_t size)
    {
        indices.resize(size);
        values.resize(size);
    }

    std::string name; //!< Associated name
    std::vector<int> indices; //!< Cell indices for values
    std::vector<Scalar> values; //!< Values
};

} // namespace Opm

#endif // OPM_FLOWS_DATA_HEADER_INCLUDED
