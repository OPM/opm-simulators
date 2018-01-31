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

#include "config.h"
#include <opm/core/props/rock/RockBasic.hpp>

namespace Opm
{


    /// Default constructor.
    RockBasic::RockBasic()
        : dimensions_(-1)
    {
    }


    /// Initialize with homogenous porosity and permeability.
    void RockBasic::init(const int dimensions,
                         const int num_cells,
                         const double poro,
                         const double perm)
    {
        dimensions_ = dimensions;
        porosity_.clear();
        porosity_.resize(num_cells, poro);
        permeability_.clear();
        const int dsq = dimensions*dimensions;
        permeability_.resize(num_cells*dsq, 0.0);
// #pragma omp parallel for
        for (int i = 0; i < num_cells; ++i) {
            for (int d = 0; d < dimensions; ++d) {
                permeability_[dsq*i + dimensions*d + d] = perm;
            }
        }
    }


} // namespace Opm
