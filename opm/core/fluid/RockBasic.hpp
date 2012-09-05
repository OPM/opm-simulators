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

#ifndef OPM_ROCKBASIC_HEADER_INCLUDED
#define OPM_ROCKBASIC_HEADER_INCLUDED


#include <vector>


namespace Opm
{

    class RockBasic
    {
    public:
        /// Default constructor.
        RockBasic();

        /// Initialize with homogenous porosity and permeability.
        void init(const int dimensions,
                  const int num_cells,
                  const double poro,
                  const double perm);

        /// \return   D, the number of spatial dimensions.
        int numDimensions() const
        {
            return dimensions_;
        }

        /// \return   N, the number of cells.
        int numCells() const
        {
            return porosity_.size();
        }

        /// \return   Array of N porosity values.
        const double* porosity() const
        {
            return &porosity_[0];
        }

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        const double* permeability() const
        {
            return &permeability_[0];
        }

    private:
        int dimensions_;
        std::vector<double> porosity_;
        std::vector<double> permeability_;
    };



} // namespace Opm


#endif // OPM_ROCKBASIC_HEADER_INCLUDED
