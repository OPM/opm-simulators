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

#ifndef OPM_ROCKFROMDECK_HEADER_INCLUDED
#define OPM_ROCKFROMDECK_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class RockFromDeck
    {
        // BlackoilPropsDataHandle needs mutable
        // access to porosity and permeability
        friend class BlackoilPropsDataHandle;

    public:
        /// Default constructor.
        RockFromDeck();
        /// Creates rock properties with zero porosity and permeability
        /// \param number_of_cells The number of cells
        explicit RockFromDeck(std::size_t number_of_cells);
        /// Initialize from deck and cell mapping.
        /// \param  eclState        The EclipseState (processed deck) produced by the opm-parser code
        /// \param  number_of_cells The number of cells in the grid.
        /// \param  global_cell     The mapping fom local to global cell indices.
        ///                         global_cell[i] is the corresponding global index of i.
        /// \param  cart_dims       The size of the underlying cartesian grid.
        void init(const Opm::EclipseState& eclState,
                  int number_of_cells, const int* global_cell,
                  const int* cart_dims);

        /// \return   D, the number of spatial dimensions. Always 3 for deck input.
        int numDimensions() const
        {
            return 3;
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

        /// Convert the permeabilites for the logically Cartesian grid in EclipseState to
        /// an array of size number_of_cells*dim*dim for the compressed array.
        /// \param  eclState        The EclipseState (processed deck) produced by the opm-parser code
        /// \param  number_of_cells The number of cells in the grid.
        /// \param  global_cell     The mapping fom local to global cell indices.
        ///                         global_cell[i] is the corresponding global index of i.
        /// \param  cart_dims       The size of the underlying cartesian grid.
        /// \param  perm_threshold  The threshold for permeability
        /// \param  permeability    The result array
        static
        void extractInterleavedPermeability(const Opm::EclipseState& eclState,
                                            const int number_of_cells,
                                            const int* global_cell,
                                            const int* cart_dims,
                                            const double perm_threshold,
                                            std::vector<double>& permeability);

    private:
        void assignPorosity(const Opm::EclipseState& eclState,
                            int number_of_cells,
                            const int* global_cell);

        std::vector<double> porosity_;
        std::vector<double> permeability_;
        std::vector<unsigned char> permfield_valid_;
    };



} // namespace Opm


#endif // OPM_ROCKFROMDECK_HEADER_INCLUDED
