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


#include <opm/core/io/eclipse/EclipseGridParser.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class RockFromDeck
    {
    public:
        /// Default constructor.
        RockFromDeck();

        /// Initialize from deck and grid.
        /// \param  deck         Deck input parser
        /// \param  grid         Grid to which property object applies, needed for the
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        void init(const EclipseGridParser& deck,
                  const UnstructuredGrid& grid);

        /// Initialize from deck and cell mapping.
        /// \param  deck            Deck input parser
        /// \param  number_of_cells The number of cells in the grid.
        /// \param  global_cell     The mapping fom local to global cell indices.
        ///                         global_cell[i] is the corresponding global index of i.
        /// \param  cart_dims       The size of the underlying cartesian grid.
        void init(const EclipseGridParser& deck,
                  int number_of_cells, const int* global_cell,
                  const int* cart_dims);
         /// Initialize from deck and cell mapping.
        /// \param  newParserDeck            Deck produced by the opm-parser code
        /// \param  number_of_cells The number of cells in the grid.
        /// \param  global_cell     The mapping fom local to global cell indices.
        ///                         global_cell[i] is the corresponding global index of i.
        /// \param  cart_dims       The size of the underlying cartesian grid.
        void init(Opm::DeckConstPtr newParserDeck,
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

    private:
        void assignPorosity(const EclipseGridParser& parser,
                            int number_of_cells,
                            const int* global_cell);
        void assignPorosity(Opm::DeckConstPtr newParserDeck,
                            int number_of_cells,
                            const int* global_cell);
        void assignPermeability(const EclipseGridParser& parser,
                                int number_of_cells,
                                const int* global_cell,
                                const int* cart_dims,
                                const double perm_threshold);
        void assignPermeability(Opm::DeckConstPtr newParserDeck,
                                int number_of_cells,
                                const int* global_cell,
                                const int* cart_dims,
                                double perm_threshold);

        std::vector<double> porosity_;
        std::vector<double> permeability_;
        std::vector<unsigned char> permfield_valid_;
    };



} // namespace Opm


#endif // OPM_ROCKFROMDECK_HEADER_INCLUDED
