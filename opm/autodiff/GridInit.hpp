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

#ifndef OPM_GRIDINIT_HEADER_INCLUDED
#define OPM_GRIDINIT_HEADER_INCLUDED

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/core/grid/GridManager.hpp>

#if HAVE_DUNE_CORNERPOINT
#include <dune/grid/CpGrid.hpp>
#endif


namespace Opm
{

    /// A class intended to give a generic interface to
    /// initializing and accessing UnstructuredGrid and CpGrid,
    /// using specialized templates to accomplish this.
    template <class Grid>
    class GridInit
    {
    public:
        /// Initialize from a deck and/or an eclipse state and (logical cartesian) specified pore volumes.
        GridInit(DeckConstPtr, EclipseStateConstPtr, const std::vector<double>&)
        {
            OPM_THROW(std::logic_error, "Found no specialization for GridInit for the requested Grid class.");
        }
    };


    /// Specialization for UnstructuredGrid.
    template <>
    class GridInit<UnstructuredGrid>
    {
    public:
        /// Initialize from a deck and/or an eclipse state and (logical cartesian) specified pore volumes.
        GridInit(DeckConstPtr, EclipseStateConstPtr eclipse_state, const std::vector<double>& porv)
            : grid_manager_(eclipse_state->getEclipseGrid(), porv)
        {
        }
        /// Access the created grid.
        const UnstructuredGrid& grid()
        {
            return *grid_manager_.c_grid();
        }
    private:
        GridManager grid_manager_;
    };


#if HAVE_DUNE_CORNERPOINT
    /// Specialization for CpGrid.
    template <>
    class GridInit<Dune::CpGrid>
    {
    public:
        /// Initialize from a deck and/or an eclipse state and (logical cartesian) specified pore volumes.
        GridInit(DeckConstPtr deck, EclipseStateConstPtr, const std::vector<double>& porv)
        {
            grid_.processEclipseFormat(deck, false, false, false, porv);
        }
        /// Access the created grid. Note that mutable access may be required for load balancing.
        Dune::CpGrid& grid()
        {
            return grid_;
        }
    private:
        Dune::CpGrid grid_;
    };
#endif // HAVE_DUNE_CORNERPOINT


} // namespace Opm

#endif // OPM_GRIDINIT_HEADER_INCLUDED
