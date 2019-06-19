/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015 IRIS AS

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

#ifndef OPM_BLACKOILDETAILS_HEADER_INCLUDED
#define OPM_BLACKOILDETAILS_HEADER_INCLUDED

#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

namespace Opm {
namespace detail {


    inline
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }



    template <class PU>
    std::vector<bool>
    activePhases(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<bool> active(maxnp, false);

        for (int p = 0; p < pu.MaxNumPhases; ++p) {
            active[ p ] = pu.phase_used[ p ] != 0;
        }

        return active;
    }



    template <class PU>
    std::vector<int>
    active2Canonical(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<int> act2can(maxnp, -1);

        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[ phase ]) {
                act2can[ pu.phase_pos[ phase ] ] = phase;
            }
        }

        return act2can;
    }



    inline
    double getGravity(const double* g, const int dim) {
        double grav = 0.0;
        if (g) {
            // Guard against gravity in anything but last dimension.
            for (int dd = 0; dd < dim - 1; ++dd) {
                assert(g[dd] == 0.0);
            }
            grav = g[dim - 1];
        }
        return grav;
    }


        /// \brief Compute the Euclidian norm of a vector
        /// \warning In the case that num_components is greater than 1
        ///          an interleaved ordering is assumed. E.g. for each cell
        ///          all phases of that cell are stored consecutively. First
        ///          the ones for cell 0, then the ones for cell 1, ... .
        /// \param it              begin iterator for the given vector
        /// \param end             end iterator for the given vector
        /// \param num_components  number of components (i.e. phases) in the vector
        /// \param pinfo           In a parallel this holds the information about the data distribution.
        template <class Iterator>
        inline
        double euclidianNormSquared( Iterator it, const Iterator end, int num_components, const boost::any& pinfo = boost::any() )
        {
            static_cast<void>(num_components); // Suppress warning in the serial case.
            static_cast<void>(pinfo); // Suppress warning in non-MPI case.
#if HAVE_MPI
            if ( pinfo.type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(pinfo);
                typedef typename Iterator::value_type Scalar;
                Scalar product = 0.0;
                int size_per_component = (end - it);
                size_per_component /= num_components; // two lines to supresse unused warning.
                assert((end - it) == num_components * size_per_component);

                if( num_components == 1 )
                {
                    auto component_container =
                        boost::make_iterator_range(it, end);
                    info.computeReduction(component_container,
                                           Opm::Reduction::makeInnerProductFunctor<double>(),
                                           product);
                }
                else
                {
                    auto& maskContainer = info.getOwnerMask();
                    auto mask = maskContainer.begin();
                    assert(static_cast<int>(maskContainer.size()) == size_per_component);

                    for(int cell = 0; cell < size_per_component; ++cell, ++mask)
                    {
                        Scalar cell_product = (*it) * (*it);
                        ++it;
                        for(int component=1; component < num_components;
                            ++component, ++it)
                        {
                            cell_product += (*it) * (*it);
                        }
                        product += cell_product * (*mask);
                    }
                }
                return info.communicator().sum(product);
            }
            else
#endif
            {
                double product = 0.0 ;
                for( ; it != end; ++it ) {
                    product += ( *it * *it );
                }
                return product;
            }
        }
        /// \brief Get the number of local interior cells in a grid.
        /// \tparam The type of the DUNE grid.
        /// \param grid The grid which cells we count
        /// \return The number of interior cell in the partition of the
        /// grid stored on this process.
        template<class Grid>
        std::size_t countLocalInteriorCells(const Grid& grid)
        {
            if ( grid.comm().size() == 1)
            {
                return grid.size(0);
            }
            std::size_t count = 0;
            const auto& gridView = grid.leafGridView();
            for(auto cell = gridView.template begin<0, Dune::Interior_Partition>(),
                    endCell = gridView.template end<0, Dune::Interior_Partition>();
                cell != endCell; ++cell)
            {
                    ++count;
            }
            return count;
        }

        /// \brief Get the number of cells of a global grid.
        ///
        /// In a parallel run this is the number of cells that a grid would
        /// have if the whole grid was stored on one process only.
        /// \tparam The type of the DUNE grid.
        /// \param grid The grid which cells we count
        /// \return The global number of cells.
        template<class Grid>
        std::size_t countGlobalCells(const Grid& grid)
        {
            if ( grid.comm().size() == 1)
            {
                return grid.size(0);
            }
            std::size_t count = countLocalInteriorCells(grid);
            return grid.comm().sum(count);
        }

    } // namespace detail
} // namespace Opm

#endif // OPM_BLACKOILDETAILS_HEADER_INCLUDED
