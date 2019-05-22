/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
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

#ifndef OPM_WELLCONNECTIONAUXILIARYMODULE_HEADER_INCLUDED
#define OPM_WELLCONNECTIONAUXILIARYMODULE_HEADER_INCLUDED

#include <ewoms/disc/common/baseauxiliarymodule.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/core/wells.h>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

namespace Opm
{
template<class TypeTag>
class WellConnectionAuxiliaryModule
    : public Ewoms::BaseAuxiliaryModule<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

public:

    using NeighborSet = typename
        Ewoms::BaseAuxiliaryModule<TypeTag>::NeighborSet;

    WellConnectionAuxiliaryModule(const Schedule& schedule,
                                  const Dune::CpGrid& grid)
    {
        // Create cartesian to compressed mapping
        const auto& globalCell = grid.globalCell();
        const auto& cartesianSize = grid.logicalCartesianSize();

        auto size = cartesianSize[0]*cartesianSize[1]*cartesianSize[2];

        std::vector<int> cartesianToCompressed(size, -1);
        auto begin = globalCell.begin();

        for ( auto cell = begin, end= globalCell.end(); cell != end; ++cell )
        {
          cartesianToCompressed[ *cell ] = cell - begin;
        }

        int last_time_step = schedule.getTimeMap().size() - 1;
        const auto& schedule_wells = schedule.getWells2atEnd();
        wells_.reserve(schedule_wells.size());

        // initialize the additional cell connections introduced by wells.
        for ( const auto well : schedule_wells )
        {
            std::vector<int> compressed_well_perforations;
            // All possible completions of the well
            const auto& completionSet = well.getConnections();
            compressed_well_perforations.reserve(completionSet.size());

            for ( size_t c=0; c < completionSet.size(); c++ )
            {
                const auto& completion = completionSet.get(c);
                int i = completion.getI();
                int j = completion.getJ();
                int k = completion.getK();
                int cart_grid_idx = i + cartesianSize[0]*(j + cartesianSize[1]*k);
                int compressed_idx = cartesianToCompressed[cart_grid_idx];

                if ( compressed_idx >= 0 ) // Ignore completions in inactive/remote cells.
                {
                    compressed_well_perforations.push_back(compressed_idx);
                }
            }

            if ( ! compressed_well_perforations.empty() )
            {
                std::sort(compressed_well_perforations.begin(),
                          compressed_well_perforations.end());

                wells_.push_back(compressed_well_perforations);
            }
        }
    }

    unsigned numDofs() const
    {
        // No extra dofs are inserted for wells.
        return 0;
    }

    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        for(const auto well_perforations : wells_)
        {
            for(const auto& perforation : well_perforations)
                neighbors[perforation].insert(well_perforations.begin(),
                                              well_perforations.end());
        }
    }

    void applyInitial()
    {}

    void linearize(SparseMatrixAdapter& , GlobalEqVector&)
    {
        // Linearization is done in StandardDenseWells
    }

    private:

    std::vector<std::vector<int> > wells_;
};

} // end namespace OPM
#endif
