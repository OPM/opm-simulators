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

#include <ewoms/aux/baseauxiliarymodule.hh>

#include <opm/core/wells.h>

namespace Opm
{
template<class TypeTag>
class WellConnectionAuxiliaryModule
    : public Ewoms::BaseAuxiliaryModule<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;

public:

    using NeighborSet = typename
        Ewoms::BaseAuxiliaryModule<TypeTag>::NeighborSet;

    WellConnectionAuxiliaryModule(const Wells* wells)
        : wells_(wells)
    {
    }

    unsigned numDofs() const
    {
        // No extra dofs are inserted for wells.
        return 0;
    }

    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        const int nw = wells().number_of_wells;

        for (int w = 0; w < nw; ++w)
        {
            const int nperfs = wells().well_connpos[w+1];
            for (int perf = wells().well_connpos[w] ; perf < nperfs; ++perf) {
                const auto cell1_idx = wells().well_cells[perf];
                for(int perf1 = perf; perf1 < nperfs; ++perf1)
                {
                    const auto cell2_idx = wells().well_cells[perf1];
                    neighbors[cell1_idx].insert(cell2_idx);
                    neighbors[cell2_idx].insert(cell1_idx);
                }
            }
        }
    }

    void applyInitial()
    {}

    void linearize(JacobianMatrix& matrix, GlobalEqVector& residual)
    {
        // Linearization is done in StandardDenseWells
    }

    private:

    const Wells& wells() const
    {
        return *wells_;
    }

    const Wells* wells_;
};

} // end namespace OPM
#endif
