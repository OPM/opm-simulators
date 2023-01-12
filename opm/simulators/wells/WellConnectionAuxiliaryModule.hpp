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

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <vector>

namespace Dune { class CpGrid; }

namespace Opm
{

class Schedule;

class WellConnectionAuxiliaryModuleGeneric
{
protected:
    WellConnectionAuxiliaryModuleGeneric(const Schedule& schedule,
                                         const Dune::CpGrid& grid);

    std::vector<std::vector<int> > wells_;
};

template<class TypeTag>
class WellConnectionAuxiliaryModule
    : public BaseAuxiliaryModule<TypeTag>
    , private WellConnectionAuxiliaryModuleGeneric
{
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

public:

    using NeighborSet = typename
        ::Opm::BaseAuxiliaryModule<TypeTag>::NeighborSet;

    WellConnectionAuxiliaryModule(const Schedule& schedule,
                                  const Dune::CpGrid& grid)
        : WellConnectionAuxiliaryModuleGeneric(schedule, grid)
    {
    }

    unsigned numDofs() const
    {
        // No extra dofs are inserted for wells.
        return 0;
    }

    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        for (const auto& well_perforations : wells_)
        {
            for (const auto& perforation : well_perforations)
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
};

} // end namespace OPM
#endif
