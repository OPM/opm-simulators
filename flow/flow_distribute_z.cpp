/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015, 2017 IRIS AS
  Copyright 2021 OPM-OP AS

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

#include <opm/simulators/flow/Main.hpp>
#include <ebos/eclcpgridvanguard.hh>

std::vector<int> loadBalanceInZOnly(const Dune::CpGrid& grid)
{
    auto cartMapper = Dune::CartesianIndexMapper<Dune::CpGrid>(grid);
    auto dims = cartMapper.cartesianDimensions();
    std::vector<int> parts(grid.leafGridView().size(0));
    auto numCellsPerProc = dims[2]/grid.comm().size();

    if (grid.size(0)>0 && numCellsPerProc == 0)
    {
        OPM_THROW(std::logic_error,
                  "cartesian grid must have more cells in z direction than number of processes.");
    }

    using ElementMapper =
        Dune::MultipleCodimMultipleGeomTypeMapper<typename Dune::CpGrid::LeafGridView>;
    const auto& gridView = grid.leafGridView();
    const auto& idSet = grid.localIdSet();
    ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

    for( const auto &element : elements(gridView) )
    {
        const auto& id = idSet.id(element);
        unsigned elemIdx = elemMapper.index(element);
        const auto& cartIndex = cartMapper.cartesianIndex(elemIdx);
        std::array<int,3> cartCoord;
        cartMapper.cartesianCoordinate(cartIndex, cartCoord);
        using std::min;
        auto rank = min(grid.comm().size() -1, cartCoord[2] / numCellsPerProc);
        parts[id] = rank;
    }
    return parts;
}

int main(int argc, char** argv)
{
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    Opm::EclCpGridVanguard<Opm::Properties::TTag::FlowProblem>::setExternalLoadBalancer(loadBalanceInZOnly);
    auto ret = mainObject->runDynamic();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

