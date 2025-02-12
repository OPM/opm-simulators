/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BOUNDINGBOXTREE_INCLUDED
#define OPM_BOUNDINGBOXTREE_INCLUDED
#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractor.hpp>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>
#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>
#include <external/resinsight/CommonCode/cvfStructGrid.h>
#include <external/resinsight/LibGeometry/cvfBoundingBox.h>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>

#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractorGrid.hpp>
//#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractorGrid_impl.hpp>
namespace external
{
void
buildBoundingBoxTree(cvf::ref<cvf::BoundingBoxTree>& m_cellSearchTree, const ::Dune::CpGrid& grid)
{
    using Grid = ::Dune::CpGrid;
    using GridView = typename Grid::LeafGridView;
    const auto& gv = grid.leafGridView();
    size_t cellCount = gv.size(0);
    std::vector<size_t> cellIndicesForBoundingBoxes;
    std::vector<cvf::BoundingBox> cellBoundingBoxes;

    std::array<double, 3> cornerPointArray;
    cvf::Vec3d cornerPoint;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interally
    for (const auto& element : Dune::elements(gv)) {
        int index = mapper.index(element);
        auto geom = element.geometry();
        cvf::BoundingBox cellBB;
        cvf::Vec3d cornerPoint;
        // NB order should not matter when adding to bounding box: dune ordring and resinsight ordering is different
        //  dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
        for (std::size_t l = 0; l < geom.corners(); l++) {
            auto cornerPointArray = geom.corner(l);
            cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
            cellBB.add(cornerPoint);
        }
        cellIndicesForBoundingBoxes.emplace_back(index);
        cellBoundingBoxes.emplace_back(cellBB);
    }
    m_cellSearchTree = new cvf::BoundingBoxTree;
    m_cellSearchTree->buildTreeFromBoundingBoxes(cellBoundingBoxes, &cellIndicesForBoundingBoxes);
}
}
#endif // OPM_BLACKOILWELLMODEL_WBP_HEADER_INCLUDED