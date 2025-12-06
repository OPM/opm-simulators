// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef TPSA_FACE_PROPERTIES_IMPL_HPP
#define TPSA_FACE_PROPERTIES_IMPL_HPP

#ifndef TPSA_FACE_PROPERTIES_HPP
#include <config.h>
#include <opm/simulators/flow/FacePropertiesTPSA.hpp>
#endif

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/common/utility/ThreadSafeMapBuilder.hpp>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/models/parallel/threadmanager.hpp>

#include <algorithm>
#include <stdexcept>
#include <utility>


namespace Opm {

// Copied from Opm::Transmissibility class
namespace details {
    constexpr unsigned elemIdxShift = 32; // bits

    std::uint64_t isIdTPSA(std::uint32_t elemIdx1, std::uint32_t elemIdx2)
    {
        const std::uint32_t elemAIdx = std::min(elemIdx1, elemIdx2);
        const std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx << elemIdxShift) + elemAIdx;
    }
}  // namespace Opm::details

// /////
// Public functions
// ////
/*!
* \brief Constructor
*
* \param eclState Eclipse state made from a deck
* \param gridView The view on the DUNE grid which ought to be used (normally the leaf grid view)
* \param cartMapper The cartesian index mapper for lookup of cartesian indices
* \param grid The grid to lookup cell properties
* \param centroids Function to lookup cell centroids
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
FacePropertiesTPSA(const EclipseState& eclState,
                   const GridView& gridView,
                   const CartesianIndexMapper& cartMapper,
                   const Grid& grid,
                   std::function<std::array<double,dimWorld>(int)> centroids)
    : eclState_(eclState)
    , gridView_(gridView)
    , cartMapper_(cartMapper)
    , grid_(grid)
    , centroids_(centroids)
    , lookUpData_(gridView)
    , lookUpCartesianData_(gridView, cartMapper)
{ }

/*!
* \brief Compute TPSA face properties
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
finishInit()
{
    update();
}

/*!
* \brief Compute TPSA face properties
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
update()
{
    // Number of elements
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());
    unsigned numElements = elemMapper.size();

    // Extract shear modulus (Lame's sec. params)
    extractSModulus_();

    // Init. containers
    // Note (from Transmissibility::update): Reserving some space in the hashmap upfront saves quite a bit of time
    // because resizes are costly for hashmaps and there would be quite a few of them if we would not have a rough idea
    // of how large the final map will be (the rough idea is a conforming Cartesian grid).
    const int num_threads = ThreadManager::maxThreads();
    weightsAvg_.clear();
    if (num_threads == 1) {
        weightsAvg_.reserve(numElements * 3 * 1.05);
        weightsProd_.reserve(numElements * 3 * 1.05);
        distance_.reserve(numElements * 3 * 1.05);
        faceNormal_.reserve(numElements * 3 * 1.05);
    }
    weightsAvgBoundary_.clear();
    weightsProdBoundary_.clear();
    distanceBoundary_.clear();
    faceNormalBoundary_.clear();

    // Initialize thread safe insert_or_assign for face properties in the grid and separate for boundaries
    ThreadSafeMapBuilder weightsAvgMap(weightsAvg_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder weightsProdMap(weightsProd_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder distanceMap(distance_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder faceNormalMap(faceNormal_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);

    ThreadSafeMapBuilder weightsAvgBoundaryMap(weightsAvgBoundary_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder weightsProdBoundaryMap(weightsProdBoundary_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder distanceBoundaryMap(distanceBoundary_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder faceNormalBoundaryMap(faceNormalBoundary_, num_threads, MapBuilderInsertionMode::Insert_Or_Assign);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Loop over grid element an compute face properties
    for (const auto& chunk : ElementChunks(gridView_, Dune::Partitions::all, num_threads)) {
        for (const auto& elem : chunk) {
            // Init. face info for inside/outside cells
            FaceInfo inside;
            FaceInfo outside;
            DimVector faceNormal;

            // Set inside info
            inside.elemIdx = elemMapper.index(elem);
            inside.cartElemIdx = lookUpCartesianData_.
                template getFieldPropCartesianIdx<Grid>(inside.elemIdx);

            // Loop over intersection to neighboring cells
            unsigned boundaryIsIdx = 0;
            for (const auto& intersection : intersections(gridView_, elem)) {
                // Handle grid boundary
                if (intersection.boundary()) {
                    // One-sided cell properties
                    const auto& geometry = intersection.geometry();
                    inside.faceCenter = geometry.center();
                    faceNormal = intersection.centerUnitOuterNormal();

                    // Face properties on boundary
                    const auto index_pair = std::make_pair(inside.elemIdx, boundaryIsIdx);
                    Scalar distBound = computeDistance_(distanceVector_(inside.faceCenter, inside.elemIdx), faceNormal);
                    distanceBoundaryMap.insert_or_assign(index_pair, distBound);

                    Scalar weightsAvgBound = 1.0;  // w_j = 0 -> w_avg = wi / (wi + 0)
                    Scalar weightsProdBound = 0.0;  // w_j = 0 -> w_prod = wi * 0
                    weightsAvgBoundaryMap.insert_or_assign(index_pair, weightsAvgBound);
                    weightsProdBoundaryMap.insert_or_assign(index_pair, weightsProdBound);

                    faceNormalBoundaryMap.insert_or_assign(index_pair, faceNormal);

                    ++boundaryIsIdx;
                    continue;
                }

                // Handle intersection on process boundary (i.e., neighbor on different rank)
                if (!intersection.neighbor()) {
                    ++boundaryIsIdx;
                    continue;
                }

                // Set outside info
                const auto& outsideElem = intersection.outside();
                outside.elemIdx = elemMapper.index(outsideElem);
                outside.cartElemIdx = lookUpCartesianData_.
                    template getFieldPropCartesianIdx<Grid>(outside.elemIdx);

                // Skip intersections that have already been processed
                if (std::tie(inside.cartElemIdx, inside.elemIdx) > std::tie(outside.cartElemIdx, outside.elemIdx)) {
                    continue;
                }

                // Face indices for this intersection
                inside.faceIdx  = intersection.indexInInside();
                outside.faceIdx = intersection.indexInOutside();

                // Set NNC face properties to zero
                if (inside.faceIdx == -1) {
                    const auto id = details::isIdTPSA(inside.elemIdx, outside.elemIdx);
                    weightsAvgMap.insert_or_assign(id, 0.0);
                    distanceMap.insert_or_assign(id, 0.0);
                    faceNormalMap.insert_or_assign(id, DimVector{0.0, 0.0, 0.0});

                    continue;
                }

                // Compute cell  properties from grid
                typename std::is_same<Grid, Dune::CpGrid>::type isCpGrid;
                computeCellProperties(intersection,
                                      inside,
                                      outside,
                                      faceNormal,
                                      isCpGrid);

                // Compute face properties
                const auto id = details::isIdTPSA(inside.elemIdx, outside.elemIdx);
                Scalar dist_in = computeDistance_(distanceVector_(inside.faceCenter, inside.elemIdx), faceNormal);
                Scalar dist_out = computeDistance_(distanceVector_(outside.faceCenter, outside.elemIdx), faceNormal);
                distanceMap.insert_or_assign(id, dist_in + dist_out);

                Scalar weight_in = computeWeight_(dist_in, sModulus_[inside.elemIdx]);
                Scalar weight_out = computeWeight_(dist_out, sModulus_[outside.elemIdx]);
                Scalar weightsAvg = weight_in / (weight_in + weight_out);
                Scalar weightsProd = weight_in * weight_out;
                weightsAvgMap.insert_or_assign(id, weightsAvg);
                weightsProdMap.insert_or_assign(id, weightsProd);

                faceNormalMap.insert_or_assign(id, faceNormal);
            }
        }
    }
}

/*!
* \brief Average (half-)weight at interface between two elements
*
* \param elemIdx1 Cell index 1
* \param elemIdx2 Cell index 2
* \returns Average weight
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
weightAverage(unsigned elemIdx1, unsigned elemIdx2) const
{
    auto tmp_whgt = weightsAvg_.at(details::isIdTPSA(elemIdx1, elemIdx2));

    auto cartIdx1 = lookUpCartesianData_.
        template getFieldPropCartesianIdx<Grid>(elemIdx1);
    auto cartIdx2 = lookUpCartesianData_.
        template getFieldPropCartesianIdx<Grid>(elemIdx2);
    if (cartIdx1 < cartIdx2) {
        return tmp_whgt;
    }
    else {
        return 1.0 - tmp_whgt;
    }
}

/*!
* \brief Average (half-)weight at boundary interface
*
* \param elemIdx Cell index
* \param boundaryFaceIdx Face index
* \returns Average weight at boundary
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
weightAverageBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return weightsAvgBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

/*!
* \brief Product of weights at interface between two elements
*
* \param elemIdx1 Cell index 1
* \param elemIdx2 Cell index 2
* \returns Weight product
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
weightProduct(unsigned elemIdx1, unsigned elemIdx2) const
{
    return weightsProd_.at(details::isIdTPSA(elemIdx1, elemIdx2));
}

/*!
* \brief Product of weights at boundary interface
*
* \param elemIdx Cell index
* \param boundaryFaceIdx Face index
* \returns Weight product at boundary
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
weightProductBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return weightsProdBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

/*!
* \brief Distance between two elements
*
* \param elemIdx1 Cell index 1
* \param elemIdx2 Cell index 2
* \returns Distance
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
normalDistance(unsigned elemIdx1, unsigned elemIdx2) const
{
    return distance_.at(details::isIdTPSA(elemIdx1, elemIdx2));
}

/*!
* \brief Distance to boundary interface
*
* \param elemIdx Cell index
* \param boundaryFaceIdx Face index
* \returns Distance to boundary
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
normalDistanceBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return distanceBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

/*!
* \brief Cell face normal at interface between two elements
*
* \param elemIdx1 Cell index 1
* \param elemIdx2 Cell index 2
* \returns Face normals
*
* \note It is assumed that normals point from lower index to higher index!
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
typename FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::DimVector
FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
cellFaceNormal(unsigned elemIdx1, unsigned elemIdx2)
{
    auto cartIdx1 = lookUpCartesianData_.
        template getFieldPropCartesianIdx<Grid>(elemIdx1);
    auto cartIdx2 = lookUpCartesianData_.
        template getFieldPropCartesianIdx<Grid>(elemIdx2);
    int sign = (cartIdx1 < cartIdx2) ? 1 : -1;
    return sign * faceNormal_.at(details::isIdTPSA(elemIdx1, elemIdx2));
}

/*!
* \brief Cell face normal of boundary interface
*
* \param elemIdx Cell index
* \param boundaryFaceIdx Face index
* \returns Face normals
*
* \note Boundary normals points outwards
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
const typename FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::DimVector&
FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
cellFaceNormalBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return faceNormalBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

// /////
// Protected functions
// ////
/*!
* \brief Compute normal distance from cell center to face center
*
* \param distVec Distance vector from cell center to face center
* \param faceNormal Face (unit) normal vector
* \returns Distance cell center -> face center
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
computeDistance_(const DimVector& distVec, const DimVector& faceNormal)
{
    return std::abs(Dune::dot(faceNormal, distVec));
}

/*!
* \brief Compute face properties from general DUNE grid
*
* \param intersection Info on cell intersection
* \param inside Info describing inside face
* \param outside Info describing outside face
* \param faceNormal Face (unit) normal vector
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
computeCellProperties(const Intersection& intersection,
                      FaceInfo& inside,
                      FaceInfo& outside,
                      DimVector& faceNormal,
                      /*isCpGrid=*/std::false_type) const
{
    // default implementation for DUNE grids
    const auto& geometry = intersection.geometry();
    outside.faceCenter = inside.faceCenter = geometry.center();

    // OBS: Have not checked if this points from cell with lower to higher index!
    faceNormal = intersection.centerUnitOuterNormal();
}

/*!
* \brief Compute face properties from DUNE CpGrid
*
* \param intersection Info on cell intersection
* \param inside Info describing inside face
* \param outside Info describing outside face
* \param faceNormal Face (unit) normal
*
* \warning LGR computations not implemented!
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
computeCellProperties(const Intersection& intersection,
                      FaceInfo& inside,
                      FaceInfo& outside,
                      DimVector& faceNormal,
                      /*isCpGrid=*/std::true_type) const
{
    int faceIdx = intersection.id();
    if (grid_.maxLevel() == 0) {
        // Face center coordinates
        inside.faceCenter = grid_.faceCenterEcl(inside.elemIdx, inside.faceIdx, intersection);
        outside.faceCenter = grid_.faceCenterEcl(outside.elemIdx, outside.faceIdx, intersection);

        // Face normal, ensuring it points from cell with lower to higher (global grid) index
        faceNormal = grid_.faceNormal(faceIdx);
        auto cartFaceCell0 = lookUpCartesianData_.
            template getFieldPropCartesianIdx<Grid>(grid_.faceCell(faceIdx, 0));
        auto cartFaceCell1 = lookUpCartesianData_.
            template getFieldPropCartesianIdx<Grid>(grid_.faceCell(faceIdx, 1));
        if (cartFaceCell0 > cartFaceCell1) {
            faceNormal *= -1;
        }
    }
    else {
        throw std::runtime_error("TPSA not implemented with LGR");
    }

}

/*!
* \brief Compute weight ratio between distance and shear modulus
*
* \param distance Normal distance from cell to face center
* \param smod Shear modulus
* \returns Weight ratio
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
computeWeight_(const Scalar distance, const Scalar smod)
{
    return distance / smod;
}

/*!
* \brief Distance vector from cell center to face center
*
* \param faceCenter Face center coordinates
* \param cellIdx Cell index
* \returns Distance vector cell center -> face center
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
typename FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::DimVector
FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
distanceVector_(const DimVector& faceCenter, const unsigned& cellIdx) const
{
    const auto& cellCenter = centroids_cache_.empty() ? centroids_(cellIdx)
                                                      : centroids_cache_[cellIdx];
    DimVector x = faceCenter;
    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
        x[dimIdx] -= cellCenter[dimIdx];
    }

    return x;
}

/*!
* \brief Extract shear modulus from eclState
*
* \note (from Transmissibility::extractPorosity()):
*   All arrays provided by eclState are one-per-cell of "uncompressed" grid, whereas the simulation grid might remove a
*   few elements (e.g. because it is distributed over several processes).
*/
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>::
extractSModulus_()
{
    unsigned numElem = gridView_.size(/*codim=*/0);
    sModulus_.resize(numElem);

    const auto& fp = eclState_.fieldProps();
    std::vector<double> sModulusData;
    if (fp.has_double("SMODULUS")) {
        sModulusData = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "SMODULUS");
    }
    else if (fp.has_double("YMODULE") && fp.has_double("LAME")) {
        // Convert from Young's modulus and Lame's first parameter
        const std::vector<double>& ymodulus = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "YMODULE");
        const std::vector<double>& lameParam = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "LAME");
        sModulusData.resize(numElem);
        for (std::size_t i = 0; i < sModulusData.size(); ++i) {
            const double r = std::sqrt(ymodulus[i] * ymodulus[i] + 9 * lameParam[i] * lameParam[i]
                                       + 2 * ymodulus[i] * lameParam[i]);
            sModulusData[i] = (ymodulus[i] - 3 * lameParam[i] + r) / 4.0;
        }
    }
    else if (fp.has_double("YMODULE") && fp.has_double("PRATIO")) {
        const std::vector<double>& ymodulus = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "YMODULE");
        const std::vector<double>& pratio = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "PRATIO");
        sModulusData.resize(numElem);
        for (std::size_t i = 0; i < sModulusData.size(); ++i) {
            sModulusData[i] = ymodulus[i] / (2 * (1 + pratio[i]));
        }
    }
    else if (fp.has_double("LAME") && fp.has_double("PRATIO")) {
        const std::vector<double>& lameParam = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "LAME");
        const std::vector<double>& pratio = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "PRATIO");
        sModulusData.resize(numElem);
        for (std::size_t i = 0; i < sModulusData.size(); ++i) {
            sModulusData[i] = lameParam[i] * (1 - 2 * pratio[i]) / (2 * pratio[i]);
        }
    }
    else {
        throw std::logic_error("Cannot read shear modulus data from ecl state, SMODULUS keyword missing, "
                               "and one of the following keyword pairs are missing for conversion: "
                               "(YMODULE, LAME), (YMODULE, PRATIO) and (LAME, PRATIO)!");
    }

    // Assign shear modulus
    for (std::size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
        sModulus_[dofIdx] = sModulusData[dofIdx];
    }
}

}  // namespace Opm

#endif
