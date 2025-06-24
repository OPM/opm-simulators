// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
#ifndef OPM_TRANSMISSIBILITY_IMPL_HPP
#define OPM_TRANSMISSIBILITY_IMPL_HPP

#ifndef OPM_TRANSMISSIBILITY_HPP
#include <config.h>
#include <opm/simulators/flow/Transmissibility.hpp>
#endif

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/common/utility/ThreadSafeMapBuilder.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/models/parallel/threadmanager.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace Opm {

namespace details {

    constexpr unsigned elemIdxShift = 32; // bits

    std::uint64_t isId(std::uint32_t elemIdx1, std::uint32_t elemIdx2)
    {
        const std::uint32_t elemAIdx = std::min(elemIdx1, elemIdx2);
        const std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx << elemIdxShift) + elemAIdx;
    }

    std::pair<std::uint32_t, std::uint32_t> isIdReverse(const std::uint64_t& id)
    {
        // Assigning an unsigned integer to a narrower type discards the most significant bits.
        // See "The C programming language", section A.6.2.
        // NOTE that the ordering of element A and B may have changed
        const std::uint32_t elemAIdx = static_cast<uint32_t>(id);
        const std::uint32_t elemBIdx = (id - elemAIdx) >> elemIdxShift;

        return std::make_pair(elemAIdx, elemBIdx);
    }

    std::uint64_t directionalIsId(std::uint32_t elemIdx1, std::uint32_t elemIdx2)
    {
        return (std::uint64_t(elemIdx1) << elemIdxShift) + elemIdx2;
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
Transmissibility(const EclipseState& eclState,
                 const GridView& gridView,
                 const CartesianIndexMapper& cartMapper,
                 const Grid& grid,
                 std::function<std::array<double,dimWorld>(int)> centroids,
                 bool enableEnergy,
                 bool enableDiffusivity,
                 bool enableDispersivity)
      : eclState_(eclState)
      , gridView_(gridView)
      , cartMapper_(cartMapper)
      , grid_(grid)
      , centroids_(centroids)
      , enableEnergy_(enableEnergy)
      , enableDiffusivity_(enableDiffusivity)
      , enableDispersivity_(enableDispersivity)
      , lookUpData_(gridView)
      , lookUpCartesianData_(gridView, cartMapper)
{
    const UnitSystem& unitSystem = eclState_.getDeckUnitSystem();
    transmissibilityThreshold_  = unitSystem.parse("Transmissibility").getSIScaling() * 1e-6;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
{
    return trans_.at(details::isId(elemIdx1, elemIdx2));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
transmissibilityBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return transBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
thermalHalfTrans(unsigned insideElemIdx, unsigned outsideElemIdx) const
{
    return thermalHalfTrans_.at(details::directionalIsId(insideElemIdx, outsideElemIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const
{
    return thermalHalfTransBoundary_.at(std::make_pair(insideElemIdx, boundaryFaceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
diffusivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (diffusivity_.empty())
        return 0.0;

    return diffusivity_.at(details::isId(elemIdx1, elemIdx2));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
dispersivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (dispersivity_.empty())
        return 0.0;

    return dispersivity_.at(details::isId(elemIdx1, elemIdx2));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
update(bool global, const TransUpdateQuantities update_quantities,
       const std::function<unsigned int(unsigned int)>& map, const bool applyNncMultregT)
{
    // whether only update the permeability related transmissibility
    const bool onlyTrans = (update_quantities == TransUpdateQuantities::Trans);
    const auto& cartDims = cartMapper_.cartesianDimensions();
    const auto& transMult = eclState_.getTransMult();
    const auto& comm = gridView_.comm();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    unsigned numElements = elemMapper.size();
    // get the ntg values, the ntg values are modified for the cells merged with minpv
    const std::vector<double>& ntg = this->lookUpData_.assignFieldPropsDoubleOnLeaf(eclState_.fieldProps(), "NTG");
    const bool updateDiffusivity = eclState_.getSimulationConfig().isDiffusive();
    const bool updateDispersivity = eclState_.getSimulationConfig().rock_config().dispersion();

    const bool disableNNC = eclState_.getSimulationConfig().useNONNC();

    if (map) {
        extractPermeability_(map);
    }
    else {
        extractPermeability_();
    }

    const int num_threads = ThreadManager::maxThreads();

    // reserving some space in the hashmap upfront saves quite a bit of time because
    // resizes are costly for hashmaps and there would be quite a few of them if we
    // would not have a rough idea of how large the final map will be (the rough idea
    // is a conforming Cartesian grid).
    trans_.clear();
    if (num_threads == 1) {
        trans_.reserve(numElements*3*1.05);
    }

    transBoundary_.clear();

    // if energy is enabled, let's do the same for the "thermal half transmissibilities"
    if (enableEnergy_ && !onlyTrans) {
        thermalHalfTrans_.clear();
        if (num_threads == 1) {
            thermalHalfTrans_.reserve(numElements*6*1.05);
        }

        thermalHalfTransBoundary_.clear();
    }

    // if diffusion is enabled, let's do the same for the "diffusivity"
    if (updateDiffusivity && !onlyTrans) {
        diffusivity_.clear();
        if (num_threads == 1) {
            diffusivity_.reserve(numElements*3*1.05);
        }
        extractPorosity_();
    }

    // if dispersion is enabled, let's do the same for the "dispersivity"
    if (updateDispersivity && !onlyTrans) {
        dispersivity_.clear();
        if (num_threads == 1) {
            dispersivity_.reserve(numElements*3*1.05);
        }
        extractDispersion_();
    }

    // The MULTZ needs special case if the option is ALL
    // Then the smallest multiplier is applied.
    // Default is to apply the top and bottom multiplier
    bool useSmallestMultiplier;
    bool pinchOption4ALL;
    bool pinchActive;
    if (comm.rank() == 0) {
        const auto& eclGrid = eclState_.getInputGrid();
        pinchActive = eclGrid.isPinchActive();
        auto pinchTransCalcMode = eclGrid.getPinchOption();
        useSmallestMultiplier = eclGrid.getMultzOption() == PinchMode::ALL;
        pinchOption4ALL = (pinchTransCalcMode == PinchMode::ALL);
        if (pinchOption4ALL) {
            useSmallestMultiplier = false;
        }
    }
    if (global && comm.size() > 1) {
        comm.broadcast(&useSmallestMultiplier, 1, 0);
        comm.broadcast(&pinchOption4ALL, 1, 0);
        comm.broadcast(&pinchActive, 1, 0);
    }

    // fill the centroids cache to avoid repeated calculations in loops below
    centroids_cache_.resize(gridView_.size(0));
    for (const auto& elem : elements(gridView_)) {
        const unsigned elemIdx = elemMapper.index(elem);
        centroids_cache_[elemIdx] = centroids_(elemIdx);
    }

    auto harmonicMean = [](const Scalar x1, const Scalar x2)
    {
        return (std::abs(x1) < 1e-30 || std::abs(x2) < 1e-30)
            ? 0.0
            : 1.0 / (1.0 / x1 + 1.0 / x2);
    };

    auto faceIdToDir = [](int insideFaceIdx)
    {
        switch (insideFaceIdx) {
        case 0:
        case 1:
            return FaceDir::XPlus;
        case 2:
        case 3:
            return FaceDir::YPlus;
            break;
        case 4:
        case 5:
            return FaceDir::ZPlus;
        default:
            throw std::logic_error("Could not determine a face direction");
        }
    };

    auto halfDiff = [](const DimVector& faceAreaNormal,
                       const unsigned,
                       const DimVector& distVector,
                       const Scalar prop)
    {
        return computeHalfDiffusivity_(faceAreaNormal,
                                       distVector,
                                       prop);
    };

    ThreadSafeMapBuilder transBoundary(transBoundary_, num_threads,
                                       MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder transMap(trans_, num_threads,
                                  MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder thermalHalfTransBoundary(thermalHalfTransBoundary_, num_threads,
                                                  MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder thermalHalfTrans(thermalHalfTrans_, num_threads,
                                          MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder diffusivity(diffusivity_, num_threads,
                                     MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder dispersivity(dispersivity_, num_threads,
                                      MapBuilderInsertionMode::Insert_Or_Assign);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (const auto& chunk : ElementChunks(gridView_, Dune::Partitions::all, num_threads)) {
        for (const auto& elem : chunk) {
            FaceInfo inside;
            FaceInfo outside;
            DimVector faceAreaNormal;

            inside.elemIdx = elemMapper.index(elem);
            // Get the Cartesian index of the origin cells (parent or equivalent cell on level zero),
            // for CpGrid with LGRs. For general grids and no LGRs, get the usual Cartesian Index.
            inside.cartElemIdx = this->lookUpCartesianData_.
                template getFieldPropCartesianIdx<Grid>(inside.elemIdx);

            auto computeHalf = [this, &faceAreaNormal, &inside, &outside]
                               (const auto& halfComputer,
                                const auto& prop1, const auto& prop2) -> std::array<Scalar,2>
            {
                return {
                    halfComputer(faceAreaNormal,
                                 inside.faceIdx,
                                 distanceVector_(inside.faceCenter, inside.elemIdx),
                                 prop1),
                    halfComputer(faceAreaNormal,
                                 outside.faceIdx,
                                 distanceVector_(outside.faceCenter, outside.elemIdx),
                                 prop2)
                };
            };

            auto computeHalfMean = [&inside, &outside, &computeHalf, &ntg, &harmonicMean]
                                   (const auto& halfComputer, const auto& prop)
            {
                auto onesided = computeHalf(halfComputer, prop[inside.elemIdx], prop[outside.elemIdx]);
                applyNtg_(onesided[0], inside, ntg);
                applyNtg_(onesided[1], outside, ntg);

                //TODO Add support for multipliers
                return harmonicMean(onesided[0], onesided[1]);
            };

            unsigned boundaryIsIdx = 0;
            for (const auto& intersection : intersections(gridView_, elem)) {
                // deal with grid boundaries
                if (intersection.boundary()) {
                    // compute the transmissibilty for the boundary intersection
                    const auto& geometry = intersection.geometry();
                    inside.faceCenter = geometry.center();

                    faceAreaNormal = intersection.centerUnitOuterNormal();
                    faceAreaNormal *= geometry.volume();

                    Scalar transBoundaryIs =
                        computeHalfTrans_(faceAreaNormal,
                                          intersection.indexInInside(),
                                          distanceVector_(inside.faceCenter, inside.elemIdx),
                                          permeability_[inside.elemIdx]);

                    // normally there would be two half-transmissibilities that would be
                    // averaged. on the grid boundary there only is the half
                    // transmissibility of the interior element.
                    applyMultipliers_(transBoundaryIs, intersection.indexInInside(), inside.cartElemIdx, transMult);
                    transBoundary.insert_or_assign(std::make_pair(inside.elemIdx, boundaryIsIdx), transBoundaryIs);

                    // for boundary intersections we also need to compute the thermal
                    // half transmissibilities
                    if (enableEnergy_ && !onlyTrans) {
                        Scalar transBoundaryEnergyIs =
                            computeHalfDiffusivity_(faceAreaNormal,
                                                    distanceVector_(inside.faceCenter, inside.elemIdx),
                                                    1.0);
                        thermalHalfTransBoundary.insert_or_assign(std::make_pair(inside.elemIdx, boundaryIsIdx),
                                                                   transBoundaryEnergyIs);
                    }

                    ++boundaryIsIdx;
                    continue;
                }

                if (!intersection.neighbor()) {
                    // elements can be on process boundaries, i.e. they are not on the
                    // domain boundary yet they don't have neighbors.
                    ++boundaryIsIdx;
                    continue;
                }

                const auto& outsideElem = intersection.outside();
                outside.elemIdx = elemMapper.index(outsideElem);

                // Get the Cartesian index of the origin cells (parent or equivalent cell on level zero),
                // for CpGrid with LGRs. For general grids and no LGRs, get the usual Cartesian Index.
                outside.cartElemIdx =  this->lookUpCartesianData_.
                    template getFieldPropCartesianIdx<Grid>(outside.elemIdx);

                // we only need to calculate a face's transmissibility
                // once...
                // In a parallel run inside.cartElemIdx > outside.cartElemIdx does not imply inside.elemIdx > outside.elemIdx for
                // ghost cells and we need to use the cartesian index as this will be used when applying Z multipliers
                // To cover the case where both cells are part of an LGR and as a consequence might have
                // the same cartesian index, we tie their Cartesian indices and the ones on the leaf grid view.
                if (std::tie(inside.cartElemIdx, inside.elemIdx) > std::tie(outside.cartElemIdx, outside.elemIdx)) {
                    continue;
                }

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                inside.faceIdx  = intersection.indexInInside();
                outside.faceIdx = intersection.indexInOutside();

                if (inside.faceIdx == -1) {
                    // NNC. Set zero transmissibility, as it will be
                    // *added to* by applyNncToGridTrans_() later.
                    assert(outside.faceIdx == -1);
                    transMap.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx), 0.0);
                    if (enableEnergy_ && !onlyTrans) {
                        thermalHalfTrans.insert_or_assign(details::directionalIsId(inside.elemIdx, outside.elemIdx), 0.0);
                        thermalHalfTrans.insert_or_assign(details::directionalIsId(outside.elemIdx, inside.elemIdx), 0.0);
                    }

                    if (updateDiffusivity && !onlyTrans) {
                        diffusivity.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx), 0.0);
                    }
                    if (updateDispersivity && !onlyTrans) {
                        dispersivity.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx), 0.0);
                    }
                    continue;
                }

                typename std::is_same<Grid, Dune::CpGrid>::type isCpGrid;
                computeFaceProperties(intersection,
                                      inside,
                                      outside,
                                      faceAreaNormal,
                                      isCpGrid);

                Scalar trans = computeHalfMean(computeHalfTrans_, permeability_);

                // apply the full face transmissibility multipliers
                // for the inside ...
                if (!pinchActive) {
                    if (inside.faceIdx > 3) { // top or bottom
                         auto find_layer = [&cartDims](std::size_t cell) {
                            cell /= cartDims[0];
                            auto k = cell / cartDims[1];
                            return k;
                        };
                        int kup = find_layer(inside.cartElemIdx);
                        int kdown = find_layer(outside.cartElemIdx);
                        // When a grid is a CpGrid with LGRs, insideCartElemIdx coincides with outsideCartElemIdx
                        // for cells on the leaf with the same parent cell on level zero.
                        assert((kup != kdown) || (inside.cartElemIdx == outside.cartElemIdx));
                        if (std::abs(kup -kdown) > 1) {
                            trans = 0.0;
                        }
                    }
                }

                if (useSmallestMultiplier) {
                    //  PINCH(4) == TOPBOT is assumed here as we set useSmallestMultipliers
                    // to false if  PINCH(4) == ALL holds
                    // In contrast to the name this will also apply
                    applyAllZMultipliers_(trans, inside, outside, transMult, cartDims);
                }
                else {
                    applyMultipliers_(trans, inside.faceIdx, inside.cartElemIdx, transMult);
                    // ... and outside elements
                    applyMultipliers_(trans, outside.faceIdx, outside.cartElemIdx, transMult);
                }

                // apply the region multipliers (cf. the MULTREGT keyword)
                trans *= transMult.getRegionMultiplier(inside.cartElemIdx,
                                                       outside.cartElemIdx,
                                                       faceIdToDir(inside.faceIdx));

                transMap.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx), trans);

                // update the "thermal half transmissibility" for the intersection
                if (enableEnergy_ && !onlyTrans) {
                    const auto half = computeHalf(halfDiff, 1.0, 1.0);
                    // TODO Add support for multipliers
                    thermalHalfTrans.insert_or_assign(details::directionalIsId(inside.elemIdx, outside.elemIdx),
                                                      half[0]);
                    thermalHalfTrans.insert_or_assign(details::directionalIsId(outside.elemIdx, inside.elemIdx),
                                                      half[1]);
                }

                // update the "diffusive half transmissibility" for the intersection
                if (updateDiffusivity && !onlyTrans) {
                    diffusivity.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx),
                                                 computeHalfMean(halfDiff, porosity_));
                }

                // update the "dispersivity half transmissibility" for the intersection
                if (updateDispersivity && !onlyTrans) {
                    dispersivity.insert_or_assign(details::isId(inside.elemIdx, outside.elemIdx),
                                                  computeHalfMean(halfDiff, dispersion_));
                }
            }
        }
    }
    centroids_cache_.clear();

#ifdef _OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef _OPENMP
#pragma omp section
#endif
        transMap.finalize();
#ifdef _OPENMP
#pragma omp section
#endif
        transBoundary.finalize();
#ifdef _OPENMP
#pragma omp section
#endif
        thermalHalfTransBoundary.finalize();
#ifdef _OPENMP
#pragma omp section
#endif
        thermalHalfTrans.finalize();
#ifdef _OPENMP
#pragma omp section
#endif
        diffusivity.finalize();
#ifdef _OPENMP
#pragma omp section
#endif
        dispersivity.finalize();
    }

    // Potentially overwrite and/or modify transmissibilities based on input from deck
    this->updateFromEclState_(global);

    // Create mapping from global to local index
    std::unordered_map<std::size_t,int> globalToLocal;

    // Loop over all elements (global grid) and store Cartesian index
    for (const auto& elem : elements(grid_.leafGridView())) {
        int elemIdx = elemMapper.index(elem);
        int cartElemIdx =  cartMapper_.cartesianIndex(elemIdx);
        globalToLocal[cartElemIdx] = elemIdx;
    }

    if (!disableNNC) {
        // For EDITNNC and EDITNNCR we warn only once
        // If transmissibility is used for load balancing this will be done
        // when computing the gobal transmissibilities and all warnings will
        // be seen in a parallel. Unfortunately, when we do not use transmissibilities
        // we will only see warnings for the partition of process 0 and also false positives.
        this->applyEditNncToGridTrans_(globalToLocal);
        this->applyPinchNncToGridTrans_(globalToLocal);
        this->applyNncToGridTrans_(globalToLocal);
        this->applyEditNncrToGridTrans_(globalToLocal);
        if (applyNncMultregT) {
            this->applyNncMultreg_(globalToLocal);
        }
        warnEditNNC_ = false;
    }

    // If disableNNC == true, remove all non-neighbouring transmissibilities.
    // If disableNNC == false, remove very small non-neighbouring transmissibilities.
    this->removeNonCartesianTransmissibilities_(disableNNC);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractPermeability_()
{
    unsigned numElem = gridView_.size(/*codim=*/0);
    permeability_.resize(numElem);

    // read the intrinsic permeabilities from the eclState. Note that all arrays
    // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
    // simulation grid might remove a few elements. (e.g. because it is distributed
    // over several processes.)
    const auto& fp = eclState_.fieldProps();
    if (fp.has_double("PERMX")) {
        const std::vector<double>& permxData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "PERMX");

        std::vector<double> permyData;
        if (fp.has_double("PERMY"))
            permyData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMY");
        else
            permyData = permxData;

        std::vector<double> permzData;
        if (fp.has_double("PERMZ"))
            permzData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMZ");
        else
            permzData = permxData;

        for (std::size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
            permeability_[dofIdx] = 0.0;
            permeability_[dofIdx][0][0] = permxData[dofIdx];
            permeability_[dofIdx][1][1] = permyData[dofIdx];
            permeability_[dofIdx][2][2] = permzData[dofIdx];
        }

        // for now we don't care about non-diagonal entries

    }
    else
        throw std::logic_error("Can't read the intrinsic permeability from the ecl state. "
                               "(The PERM{X,Y,Z} keywords are missing)");
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractPermeability_(const std::function<unsigned int(unsigned int)>& map)
{
    unsigned numElem = gridView_.size(/*codim=*/0);
    permeability_.resize(numElem);

    // read the intrinsic permeabilities from the eclState. Note that all arrays
    // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
    // simulation grid might remove a few elements. (e.g. because it is distributed
    // over several processes.)
    const auto& fp = eclState_.fieldProps();
    if (fp.has_double("PERMX")) {
        const std::vector<double>& permxData =
            this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMX");

        std::vector<double> permyData;
        if (fp.has_double("PERMY")){
            permyData = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMY");
        }
        else {
            permyData = permxData;
        }

        std::vector<double> permzData;
        if (fp.has_double("PERMZ")) {
            permzData = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMZ");
        }
        else {
            permzData = permxData;
        }

        for (std::size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
            permeability_[dofIdx] = 0.0;
            std::size_t inputDofIdx = map(dofIdx);
            permeability_[dofIdx][0][0] = permxData[inputDofIdx];
            permeability_[dofIdx][1][1] = permyData[inputDofIdx];
            permeability_[dofIdx][2][2] = permzData[inputDofIdx];
        }

        // for now we don't care about non-diagonal entries
    }
    else {
        throw std::logic_error("Can't read the intrinsic permeability from the ecl state. "
                               "(The PERM{X,Y,Z} keywords are missing)");
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractPorosity_()
{
    // read the intrinsic porosity from the eclState. Note that all arrays
    // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
    // simulation grid might remove a few elements. (e.g. because it is distributed
    // over several processes.)
    const auto& fp = eclState_.fieldProps();
    if (fp.has_double("PORO")) {
        if constexpr (std::is_same_v<Scalar,double>) {
            porosity_ = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PORO");
        }
        else {
            const auto por = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PORO");
            porosity_.resize(por.size());
            std::copy(por.begin(), por.end(), porosity_.begin());
        }
    }
    else {
        throw std::logic_error("Can't read the porosity from the ecl state. "
                               "(The PORO keywords are missing)");
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractDispersion_()
{
    if (!enableDispersivity_) {
        throw std::runtime_error("Dispersion disabled at compile time, but the deck "
                                 "contains the DISPERC keyword.");
    }
    const auto& fp = eclState_.fieldProps();
    if constexpr (std::is_same_v<Scalar,double>) {
        dispersion_ = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"DISPERC");
    }
    else {
        const auto disp = this->lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"DISPERC");
        dispersion_.resize(disp.size());
        std::copy(disp.begin(), disp.end(), dispersion_.begin());
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
removeNonCartesianTransmissibilities_(bool removeAll)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    for (auto&& trans: trans_) {
        //either remove all NNC transmissibilities or those less than the threshold (by default 1e-6 in the deck's unit system)
        if (removeAll || trans.second < transmissibilityThreshold_) {
            const auto& id = trans.first;
            const auto& elements = details::isIdReverse(id);
            int gc1 = std::min(cartMapper_.cartesianIndex(elements.first), cartMapper_.cartesianIndex(elements.second));
            int gc2 = std::max(cartMapper_.cartesianIndex(elements.first), cartMapper_.cartesianIndex(elements.second));

            // only adjust the NNCs
            // When LGRs, all neighbors in the LGR are cartesian neighbours on the level grid representing the LGR.
            // When elements on the leaf grid view have the same parent cell, gc1 and gc2 coincide.
            if (gc2 - gc1 == 1 || gc2 - gc1 == cartDims[0] || gc2 - gc1 == cartDims[0]*cartDims[1] || gc2 - gc1 == 0) {
                continue;
            }

            trans.second = 0.0;
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper, Scalar>::
applyAllZMultipliers_(Scalar& trans,
                      const FaceInfo& inside,
                      const FaceInfo& outside,
                      const TransMult& transMult,
                      const std::array<int, dimWorld>& cartDims)
{
    if (grid_.maxLevel() > 0) {
        OPM_THROW(std::invalid_argument, "MULTZ not support with LGRS, yet.");
    }
    if (inside.faceIdx > 3) { // top or or bottom
        assert(inside.faceIdx == 5); // as insideCartElemIdx < outsideCartElemIdx holds for the Z column
        // For CpGrid with LGRs, insideCartElemIdx == outsideCartElemIdx when cells on the leaf have the same parent cell on level zero.
        assert(outside.cartElemIdx >= inside.cartElemIdx);
        unsigned lastCartElemIdx;
        if (outside.cartElemIdx == inside.cartElemIdx) {
            lastCartElemIdx = outside.cartElemIdx;
        }
        else {
            lastCartElemIdx = outside.cartElemIdx - cartDims[0]*cartDims[1];
        }
        // Last multiplier using (Z+)*(Z-)
        Scalar mult = transMult.getMultiplier(lastCartElemIdx , FaceDir::ZPlus) *
            transMult.getMultiplier(outside.cartElemIdx , FaceDir::ZMinus);

        // pick the smallest multiplier using (Z+)*(Z-) while looking down
        // the pillar until reaching the other end of the connection
        for (auto cartElemIdx = inside.cartElemIdx; cartElemIdx < lastCartElemIdx;) {
            auto multiplier = transMult.getMultiplier(cartElemIdx, FaceDir::ZPlus);
            cartElemIdx += cartDims[0]*cartDims[1];
            multiplier *= transMult.getMultiplier(cartElemIdx, FaceDir::ZMinus);
            mult = std::min(mult, static_cast<Scalar>(multiplier));
        }

        trans *= mult;
    }
    else {
        applyMultipliers_(trans, inside.faceIdx, inside.cartElemIdx, transMult);
        applyMultipliers_(trans, outside.faceIdx, outside.cartElemIdx, transMult);
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
updateFromEclState_(bool global)
{
    const FieldPropsManager* fp =
        (global) ? &(eclState_.fieldProps()) :
        &(eclState_.globalFieldProps());

    std::array<bool,3> is_tran {fp->tran_active("TRANX"),
                                fp->tran_active("TRANY"),
                                fp->tran_active("TRANZ")};

    if (!(is_tran[0] || is_tran[1] || is_tran[2])) {
        // Skip unneeded expensive traversals
        return;
    }

    std::array<std::string, 3> keywords {"TRANX", "TRANY", "TRANZ"};
    std::array<std::vector<double>,3> trans = createTransmissibilityArrays_(is_tran);
    auto key = keywords.begin();
    auto perform = is_tran.begin();

    for (auto it = trans.begin(); it != trans.end(); ++it, ++key, ++perform) {
        if (*perform) {
            if (grid_.maxLevel() > 0) {
                OPM_THROW(std::invalid_argument, "Calculations on TRANX/TRANY/TRANZ arrays are not support with LGRS, yet.");
            }
            fp->apply_tran(*key, *it);
        }
    }

    resetTransmissibilityFromArrays_(is_tran, trans);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::array<std::vector<double>,3>
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
createTransmissibilityArrays_(const std::array<bool,3>& is_tran)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    auto numElem = gridView_.size(/*codim=*/0);
    std::array<std::vector<double>,3> trans = {
          std::vector<double>(is_tran[0] ? numElem : 0, 0),
          std::vector<double>(is_tran[1] ? numElem : 0, 0),
          std::vector<double>(is_tran[2] ? numElem : 0, 0)
    };

    // compute the transmissibilities for all intersections
    for (const auto& elem : elements(gridView_)) {
        for (const auto& intersection : intersections(gridView_, elem)) {
            // store intersection, this might be costly
            if (!intersection.neighbor()) {
                continue; // intersection is on the domain boundary
            }

            // In the EclState TRANX[c1] is transmissibility in X+
            // direction. we only store transmissibilities in the +
            // direction. Same for Y and Z. Ordering of compressed (c1,c2) and cartesian index
            // (gc1, gc2) is coherent (c1 < c2 <=> gc1 < gc2) only in a serial run.
            // In a parallel run this only holds in the interior as elements in the
            // ghost overlap region might be ordered after the others. Hence we need
            // to use the cartesian index to select the compressed index where to store
            // the transmissibility value.
            // c1 < c2 <=> gc1 < gc2 is no longer true (even in serial) when the grid is a
            // CpGrid with LGRs. When cells c1 and c2 have the same parent
            // cell on level zero, then gc1 == gc2.
            unsigned c1 = elemMapper.index(intersection.inside());
            unsigned c2 = elemMapper.index(intersection.outside());
            int gc1 = cartMapper_.cartesianIndex(c1);
            int gc2 = cartMapper_.cartesianIndex(c2);
            if (std::tie(gc1, c1) > std::tie(gc2, c2)) {
                // we only need to handle each connection once, thank you.
                // We do this when gc1 is smaller than the other to find the
                // correct place to store in parallel when ghost/overlap elements
                // are ordered last
                continue;
            }

            auto isID = details::isId(c1, c2);

            // For CpGrid with LGRs, when leaf grid view cells with indices c1 and c2
            // have the same parent cell on level zero, then gc2 - gc1 == 0. In that case,
            // 'intersection.indexInSIde()' needed to be checked to determine the direction, i.e.
            // add in the if/else-if  'gc2 == gc1 && intersection.indexInInside() == ... '
            if ((gc2 - gc1 == 1 || (gc2 == gc1 && (intersection.indexInInside() == 0 || intersection.indexInInside() == 1)))
                && cartDims[0] > 1)
            {
                if (is_tran[0]) {
                    // set simulator internal transmissibilities to values from inputTranx
                     trans[0][c1] = trans_[isID];
                }
            }
            else if ((gc2 - gc1 == cartDims[0] || (gc2 == gc1 && (intersection.indexInInside() == 2 || intersection.indexInInside() == 3)))
                     && cartDims[1] > 1)
            {
                if (is_tran[1]) {
                    // set simulator internal transmissibilities to values from inputTrany
                     trans[1][c1] = trans_[isID];
                }
            }
            else if (gc2 - gc1 == cartDims[0]*cartDims[1] ||
                     (gc2 == gc1 && (intersection.indexInInside() == 4 || intersection.indexInInside() == 5)))
            {
                if (is_tran[2]) {
                    // set simulator internal transmissibilities to values from inputTranz
                    trans[2][c1] = trans_[isID];
                }
            }
            // else.. We don't support modification of NNC at the moment.
        }
    }

    return trans;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
resetTransmissibilityFromArrays_(const std::array<bool,3>& is_tran,
                                 const std::array<std::vector<double>,3>& trans)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    // compute the transmissibilities for all intersections
    for (const auto& elem : elements(gridView_)) {
        for (const auto& intersection : intersections(gridView_, elem)) {
            if (!intersection.neighbor()) {
                continue; // intersection is on the domain boundary
            }

            // In the EclState TRANX[c1] is transmissibility in X+
            // direction. we only store transmissibilities in the +
            // direction. Same for Y and Z. Ordering of compressed (c1,c2) and cartesian index
            // (gc1, gc2) is coherent (c1 < c2 <=> gc1 < gc2) only in a serial run.
            // In a parallel run this only holds in the interior as elements in the
            // ghost overlap region might be ordered after the others. Hence we need
            // to use the cartesian index to select the compressed index where to store
            // the transmissibility value.
            // c1 < c2 <=> gc1 < gc2 is no longer true (even in serial) when the grid is a
            // CpGrid with LGRs. When cells c1 and c2 have the same parent
            // cell on level zero, then gc1 == gc2.
            unsigned c1 = elemMapper.index(intersection.inside());
            unsigned c2 = elemMapper.index(intersection.outside());
            int gc1 = cartMapper_.cartesianIndex(c1);
            int gc2 = cartMapper_.cartesianIndex(c2);
            if (std::tie(gc1, c1) > std::tie(gc2, c2)) {
                // we only need to handle each connection once, thank you.
                // We do this when gc1 is smaller than the other to find the
                // correct place to read in parallel when ghost/overlap elements
                // are ordered last
                continue;
            }

            auto isID = details::isId(c1, c2);

            // For CpGrid with LGRs, when leaf grid view cells with indices c1 and c2
            // have the same parent cell on level zero, then gc2 - gc1 == 0. In that case,
            // 'intersection.indexInSIde()' needed to be checked to determine the direction, i.e.
            // add in the if/else-if  'gc2 == gc1 && intersection.indexInInside() == ... '
            if ((gc2 - gc1 == 1  || (gc2 == gc1 && (intersection.indexInInside() == 0 || intersection.indexInInside() == 1)))
                 && cartDims[0] > 1)
            {
                if (is_tran[0]) {
                    // set simulator internal transmissibilities to values from inputTranx
                    trans_[isID] = trans[0][c1];
                }
            }
            else if ((gc2 - gc1 == cartDims[0] || (gc2 == gc1 && (intersection.indexInInside() == 2|| intersection.indexInInside() == 3)))
                     && cartDims[1] > 1)
            {
                if (is_tran[1]) {
                    // set simulator internal transmissibilities to values from inputTrany
                    trans_[isID] = trans[1][c1];
                }
            }
            else if (gc2 - gc1 == cartDims[0]*cartDims[1] ||
                     (gc2 == gc1 && (intersection.indexInInside() == 4 || intersection.indexInInside() == 5)))
            {
                if (is_tran[2]) {
                    // set simulator internal transmissibilities to values from inputTranz
                    trans_[isID] = trans[2][c1];
                }
            }

            // else.. We don't support modification of NNC at the moment.
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeFaceProperties(const Intersection& intersection,
                      FaceInfo& inside,
                      FaceInfo& outside,
                      DimVector& faceAreaNormal,
                      /*isCpGrid=*/std::false_type) const
{
    // default implementation for DUNE grids
    const auto& geometry = intersection.geometry();
    outside.faceCenter = inside.faceCenter = geometry.center();

    faceAreaNormal = intersection.centerUnitOuterNormal();
    faceAreaNormal *= geometry.volume();
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeFaceProperties(const Intersection& intersection,
                      FaceInfo& inside,
                      FaceInfo& outside,
                      DimVector& faceAreaNormal,
                      /*isCpGrid=*/std::true_type) const
{
    int faceIdx = intersection.id();

    if (grid_.maxLevel() == 0) {
        inside.faceCenter = grid_.faceCenterEcl(inside.elemIdx, inside.faceIdx, intersection);
        outside.faceCenter = grid_.faceCenterEcl(outside.elemIdx, outside.faceIdx, intersection);
        faceAreaNormal = grid_.faceAreaNormalEcl(faceIdx);
    }
    else {
        if ((intersection.inside().level() != intersection.outside().level())) {
            // For CpGrid with LGRs, intersection laying on the boundary of an LGR, having two neighboring cells:
            // one coarse neighboring cell and one refined neighboring cell, we get the corresponding parent
            // intersection (from level 0), and use the center of the parent intersection for the coarse
            // neighboring cell.

            // Get parent intersection and its geometry
            const auto& parentIntersection =
                grid_.getParentIntersectionFromLgrBoundaryFace(intersection);
            const auto& parentIntersectionGeometry = parentIntersection.geometry();

            // For the coarse neighboring cell, take the center of the parent intersection.
            // For the refined neighboring cell, take the 'usual' center.
            inside.faceCenter =  (intersection.inside().level() == 0)
                ? parentIntersectionGeometry.center()
                : grid_.faceCenterEcl(inside.elemIdx, inside.faceIdx, intersection);
            outside.faceCenter = (intersection.outside().level() == 0)
                ?  parentIntersectionGeometry.center()
                : grid_.faceCenterEcl(outside.elemIdx, outside.faceIdx, intersection);

            // For some computations, it seems to be benefitial to replace the actual area of the refined face, by
            // the area of its parent face.
            // faceAreaNormal = parentIntersection.centerUnitOuterNormal();
            // faceAreaNormal *= parentIntersectionGeometry.volume();

            /// Alternatively, the actual area of the refined face can be computed as follows:
            faceAreaNormal = intersection.centerUnitOuterNormal();
            faceAreaNormal *= intersection.geometry().volume();
        }
        else {
            assert(intersection.inside().level() == intersection.outside().level());

            inside.faceCenter = grid_.faceCenterEcl(inside.elemIdx, inside.faceIdx, intersection);
            outside.faceCenter = grid_.faceCenterEcl(outside.elemIdx, outside.faceIdx, intersection);

            // When the CpGrid has LGRs, we compute the face area normal differently.
            if (intersection.inside().level() > 0) {  // remove intersection.inside().level() > 0
                faceAreaNormal = intersection.centerUnitOuterNormal();
                faceAreaNormal *= intersection.geometry().volume();
            }
            else {
                faceAreaNormal = grid_.faceAreaNormalEcl(faceIdx);
            }
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyPinchNncToGridTrans_(const std::unordered_map<std::size_t,int>& cartesianToCompressed)
{
    // First scale NNCs with EDITNNC.
    const auto& nnc_input = eclState_.getPinchNNC();

    for (const auto& nncEntry : nnc_input) {
        auto c1 = nncEntry.cell1;
        auto c2 = nncEntry.cell2;
        auto lowIt = cartesianToCompressed.find(c1);
        auto highIt = cartesianToCompressed.find(c2);
        int low = (lowIt == cartesianToCompressed.end())? -1 : lowIt->second;
        int high = (highIt == cartesianToCompressed.end())? -1 : highIt->second;

        if (low > high) {
            std::swap(low, high);
        }

        if (low == -1 && high == -1) {
            // Silently discard as it is not between active cells
            continue;
        }

        if (low == -1 || high == -1) {
            // We can end up here if one of the cells is overlap/ghost, because those
            // are lacking connections to other cells in the ghost/overlap.
            // Hence discard the NNC if it is between active cell and inactive cell
            continue;
        }

        {
            auto candidate = trans_.find(details::isId(low, high));
            if (candidate != trans_.end()) {
                // the correctly calculated transmissibility is stored in
                // the NNC. Overwrite previous value with it.
               candidate->second = nncEntry.trans;
            }
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNncToGridTrans_(const std::unordered_map<std::size_t,int>& cartesianToCompressed)
{
    // First scale NNCs with EDITNNC.
    const auto& nnc_input = eclState_.getInputNNC().input();

    for (const auto& nncEntry : nnc_input) {
        auto c1 = nncEntry.cell1;
        auto c2 = nncEntry.cell2;
        auto lowIt = cartesianToCompressed.find(c1);
        auto highIt = cartesianToCompressed.find(c2);
        int low = (lowIt == cartesianToCompressed.end())? -1 : lowIt->second;
        int high = (highIt == cartesianToCompressed.end())? -1 : highIt->second;

        if (low > high) {
            std::swap(low, high);
        }

        if (low == -1 && high == -1) {
            // Silently discard as it is not between active cells
            continue;
        }

        if (low == -1 || high == -1) {
            // Discard the NNC if it is between active cell and inactive cell
            std::ostringstream sstr;
            sstr << "NNC between active and inactive cells ("
                 << low << " -> " << high << ") with globalcell is (" << c1 << "->" << c2 <<")";
            OpmLog::warning(sstr.str());
            continue;
        }

        if (auto candidate = trans_.find(details::isId(low, high)); candidate != trans_.end()) {
            // NNC is represented by the grid and might be a neighboring connection
            // In this case the transmissibilty is added to the value already
            // set or computed.
            candidate->second += nncEntry.trans;
        }
        // if (enableEnergy_) {
        //     auto candidate = thermalHalfTrans_.find(details::directionalIsId(low, high));
        //     if (candidate != trans_.end()) {
        //         // NNC is represented by the grid and might be a neighboring connection
        //         // In this case the transmissibilty is added to the value already
        //         // set or computed.
        //         candidate->second += nncEntry.transEnergy1;
        //     }
        //     auto candidate = thermalHalfTrans_.find(details::directionalIsId(high, low));
        //     if (candidate != trans_.end()) {
        //         // NNC is represented by the grid and might be a neighboring connection
        //         // In this case the transmissibilty is added to the value already
        //         // set or computed.
        //         candidate->second += nncEntry.transEnergy2;
        //     }
        // }
        // if (enableDiffusivity_) {
        //     auto candidate = diffusivity_.find(details::isId(low, high));
        //     if (candidate != trans_.end()) {
        //         // NNC is represented by the grid and might be a neighboring connection
        //         // In this case the transmissibilty is added to the value already
        //         // set or computed.
        //         candidate->second += nncEntry.transDiffusion;
        //     }
        // }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal)
{
    const auto& input = eclState_.getInputNNC();
    applyEditNncToGridTransHelper_(globalToLocal, "EDITNNC",
                                   input.edit(),
                                   [&input](const NNCdata& nnc){
                                       return input.edit_location(nnc);},
                                   // Multiply transmissibility with EDITNNC value
                                   [](Scalar& trans, const Scalar& rhs){ trans *= rhs;});
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncrToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal)
{
    const auto& input = eclState_.getInputNNC();
    applyEditNncToGridTransHelper_(globalToLocal, "EDITNNCR",
                                   input.editr(),
                                   [&input](const NNCdata& nnc){
                                       return input.editr_location(nnc);},
                                   // Replace Transmissibility with EDITNNCR value
                                   [](Scalar& trans, const Scalar& rhs){ trans = rhs;});
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncToGridTransHelper_(const std::unordered_map<std::size_t,int>& globalToLocal,
                               const std::string& keyword,
                               const std::vector<NNCdata>& nncs,
                               const std::function<KeywordLocation(const NNCdata&)>& getLocation,
                               const std::function<void(Scalar&, const Scalar&)>& apply)
{
    if (nncs.empty()) {
        return;
    }
    const auto& cartDims = cartMapper_.cartesianDimensions();

    auto format_ijk = [&cartDims](std::size_t cell) -> std::string
    {
        auto i = cell % cartDims[0]; cell /= cartDims[0];
        auto j = cell % cartDims[1];
        auto k = cell / cartDims[1];

        return fmt::format("({},{},{})", i + 1,j + 1,k + 1);
    };

    auto print_warning = [&format_ijk, &getLocation, &keyword] (const NNCdata& nnc)
    {
        const auto& location = getLocation( nnc );
        auto warning =  fmt::format("Problem with {} keyword\n"
                                    "In {} line {} \n"
                                    "No NNC defined for connection {} -> {}", keyword, location.filename,
                                    location.lineno, format_ijk(nnc.cell1), format_ijk(nnc.cell2));
        OpmLog::warning(keyword, warning);
    };

    // editNnc is supposed to only reference non-neighboring connections and not
    // neighboring connections. Use all entries for scaling if there is an NNC.
    // variable nnc incremented in loop body.
    auto nnc = nncs.begin();
    auto end = nncs.end();
    std::size_t warning_count = 0;
    while (nnc != end) {
        auto c1 = nnc->cell1;
        auto c2 = nnc->cell2;
        auto lowIt = globalToLocal.find(c1);
        auto highIt = globalToLocal.find(c2);

        if (lowIt == globalToLocal.end() || highIt == globalToLocal.end()) {
            // Prevent warnings for NNCs stored on other processes in parallel (both cells inactive)
            if (lowIt != highIt && warnEditNNC_) {
                print_warning(*nnc);
                warning_count++;
            }
            ++nnc;
            continue;
        }

        auto low = lowIt->second, high = highIt->second;

        if (low > high) {
            std::swap(low, high);
        }

        auto candidate = trans_.find(details::isId(low, high));
        if (candidate == trans_.end() && warnEditNNC_) {
            print_warning(*nnc);
            ++nnc;
            warning_count++;
        }
        else {
            // NNC exists
            while (nnc != end && c1 == nnc->cell1 && c2 == nnc->cell2) {
                apply(candidate->second, nnc->trans);
                ++nnc;
            }
        }
    }

    if (warning_count > 0) {
        auto warning = fmt::format("Problems with {} keyword\n"
                                   "A total of {} connections not defined in grid", keyword, warning_count);
        OpmLog::warning(warning);
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNncMultreg_(const std::unordered_map<std::size_t,int>& cartesianToCompressed)
{
    const auto& inputNNC = this->eclState_.getInputNNC();
    const auto& transMult = this->eclState_.getTransMult();

    auto compressedIdx = [&cartesianToCompressed](const std::size_t globIdx)
    {
        auto ixPos = cartesianToCompressed.find(globIdx);
        return (ixPos == cartesianToCompressed.end()) ? -1 : ixPos->second;
    };

    // Apply region-based transmissibility multipliers (i.e., the MULTREGT
    // keyword) to those transmissibilities that are directly assigned from
    // the input.
    //
    //  * NNC::input() covers the NNC keyword and any numerical aquifers
    //  * NNC::editr() covers the EDITNNCR keyword
    //
    // Note: We do not apply MULTREGT to the entries in NNC::edit() since
    // those act as regular multipliers and have already been fully
    // accounted for in the multiplier part of the main loop of update() and
    // the applyEditNncToGridTrans_() member function.
    for (const auto& nncList : {&NNC::input, &NNC::editr}) {
        for (const auto& nncEntry : (inputNNC.*nncList)()) {
            const auto c1 = nncEntry.cell1;
            const auto c2 = nncEntry.cell2;

            auto low = compressedIdx(c1);
            auto high = compressedIdx(c2);

            if ((low == -1) || (high == -1)) {
                continue;
            }

            if (low > high) {
                std::swap(low, high);
            }

            auto candidate = this->trans_.find(details::isId(low, high));
            if (candidate != this->trans_.end()) {
                candidate->second *= transMult.getRegionMultiplierNNC(c1, c2);
            }
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeHalfTrans_(const DimVector& areaNormal,
                  int faceIdx, // in the reference element that contains the intersection
                  const DimVector& distance,
                  const DimMatrix& perm)
{
    assert(faceIdx >= 0);
    unsigned dimIdx = faceIdx / 2;
    assert(dimIdx < dimWorld);
    Scalar halfTrans = perm[dimIdx][dimIdx];
    halfTrans *= std::abs(Dune::dot(areaNormal, distance));
    halfTrans /= distance.two_norm2();

    return halfTrans;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeHalfDiffusivity_(const DimVector& areaNormal,
                        const DimVector& distance,
                        const Scalar poro)
{
    Scalar halfDiff = poro;
    halfDiff *= std::abs(Dune::dot(areaNormal, distance));
    halfDiff /= distance.two_norm2();

    return halfDiff;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
typename Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::DimVector
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
distanceVector_(const DimVector& faceCenter,
                const unsigned& cellIdx) const
{
    const auto& cellCenter = centroids_cache_.empty() ? centroids_(cellIdx)
                                                      : centroids_cache_[cellIdx];
    DimVector x = faceCenter;
    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
        x[dimIdx] -= cellCenter[dimIdx];
    }

    return x;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyMultipliers_(Scalar& trans,
                  unsigned faceIdx,
                  unsigned cartElemIdx,
                  const TransMult& transMult) const
{
    // apply multiplier for the transmissibility of the face. (the
    // face index is the index of the reference-element face which
    // contains the intersection of interest.)
    trans *= transMult.getMultiplier(cartElemIdx,
                                     FaceDir::FromIntersectionIndex(faceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNtg_(Scalar& trans,
          const FaceInfo& face,
          const std::vector<double>& ntg)
{
    // apply multiplier for the transmissibility of the face. (the
    // face index is the index of the reference-element face which
    // contains the intersection of interest.)
    // NTG does not apply to top and bottom faces
    if (face.faceIdx >= 0 && face.faceIdx <= 3) {
        trans *= ntg[face.elemIdx];
    }
}

} // namespace Opm

#endif // OPM_TRANSMISSIBILITY_IMPL_HPP
