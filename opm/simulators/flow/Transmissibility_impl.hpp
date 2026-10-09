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
#include <opm/input/eclipse/EclipseState/Grid/TranCalculator.hpp>
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
#include <string>
#include <string_view>
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
bool Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
gridJoins_(unsigned elemIdx1, unsigned elemIdx2) const
{
    using Ix = typename ElementMapper::Index;

    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    for (const auto& elem : elements(gridView_)) {
        if (elemMapper.index(elem) != static_cast<Ix>(elemIdx1)) {
            continue;
        }

        for (const auto& intersection : intersections(gridView_, elem)) {
            if (intersection.neighbor() &&
                (elemMapper.index(intersection.outside()) == static_cast<Ix>(elemIdx2)))
            {
                return true;
            }
        }

        return false;
    }

    return false;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::string Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
describeCell_(unsigned elemIdx) const
{
    auto text = fmt::format("cell {}", elemIdx);

    if (elemIdx >= static_cast<unsigned>(gridView_.size(/*codim=*/0))) {
        return text + " (out of range)";
    }

    const auto cartIdx = cartMapper_.cartesianIndex(elemIdx);
    std::array<int,dimWorld> ijk{};
    cartMapper_.cartesianCoordinate(elemIdx, ijk);

    text += fmt::format(" (Cartesian {} = [{},{},{}]", cartIdx, ijk[0], ijk[1], ijk[2]);

    if constexpr (requires { grid_.maxLevel(); }) {
        if (grid_.maxLevel() > 0) {
            using Ix = typename ElementMapper::Index;

            // Which grid the cell belongs to is the first thing worth knowing:
            // a pair spanning a refinement boundary, or two cells refining one
            // coarse cell, is where these lookups go wrong. Found through the
            // mapper -- iteration order is not index order.
            ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());
            for (const auto& elem : elements(gridView_)) {
                if (elemMapper.index(elem) == static_cast<Ix>(elemIdx)) {
                    text += fmt::format(", level {}", elem.level());
                    break;
                }
            }
        }
    }

    return text + ")";
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
lookupTrans_(const std::unordered_map<std::uint64_t, Scalar>& map,
             unsigned elemIdx1, unsigned elemIdx2, std::string_view what) const
{
    const auto entry = map.find(details::isId(elemIdx1, elemIdx2));
    if (entry == map.end()) {
        // "unordered_map::at: key not found" on its own says nothing about
        // which connection the model asked for.
        // Whether the pair is Cartesian-adjacent says which kind of connection
        // went missing: a face the grid should have, or a non-neighbour one.
        const auto& cartDims = cartMapper_.cartesianDimensions();
        const auto gc1 = cartMapper_.cartesianIndex(std::min(elemIdx1, elemIdx2));
        const auto gc2 = cartMapper_.cartesianIndex(std::max(elemIdx1, elemIdx2));
        const auto delta = gc2 - gc1;
        const auto kind =
            (delta == 0)                          ? "same coarse cell (refined siblings)" :
            (delta == 1)                          ? "neighbours in I" :
            (delta == cartDims[0])                ? "neighbours in J" :
            (delta == cartDims[0]*cartDims[1])    ? "neighbours in K" :
            ((delta % (cartDims[0]*cartDims[1])) == 0) ? "same column, several layers apart (a pinch-out or vertical NNC)"
                                                  : "not Cartesian neighbours (a non-neighbour connection)";

        // Does the grid actually join these two cells? That separates a caller
        // asking about a connection that does not exist from this calculation
        // having missed one that does.
        const auto joined = this->gridJoins_(elemIdx1, elemIdx2);

        OPM_THROW(std::out_of_range,
                  fmt::format("No {} between {} and {}: {}. {}",
                              what, this->describeCell_(elemIdx1), this->describeCell_(elemIdx2), kind,
                              joined
                              ? "The grid does join them, so the transmissibility "
                                "calculation skipped a face it should have computed."
                              : "The grid does not join them either, so the caller asked "
                                "about a connection that does not exist -- the cell "
                                "indices it used are not the ones this grid knows."));
    }

    return entry->second;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
{
    return this->lookupTrans_(trans_, elemIdx1, elemIdx2, "transmissibility");
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::optional<Scalar> Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
findTransmissibility(unsigned elemIdx1, unsigned elemIdx2) const
{
    const auto entry = trans_.find(details::isId(elemIdx1, elemIdx2));
    return (entry == trans_.end()) ? std::nullopt : std::optional<Scalar>{entry->second};
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
halfTransmissibility(unsigned insideElemIdx, unsigned outsideElemIdx) const
{
    if (!storeHalfTrans_) {
        // Without this the caller gets a bare std::out_of_range from the empty
        // map.  One branch on a bool, and this is not a hot path: only the
        // adjoint code asks for these today.
        OPM_THROW(std::logic_error,
                  "One-sided half transmissibilities were not stored. "
                  "Call setStoreHalfTrans(true) before update() to have them "
                  "computed; they are off by default because only the adjoint "
                  "code needs them.");
    }

    return halfTrans_.at(details::directionalIsId(insideElemIdx, outsideElemIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const
{
    return thermalHalfTransBoundary_.at(std::make_pair(insideElemIdx, boundaryFaceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
const std::map<std::pair<unsigned, unsigned>, Scalar>& Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
getThermalHalfTransBoundary() const
{
    return thermalHalfTransBoundary_;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
diffusivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (diffusivity_.empty())
        return 0.0;

    return this->lookupTrans_(diffusivity_, elemIdx1, elemIdx2, "diffusivity");
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
dispersivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (dispersivity_.empty())
        return 0.0;

    return this->lookupTrans_(dispersivity_, elemIdx1, elemIdx2, "dispersivity");
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

    if (storeHalfTrans_) {
        halfTrans_.clear();
        if (num_threads == 1) {
            halfTrans_.reserve(numElements*6*1.05);
        }
    }

    // if energy is enabled, let's do the same for the "thermal half transmissibilities"
    if ( enableEnergy_ && !onlyTrans) {
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
    ThreadSafeMapBuilder halfTransMap(halfTrans_, num_threads,
                                      MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder diffusivity(diffusivity_, num_threads,
                                     MapBuilderInsertionMode::Insert_Or_Assign);
    ThreadSafeMapBuilder dispersivity(dispersivity_, num_threads,
                                      MapBuilderInsertionMode::Insert_Or_Assign);

    const auto& nnc_input = eclState_.getInputNNC().input();

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

                if (storeHalfTrans_) {
                    // one-sided half transmissibilities (NTG applied),
                    // for the adjoint permeability chain rule
                    auto onesided = computeHalf(computeHalfTrans_,
                                                permeability_[inside.elemIdx],
                                                permeability_[outside.elemIdx]);
                    applyNtg_(onesided[0], inside, ntg);
                    applyNtg_(onesided[1], outside, ntg);
                    halfTransMap.insert_or_assign(
                        details::directionalIsId(inside.elemIdx, outside.elemIdx),
                        onesided[0]);
                    halfTransMap.insert_or_assign(
                        details::directionalIsId(outside.elemIdx, inside.elemIdx),
                        onesided[1]);
                }

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
                else if (!this->conformingLgr_() || inside.cartElemIdx != outside.cartElemIdx) {
                    // Same reasoning as in applyAllZMultipliers_: equal Cartesian
                    // indices mean both cells refine one coarse cell, and its
                    // MULT[XYZ] belong to its own faces, not to the faces the
                    // refinement introduces inside it.
                    applyMultipliers_(trans, inside.faceIdx, inside.cartElemIdx, transMult);
                    // ... and outside elements
                    applyMultipliers_(trans, outside.faceIdx, outside.cartElemIdx, transMult);
                }

                bool foundInputNNC = false;
                if (! nnc_input.empty()) {
                    // Skip region multipliers for overlapping input NNCs (they are handled later)
                    auto it = std::lower_bound(nnc_input.begin(), nnc_input.end(),
                                               NNCdata { inside.cartElemIdx, outside.cartElemIdx, 0.0 });
                    foundInputNNC = it != nnc_input.end() && it->cell1 == inside.cartElemIdx && it->cell2 == outside.cartElemIdx;
                }
                if (! foundInputNNC) {
                    // apply the region multipliers (cf. the MULTREGT keyword)
                    trans *= transMult.getRegionMultiplier(inside.cartElemIdx,
                                                           outside.cartElemIdx,
                                                           faceIdToDir(inside.faceIdx));
                }

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
#ifdef _OPENMP
#pragma omp section
#endif
        halfTransMap.finalize();
    }

    if (this->conformingLgr_()) {
        this->computeHostLevelTrans_();
    }

    // Potentially overwrite and/or modify transmissibilities based on input from deck
    this->updateFromEclState_(global);

    // Create mapping from global to local index
    std::unordered_map<std::size_t,int> globalToLocal;

    // A refined leaf cell's Cartesian index is its father's, so all children of
    // one deck cell collide on a single key.  Keep the last-one-wins map for the
    // paths that need one cell, and a father-to-children map for EDITNNC, whose
    // multiplier belongs to every face the coarse connection became.
    CartesianToLeaf globalToChildren;

    // Loop over all elements (global grid) and store Cartesian index
    for (const auto& elem : elements(grid_.leafGridView())) {
        int elemIdx = elemMapper.index(elem);
        int cartElemIdx =  cartMapper_.cartesianIndex(elemIdx);
        globalToLocal[cartElemIdx] = elemIdx;
        globalToChildren[cartElemIdx].push_back(elemIdx);
    }

    if (!disableNNC) {
        // For EDITNNC and EDITNNCR we warn only once
        // If transmissibility is used for load balancing this will be done
        // when computing the gobal transmissibilities and all warnings will
        // be seen in a parallel. Unfortunately, when we do not use transmissibilities
        // we will only see warnings for the partition of process 0 and also false positives.
        this->applyPinchNncToGridTrans_(globalToChildren, applyNncMultregT);
        this->applyNncToGridTrans_(globalToChildren);
        this->applyEditNncToGridTrans_(globalToChildren);
        this->applyEditNncrToGridTrans_(globalToChildren);
        if (applyNncMultregT) {
            this->applyNncMultreg_(globalToChildren);
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
            std::ranges::copy(por, porosity_.begin());
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
        std::ranges::copy(disp, dispersion_.begin());
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
    // Both cells refine the same coarse cell, so this face is interior to it.
    // MULT[XYZ] describe that coarse cell's own faces, which the refinement
    // inherits on its outer boundary; an interior face has no coarse face to
    // take a multiplier from, and the pillar walk below -- which steps through
    // coarse cells -- has nothing to walk. Taking the coarse cell's own
    // Z+ times Z- here, as this did before refusing LGRs outright, would damp
    // every interior face by a multiplier meant for the coarse cell's top and
    // bottom.
    if (grid_.maxLevel() > 0 && !this->conformingLgr_()) {
        OPM_THROW(std::invalid_argument, "MULTZ not support with LGRS, yet.");
    }
    if (inside.cartElemIdx == outside.cartElemIdx) {
        return;
    }

    if (inside.faceIdx > 3) { // top or or bottom
        assert(inside.faceIdx == 5); // as insideCartElemIdx < outsideCartElemIdx holds for the Z column
        // For CpGrid with LGRs, insideCartElemIdx == outsideCartElemIdx when cells on the leaf have the same parent cell on level zero.
        assert(outside.cartElemIdx >= inside.cartElemIdx);
        const unsigned lastCartElemIdx = outside.cartElemIdx - cartDims[0]*cartDims[1];
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

    if (grid_.maxLevel() == 0) {
        std::array<std::vector<double>,3> trans = createTransmissibilityArrays_(is_tran);
        for (int dir = 0; dir < 3; ++dir) {
            if (is_tran[dir]) {
                fp->apply_tran(keywords[dir], trans[dir]);
            }
        }
        resetTransmissibilityFromArrays_(is_tran, trans);
        return;
    }

    // A refined grid has several faces per cell and direction, so work on a
    // list of faces, each with the coarse cell whose TRAN* modifier applies.
    struct RefinedFace { std::uint64_t id; int action; };
    std::array<std::vector<RefinedFace>,3> faces;
    std::array<std::vector<double>,3> values;
    std::array<std::vector<int>,3> index;
    // The coarse face each modifier belongs to, as a pair of Cartesian cells.
    std::array<std::map<int, std::pair<int,int>>,3> coarseFace;
    {
        const auto& cartDims = cartMapper_.cartesianDimensions();
        ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());
        for (const auto& elem : elements(gridView_)) {
            for (const auto& is : intersections(gridView_, elem)) {
                if (!is.neighbor()) {
                    continue;
                }
                const unsigned c1 = elemMapper.index(is.inside());
                const unsigned c2 = elemMapper.index(is.outside());
                const int gc1 = cartMapper_.cartesianIndex(c1);
                const int gc2 = cartMapper_.cartesianIndex(c2);
                if (std::tie(gc1, c1) > std::tie(gc2, c2)) {
                    continue;
                }
                const auto f = is.indexInInside();
                int dir = -1;
                if ((gc2 - gc1 == 1 || (gc2 == gc1 && (f == 0 || f == 1))) && cartDims[0] > 1) {
                    dir = 0;
                }
                else if ((gc2 - gc1 == cartDims[0] || (gc2 == gc1 && (f == 2 || f == 3))) && cartDims[1] > 1) {
                    dir = 1;
                }
                else if (gc2 - gc1 == cartDims[0]*cartDims[1] || (gc2 == gc1 && (f == 4 || f == 5))) {
                    dir = 2;
                }
                if ((dir < 0) || !is_tran[dir]) {
                    continue;
                }
                // A face interior to one coarse cell has no coarse face of its own.
                const int action = (gc1 == gc2) ? -1
                    : static_cast<int>(lookUpData_.template getFieldPropIdx<Grid>(static_cast<int>(c1)));
                const auto id = details::isId(c1, c2);
                if (action >= 0) {
                    coarseFace[dir].emplace(action, std::make_pair(gc1, gc2));
                }
                faces[dir].push_back({id, action});
                values[dir].push_back(trans_[id]);
                index[dir].push_back(action);
            }
        }
    }

    for (int dir = 0; dir < 3; ++dir) {
        if (!is_tran[dir]) {
            continue;
        }
        const auto& key = keywords[dir];
        auto& value = values[dir];
        const auto numCoarse = static_cast<std::size_t>(grid_.levelGridView(0).size(0));

        // Each coarse face's own transmissibility: the host-level value where a
        // box refines it, else the sum of its (one) leaf face.
        auto coarse = std::vector<double>(numCoarse, 0.0);
        for (std::size_t e = 0; e < value.size(); ++e) {
            if (index[dir][e] >= 0) {
                coarse[index[dir][e]] += value[e];
            }
        }
        for (const auto& [c, cells] : coarseFace[dir]) {
            if (const auto host = hostLevelTrans_.find(cells); host != hostLevelTrans_.end()) {
                coarse[c] = host->second;
            }
        }
        auto edited = coarse;
        fp->apply_tran(key, edited);
        for (const auto& [c, cells] : coarseFace[dir]) {
            if (auto host = hostLevelTrans_.find(cells); host != hostLevelTrans_.end()) {
                host->second = edited[c];
            }
        }

        const auto ops = fp->tran_operations(key);
        if (!ops.count(Fieldprops::ScalarOperation::EQUAL) &&
            !ops.count(Fieldprops::ScalarOperation::ADD))
        {
            // MUL, MIN and MAX carry over to each refined face of a coarse face.
            fp->apply_tran(key, index[dir], value);
        }
        else {
            // EQUALS and ADD state the coarse face's transmissibility; its refined
            // faces scale by the coarse face's change.
            auto count = std::vector<int>(numCoarse, 0);
            for (const auto c : index[dir]) {
                if (c >= 0) {
                    ++count[c];
                }
            }
            for (std::size_t e = 0; e < value.size(); ++e) {
                if (const auto c = index[dir][e]; c >= 0) {
                    value[e] = (coarse[c] > 0.0) ? value[e] * (edited[c] / coarse[c])
                                                 : edited[c] / count[c];
                }
            }
        }
        for (std::size_t e = 0; e < value.size(); ++e) {
            trans_[faces[dir][e].id] = value[e];
        }
    }
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
            // The coarse cell measures to its whole face's corner average, as the reference.
            auto coarseFaceCenter = [](const auto& element, const int face) {
                DimVector center(0.0);
                const auto& geometry = element.geometry();
                for (int corner = 0; corner < 8; ++corner) {
                    if (((corner >> (face / 2)) & 1) == face % 2) {
                        for (int d = 0; d < 3; ++d) {
                            center[d] += geometry.corner(corner)[d] / 4.0;
                        }
                    }
                }
                return center;
            };
            // The Trilinear backend measures the coarse cell to its parent intersection instead.
            const bool conforming = this->conformingLgr_();
            const auto parentCenter = conforming ? DimVector(0.0)
                : DimVector(grid_.getParentIntersectionFromLgrBoundaryFace(intersection).geometry().center());
            inside.faceCenter =  (intersection.inside().level() == 0)
                ? (conforming ? coarseFaceCenter(intersection.inside(), inside.faceIdx) : parentCenter)
                : grid_.faceCenterEcl(inside.elemIdx, inside.faceIdx, intersection);
            outside.faceCenter = (intersection.outside().level() == 0)
                ? (conforming ? coarseFaceCenter(intersection.outside(), outside.faceIdx) : parentCenter)
                : grid_.faceCenterEcl(outside.elemIdx, outside.faceIdx, intersection);

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
applyPinchNncToGridTrans_(const CartesianToLeaf& cartesianToCompressed,
                          const bool applyNncMultregT)
{
    const auto& pinchNnc = eclState_.getPinchNNC();
    const auto& transMult = this->eclState_.getTransMult();

    constexpr std::size_t maxReported = 5;
    std::vector<std::pair<std::size_t,std::size_t>> dropped{};
    std::size_t numDropped = 0;

    for (const auto& nncEntry : pinchNnc) {
        auto c1 = nncEntry.cell1;
        auto c2 = nncEntry.cell2;
        auto lowIt = cartesianToCompressed.find(c1);
        auto highIt = cartesianToCompressed.find(c2);

        if ((lowIt == cartesianToCompressed.end()) ||
            (highIt == cartesianToCompressed.end()))
        {
            // Not between two active cells, or one end is an overlap/ghost cell
            // that carries no connections. Discard, as before.
            continue;
        }

        // Under refinement the coarse connection became a connection between each
        // pair of the two cells' children that the grid actually joined. The
        // transmissibility is an absolute one, so it has to be divided among them:
        // n pairs each take r_d/n of it, r_d being the refinement in the connection's
        // own direction, which is the children-per-cell over the pair count. Without
        // refinement there is one child each and this is the single face as before.
        auto faces = std::vector<decltype(trans_.begin())>{};
        for (const auto childLow : lowIt->second) {
            for (const auto childHigh : highIt->second) {
                auto low = childLow, high = childHigh;
                if (low > high) {
                    std::swap(low, high);
                }
                auto candidate = trans_.find(details::isId(low, high));
                if (candidate != trans_.end()) {
                    faces.push_back(candidate);
                }
            }
        }

        if (faces.empty()) {
            if (dropped.size() < maxReported) {
                dropped.push_back(std::make_pair(c1, c2));
            }
            ++numDropped;
            continue;
        }

        const auto n = static_cast<Scalar>(faces.size());
        const auto children = static_cast<Scalar>(std::max(lowIt->second.size(),
                                                           highIt->second.size()));
        auto value = nncEntry.trans * (children / (n * n));
        if (applyNncMultregT) {
            value *= transMult.getRegionMultiplierNNC(c1, c2);
        }

        for (auto& face : faces) {
            face->second = value;
        }
    }

    if (numDropped > 0) {
        OpmLog::warning(fmt::format(
            "{} PINCH connection(s) were computed but the grid holds no face to "
            "carry them -- with refinement, no face between any pair of the two "
            "cells' children -- so the pinched-out layers they bridge are not "
            "bridged.{}",
            numDropped, this->describeDroppedNnc_(dropped, numDropped)));
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNncToGridTrans_(const CartesianToLeaf& cartesianToCompressed)
{
    // First scale NNCs with EDITNNC.
    const auto& nnc_input = eclState_.getInputNNC().input();

    constexpr std::size_t maxReported = 5;
    std::vector<std::pair<std::size_t,std::size_t>> unconnectedNnc{};
    std::size_t numUnconnectedNnc = 0;

    // A deck cell that several leaf cells descend from is inside a refinement box.
    auto refined = [&cartesianToCompressed](const std::size_t cart)
    {
        const auto it = cartesianToCompressed.find(cart);
        return (it != cartesianToCompressed.end()) && (it->second.size() > 1);
    };

    for (const auto& nncEntry : nnc_input) {
        auto c1 = nncEntry.cell1;
        auto c2 = nncEntry.cell2;
        auto lowIt = cartesianToCompressed.find(c1);
        auto highIt = cartesianToCompressed.find(c2);
        int low = (lowIt == cartesianToCompressed.end())? -1 : lowIt->second.front();
        int high = (highIt == cartesianToCompressed.end())? -1 : highIt->second.front();

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

        if (refined(c1) || refined(c2)) {
            OPM_THROW(std::invalid_argument,
                      fmt::format("An explicit connection -- the NNC keyword or a "
                                  "numerical aquifer -- names cell {} or {}, which "
                                  "a refinement box covers. Its transmissibility is "
                                  "an absolute one and there is no rule yet for "
                                  "dividing it among the faces the connection became, "
                                  "so it would be applied to one arbitrary pair of "
                                  "children or to none. Move the box off the "
                                  "connection.",
                                  ijkString_(c1), ijkString_(c2)));
        }

        if (auto candidate = trans_.find(details::isId(low, high)); candidate != trans_.end()) {
            // NNC is represented by the grid and might be a neighboring connection
            // In this case the transmissibilty is added to the value already
            // set or computed.
            candidate->second += nncEntry.trans;
        }
        else {
            // Both cells are active, yet the grid holds no connection between
            // them, so this entry's transmissibility goes nowhere.  Under
            // refinement that is the normal outcome for a connection naming a
            // cell inside a CARFIN box: the coarse cell is gone from the leaf
            // and its NNC face with it.  Silence here reads as a working
            // connection, which for a numerical aquifer means it quietly stops
            // feeding the reservoir.
            if (unconnectedNnc.size() < maxReported) {
                unconnectedNnc.push_back(std::make_pair(c1, c2));
            }
            ++numUnconnectedNnc;
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

    if (numUnconnectedNnc > 0) {
        const auto cells = this->describeDroppedNnc_(unconnectedNnc, numUnconnectedNnc);

        OpmLog::warning(fmt::format
                        ("{} explicit connection(s) -- NNC, EDITNNC or a numerical aquifer -- "
                         "name a cell pair the grid does not join, so their transmissibility is "
                         "not applied and no flow passes through them. With local grid refinement "
                         "this is what happens when a connection names a cell inside a CARFIN box: "
                         "the coarse cell is not on the leaf grid and neither is its connection. "
                         "A numerical aquifer in that position stops feeding the reservoir "
                         "entirely.{}", numUnconnectedNnc, cells));
    }
}


template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
bool Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
conformingLgr_() const
{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {
        return grid_.lgrBackend() == Opm::Refinement::Backend::Conforming;
    }
    else {
        return false;
    }
}

// The level-zero hosts' own connections as the unrefined grid has them: plain
// geometry, as in the reference's global section, and what TRAN* edits act on.
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeHostLevelTrans_()
{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {
        if (grid_.maxLevel() == 0) {
            return;
        }

        using LevelView = std::remove_const_t<decltype(grid_.levelGridView(0))>;
        const LevelView level0 = grid_.levelGridView(0);
        const auto levelMapper = Dune::MultipleCodimMultipleGeomTypeMapper<LevelView>
            { level0, Dune::mcmgElementLayout() };

        const auto& fp = eclState_.fieldProps();
        const auto permx = fp.get_double("PERMX");
        const auto permy = fp.has_double("PERMY") ? fp.get_double("PERMY") : permx;
        const auto permz = fp.has_double("PERMZ") ? fp.get_double("PERMZ") : permx;
        const auto ntgArr = fp.has_double("NTG")
            ? fp.get_double("NTG") : std::vector<double>(permx.size(), 1.0);

        const auto numHost = static_cast<std::size_t>(level0.size(/*codim=*/0));
        if (permx.size() < numHost) {
            OPM_THROW(std::invalid_argument,
                      "The level-zero grid has more cells than the field properties "
                      "describe; the host cells' transmissibilities cannot be computed.");
        }

        // Per neighbouring host pair: the half-transmissibility and face area. A
        // faulted host face is several intersections, so accumulate the pieces.
        struct Face { Scalar half{0}, area{0}; };
        auto isHost = std::vector<bool>(numHost, false);
        for (const auto& elem : elements(gridView_)) {
            if (elem.level() > 0) {
                isHost[levelMapper.index(elem.getOrigin())] = true;
            }
        }
        auto hostPair = std::map<std::tuple<std::size_t,int,std::size_t>, Face>{};


        auto level0Cart = std::vector<int>(numHost);
        for (const auto& elem : elements(level0)) {
            const auto idx = levelMapper.index(elem);
            level0Cart[idx] = elem.getLevelCartesianIdx();

            // The unrefined grid's convention: centres are corner averages.
            const auto& cellGeom = elem.geometry();
            DimVector centre(0.0);
            for (int c = 0; c < 8; ++c) {
                centre += cellGeom.corner(c);
            }
            centre /= 8.0;
            auto faceCentre = [&cellGeom](const int f) {
                DimVector fc(0.0);
                for (int c = 0; c < 8; ++c) {
                    if (((c >> (f / 2)) & 1) == (f % 2)) {
                        fc += cellGeom.corner(c);
                    }
                }
                fc /= 4.0;
                return fc;
            };

            DimMatrix K(0.0);
            K[0][0] = permx[idx]; K[1][1] = permy[idx]; K[2][2] = permz[idx];

            for (const auto& is : intersections(level0, elem)) {
                const auto f = is.indexInInside();
                if ((f < 0) || (f > 5)) {
                    continue;               // NNC: no face of the reference element
                }

                const auto& geom = is.geometry();
                // The unrefined grid's area convention, as the main loop uses on level zero.
                DimVector areaNormal = grid_.faceAreaNormalEcl(is.id(), 0);

                DimVector d = faceCentre(f);
                d -= centre;

                auto half = computeHalfTrans_(areaNormal, f, d, K);
                if (f < 4) {                // lateral faces carry NTG
                    half *= ntgArr[idx];
                }

                if (is.neighbor() && (isHost[idx] || isHost[levelMapper.index(is.outside())])) {
                    auto& pair = hostPair[{idx, f, levelMapper.index(is.outside())}];
                    pair.half += half;
                    pair.area += geom.volume();
                }
            }
        }

        // The deck's face and region multipliers, applied as the main loop does.
        const auto& transMult = eclState_.getTransMult();
        const auto& nncInput = eclState_.getInputNNC().input();
        auto multiplier = [&](const int c1, const int f1, const int c2, const int f2)
        {
            Scalar mult = 1.0;
            if (c1 != c2) {
                mult *= transMult.getMultiplier(c1, FaceDir::FromIntersectionIndex(f1))
                      * transMult.getMultiplier(c2, FaceDir::FromIntersectionIndex(f2));
            }
            const auto it = std::lower_bound(nncInput.begin(), nncInput.end(),
                                             NNCdata { static_cast<std::size_t>(c1),
                                                       static_cast<std::size_t>(c2), 0.0 });
            const bool isInputNnc = (it != nncInput.end()) &&
                (it->cell1 == static_cast<std::size_t>(c1)) &&
                (it->cell2 == static_cast<std::size_t>(c2));
            if (!isInputNnc) {
                mult *= transMult.getRegionMultiplier(c1, c2, (f1 < 2) ? FaceDir::XPlus
                                                            : (f1 < 4) ? FaceDir::YPlus
                                                                       : FaceDir::ZPlus);
            }
            return mult;
        };

        // The hosts' own connections as the unrefined grid has them, for the
        // output's global section; as in the reference, without EDITNNC.
        hostLevelTrans_.clear();
        {
            auto pinchAll = std::map<std::pair<std::size_t,std::size_t>, Scalar>{};
            for (const auto& nnc : eclState_.getPinchNNC()) {
                pinchAll[{std::min(nnc.cell1, nnc.cell2), std::max(nnc.cell1, nnc.cell2)}] = nnc.trans;
            }
            auto halves = std::map<std::pair<std::size_t,std::size_t>, std::pair<Scalar,int>>{};
            for (const auto& [key, face] : hostPair) {
                auto& half = halves[{std::get<0>(key), std::get<2>(key)}];
                half.first += face.half;
                half.second = std::get<1>(key);
            }
            for (const auto& [key, half] : halves) {
                const auto [a, b] = key;
                const auto reverse = halves.find({b, a});
                if ((a > b) || (reverse == halves.end()) ||
                    !(half.first > 0.0) || !(reverse->second.first > 0.0)) {
                    continue;
                }
                auto ca = level0Cart[a], cb = level0Cart[b];
                auto fa = half.second, fb = reverse->second.second;
                if (ca > cb) {
                    std::swap(ca, cb);
                    std::swap(fa, fb);
                }
                Scalar t = multiplier(ca, fa, cb, fb) / (1.0 / half.first + 1.0 / reverse->second.first);
                const auto pair = std::make_pair(static_cast<std::size_t>(ca), static_cast<std::size_t>(cb));
                if (const auto pinch = pinchAll.find(pair); pinch != pinchAll.end()) {
                    t = pinch->second * transMult.getRegionMultiplierNNC(pair.first, pair.second);
                }
                hostLevelTrans_[{ca, cb}] = t;
            }
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::array<int,3>
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
ijkFromCartesian_(const std::size_t cartIdx) const
{
    const auto& dims = eclState_.getInputGrid().getNXYZ();
    const auto i = cartIdx % dims[0];
    const auto j = (cartIdx / dims[0]) % dims[1];
    const auto k = cartIdx / (static_cast<std::size_t>(dims[0]) * dims[1]);

    return { static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) };
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::string
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
ijkString_(const std::size_t cartIdx) const
{
    const auto ijk = ijkFromCartesian_(cartIdx);
    return fmt::format("({},{},{})", ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::string
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
describeDroppedNnc_(const std::vector<std::pair<std::size_t,std::size_t>>& sample,
                    const std::size_t total) const
{
    auto cells = std::string{};
    for (const auto& [c1, c2] : sample) {
        const auto ijk1 = ijkFromCartesian_(c1);
        const auto ijk2 = ijkFromCartesian_(c2);
        cells += fmt::format("\n  ({},{},{}) -- ({},{},{})",
                             ijk1[0] + 1, ijk1[1] + 1, ijk1[2] + 1,
                             ijk2[0] + 1, ijk2[1] + 1, ijk2[2] + 1);
    }
    if (total > sample.size()) {
        cells += fmt::format("\n  ... and {} more", total - sample.size());
    }

    return cells;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncToGridTrans_(const CartesianToLeaf& globalToLocal)
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
applyEditNncrToGridTrans_(const CartesianToLeaf& globalToLocal)
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
applyEditNncToGridTransHelper_(const CartesianToLeaf& globalToLocal,
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

    // As in the reference, a record does not reach a connection of a refined
    // cell: it names the coarse cells, and their connection is replaced.
    std::unordered_set<std::size_t> refined;
    const ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());
    for (const auto& elem : elements(gridView_)) {
        if (elem.level() > 0) {
            refined.insert(cartMapper_.cartesianIndex(elemMapper.index(elem)));
        }
    }
    std::size_t skipped = 0;

    // editNnc is supposed to only reference non-neighboring connections and not
    // neighboring connections. Use all entries for scaling if there is an NNC.
    // variable nnc incremented in loop body.
    auto nnc = nncs.begin();
    auto end = nncs.end();
    std::size_t warning_count = 0;
    while (nnc != end) {
        auto c1 = nnc->cell1;
        auto c2 = nnc->cell2;
        if (refined.count(c1) || refined.count(c2)) {
            ++skipped;
            ++nnc;
            continue;
        }
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

        // Collect the faces first: a record may be repeated, and each repeat
        // must act on each face exactly once.
        auto faces = std::vector<decltype(trans_.begin())>{};
        for (const auto childLow : lowIt->second) {
            for (const auto childHigh : highIt->second) {
                auto low = childLow, high = childHigh;
                if (low > high) {
                    std::swap(low, high);
                }

                auto candidate = trans_.find(details::isId(low, high));
                if (candidate != trans_.end()) {
                    faces.push_back(candidate);
                }
            }
        }

        if (faces.empty()) {
            if (warnEditNNC_) {
                print_warning(*nnc);
                warning_count++;
            }
            ++nnc;
        }
        else {
            while (nnc != end && c1 == nnc->cell1 && c2 == nnc->cell2) {
                for (auto& face : faces) {
                    apply(face->second, nnc->trans);
                }
                ++nnc;
            }
        }
    }

    if (warning_count > 0) {
        auto warning = fmt::format("Problems with {} keyword\n"
                                   "A total of {} connections not defined in grid", keyword, warning_count);
        OpmLog::warning(warning);
    }
    if ((skipped > 0) && warnEditNNC_) {
        OpmLog::info(fmt::format("{} {} record(s) name a refined cell and are not applied "
                                 "to the refined connections.", skipped, keyword));
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNncMultreg_(const CartesianToLeaf& cartesianToCompressed)
{
    const auto& inputNNC = this->eclState_.getInputNNC();
    const auto& transMult = this->eclState_.getTransMult();

    auto compressedIdx = [&cartesianToCompressed](const std::size_t globIdx)
    {
        auto ixPos = cartesianToCompressed.find(globIdx);
        return (ixPos == cartesianToCompressed.end()) ? -1 : ixPos->second.front();
    };

    // A deck cell that several leaf cells descend from is inside a refinement box.
    auto refined = [&cartesianToCompressed](const std::size_t cart)
    {
        const auto it = cartesianToCompressed.find(cart);
        return (it != cartesianToCompressed.end()) && (it->second.size() > 1);
    };

    constexpr std::size_t maxReported = 5;
    std::vector<std::pair<std::size_t,std::size_t>> dropped{};
    std::size_t numDropped = 0;

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
            // As in applyEditNncrToGridTrans_, EDITNNCR does not reach a refined cell.
            if ((nncList == &NNC::editr) && (refined(c1) || refined(c2))) {
                continue;
            }

            if (low > high) {
                std::swap(low, high);
            }

            const auto mult = transMult.getRegionMultiplierNNC(c1, c2);
            if ((mult != Scalar{1}) && (refined(c1) || refined(c2))) {
                OPM_THROW(std::invalid_argument,
                          fmt::format("MULTREGT gives the connection {} -- {} a "
                                      "multiplier of {}, and a refinement box covers "
                                      "one of those cells. The connection became "
                                      "several faces there and the multiplier reaches "
                                      "at most one of them, so the region boundary "
                                      "would be left open. Move the box off the "
                                      "connection, or set the multiplier to 1.",
                                      ijkString_(c1), ijkString_(c2), mult));
            }

            auto candidate = this->trans_.find(details::isId(low, high));
            if (candidate != this->trans_.end()) {
                candidate->second *= mult;
            }
            else if (mult != Scalar{1}) {
                if (dropped.size() < maxReported) {
                    dropped.push_back(std::make_pair(c1, c2));
                }
                ++numDropped;
            }
        }
    }

    if (numDropped > 0) {
        OpmLog::warning(fmt::format(
            "{} MULTREGT multiplier(s) name a connection the grid does not hold, so "
            "the region boundary they seal is left open.{}",
            numDropped, this->describeDroppedNnc_(dropped, numDropped)));
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
