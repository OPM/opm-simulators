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
#ifndef EWOMS_ECL_TRANSMISSIBILITY_IMPL_HH
#define EWOMS_ECL_TRANSMISSIBILITY_IMPL_HH

#include <ebos/ecltransmissibility.hh>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <fmt/format.h>

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

namespace {

constexpr unsigned elemIdxShift = 32; // bits

std::uint64_t isId(std::uint32_t elemIdx1, std::uint32_t elemIdx2)
{
    std::uint32_t elemAIdx = std::min(elemIdx1, elemIdx2);
    std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

    return (elemBIdx<<elemIdxShift) + elemAIdx;
}

std::pair<std::uint32_t, std::uint32_t> isIdReverse(const std::uint64_t& id)
{
    // Assigning an unsigned integer to a narrower type discards the most significant bits.
    // See "The C programming language", section A.6.2.
    // NOTE that the ordering of element A and B may have changed
    std::uint32_t elemAIdx = id;
    std::uint32_t elemBIdx = (id - elemAIdx) >> elemIdxShift;

    return std::make_pair(elemAIdx, elemBIdx);
}

std::uint64_t directionalIsId(std::uint32_t elemIdx1, std::uint32_t elemIdx2)
{
    return (std::uint64_t(elemIdx1)<<elemIdxShift) + elemIdx2;
}

}

namespace Opm {

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
EclTransmissibility(const EclipseState& eclState,
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
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
{
    return trans_.at(isId(elemIdx1, elemIdx2));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
transmissibilityBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
{
    return transBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
thermalHalfTrans(unsigned insideElemIdx, unsigned outsideElemIdx) const
{
    return thermalHalfTrans_.at(directionalIsId(insideElemIdx, outsideElemIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const
{
    return thermalHalfTransBoundary_.at(std::make_pair(insideElemIdx, boundaryFaceIdx));
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
diffusivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (diffusivity_.empty())
        return 0.0;

    return diffusivity_.at(isId(elemIdx1, elemIdx2));

}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
Scalar EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
dispersivity(unsigned elemIdx1, unsigned elemIdx2) const
{
    if (dispersivity_.empty())
        return 0.0;

    return dispersivity_.at(isId(elemIdx1, elemIdx2));

}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
update(bool global, const std::function<unsigned int(unsigned int)>& map, const bool applyNncMultregT)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    const auto& transMult = eclState_.getTransMult();
    const auto& comm = gridView_.comm();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    unsigned numElements = elemMapper.size();
    // get the ntg values, the ntg values are modified for the cells merged with minpv
    const std::vector<double>& ntg = this->lookUpData_.assignFieldPropsDoubleOnLeaf(eclState_.fieldProps(), "NTG", numElements);
    const bool updateDiffusivity = eclState_.getSimulationConfig().isDiffusive();
    const bool updateDispersivity = eclState_.getSimulationConfig().rock_config().dispersion();

    if (map)
        extractPermeability_(map);
    else
        extractPermeability_();

    // calculate the axis specific centroids of all elements
    std::array<std::vector<DimVector>, dimWorld> axisCentroids;

    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        axisCentroids[dimIdx].resize(numElements);

    for (const auto& elem : elements(gridView_)) {
        unsigned elemIdx = elemMapper.index(elem);

        // compute the axis specific "centroids" used for the transmissibilities. for
        // consistency with the flow simulator, we use the element centers as
        // computed by opm-parser's Opm::EclipseGrid class for all axes.
        std::array<double, dimWorld> centroid = centroids_(elemIdx);

        for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                axisCentroids[axisIdx][elemIdx][dimIdx] = centroid[dimIdx];
    }

    // reserving some space in the hashmap upfront saves quite a bit of time because
    // resizes are costly for hashmaps and there would be quite a few of them if we
    // would not have a rough idea of how large the final map will be (the rough idea
    // is a conforming Cartesian grid).
    trans_.clear();
    trans_.reserve(numElements*3*1.05);

    transBoundary_.clear();

    // if energy is enabled, let's do the same for the "thermal half transmissibilities"
    if (enableEnergy_) {
        thermalHalfTrans_.clear();
        thermalHalfTrans_.reserve(numElements*6*1.05);

        thermalHalfTransBoundary_.clear();
    }

    // if diffusion is enabled, let's do the same for the "diffusivity"
    if (updateDiffusivity) {
        diffusivity_.clear();
        diffusivity_.reserve(numElements*3*1.05);
        extractPorosity_();
    }

    // if dispersion is enabled, let's do the same for the "dispersivity"
    if (updateDispersivity) {
        dispersivity_.clear();
        dispersivity_.reserve(numElements*3*1.05);
        extractDispersion_();
    }

    // The MULTZ needs special case if the option is ALL
    // Then the smallest multiplier is applied.
    // Default is to apply the top and bottom multiplier
    bool useSmallestMultiplier;
    bool pinchActive;
    if (comm.rank() == 0) {
        const auto& eclGrid = eclState_.getInputGrid();
        pinchActive = eclGrid.isPinchActive();
        useSmallestMultiplier = eclGrid.getMultzOption() == PinchMode::ALL;
    }
    if (global && comm.size() > 1) {
        comm.broadcast(&useSmallestMultiplier, 1, 0);
        comm.broadcast(&pinchActive, 1, 0);
    }

    // compute the transmissibilities for all intersections
    for (const auto& elem : elements(gridView_)) {
        unsigned elemIdx = elemMapper.index(elem);

        auto isIt = gridView_.ibegin(elem);
        const auto& isEndIt = gridView_.iend(elem);
        unsigned boundaryIsIdx = 0;
        for (; isIt != isEndIt; ++ isIt) {
            // store intersection, this might be costly
            const auto& intersection = *isIt;

            // deal with grid boundaries
            if (intersection.boundary()) {
                // compute the transmissibilty for the boundary intersection
                const auto& geometry = intersection.geometry();
                const auto& faceCenterInside = geometry.center();

                auto faceAreaNormal = intersection.centerUnitOuterNormal();
                faceAreaNormal *= geometry.volume();

                Scalar transBoundaryIs;
                computeHalfTrans_(transBoundaryIs,
                                  faceAreaNormal,
                                  intersection.indexInInside(),
                                  distanceVector_(faceCenterInside,
                                                  intersection.indexInInside(),
                                                  elemIdx,
                                                  axisCentroids),
                                  permeability_[elemIdx]);

                // normally there would be two half-transmissibilities that would be
                // averaged. on the grid boundary there only is the half
                // transmissibility of the interior element.
                unsigned insideCartElemIdx = cartMapper_.cartesianIndex(elemIdx);
                applyMultipliers_(transBoundaryIs, intersection.indexInInside(), insideCartElemIdx, transMult);
                transBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] = transBoundaryIs;

                // for boundary intersections we also need to compute the thermal
                // half transmissibilities
                if (enableEnergy_) {
                    Scalar transBoundaryEnergyIs;
                    computeHalfDiffusivity_(transBoundaryEnergyIs,
                                            faceAreaNormal,
                                            distanceVector_(faceCenterInside,
                                                            intersection.indexInInside(),
                                                            elemIdx,
                                                            axisCentroids),
                                            1.0);
                    thermalHalfTransBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] =
                        transBoundaryEnergyIs;
                }

                ++ boundaryIsIdx;
                continue;
            }

            if (!intersection.neighbor()) {
                // elements can be on process boundaries, i.e. they are not on the
                // domain boundary yet they don't have neighbors.
                ++ boundaryIsIdx;
                continue;
            }

            const auto& outsideElem = intersection.outside();
            unsigned outsideElemIdx = elemMapper.index(outsideElem);

            // Get the Cartesian indices of the origen cells (parent or equivalent cell on level zero), for CpGrid with LGRs.
            // For genral grids and no LGRs, get the usual Cartesian Index.
            unsigned insideCartElemIdx = this-> lookUpCartesianData_.template getFieldPropCartesianIdx<Grid>(elemIdx);
            unsigned outsideCartElemIdx =  this-> lookUpCartesianData_.template getFieldPropCartesianIdx<Grid>(outsideElemIdx);

            // we only need to calculate a face's transmissibility
            // once...
            // In a parallel run insideCartElemIdx>outsideCartElemIdx does not imply elemIdx>outsideElemIdx for
            // ghost cells and we need to use the cartesian index as this will be used when applying Z multipliers
            // We still need to cover the case where both cells are part of an LGR and as a consequence might have
            // the same cartesian index
            if (std::tie(insideCartElemIdx, elemIdx) > std::tie(outsideCartElemIdx, outsideElemIdx))
                continue;

            // local indices of the faces of the inside and
            // outside elements which contain the intersection
            int insideFaceIdx  = intersection.indexInInside();
            int outsideFaceIdx = intersection.indexInOutside();

            if (insideFaceIdx == -1) {
                // NNC. Set zero transmissibility, as it will be
                // *added to* by applyNncToGridTrans_() later.
                assert(outsideFaceIdx == -1);
                trans_[isId(elemIdx, outsideElemIdx)] = 0.0;
                if (enableEnergy_){
                    thermalHalfTrans_[directionalIsId(elemIdx, outsideElemIdx)] = 0.0;
                    thermalHalfTrans_[directionalIsId(outsideElemIdx, elemIdx)] = 0.0;
                }

                if (updateDiffusivity) {
                    diffusivity_[isId(elemIdx, outsideElemIdx)] = 0.0;
                }
                if (updateDispersivity) {
                    dispersivity_[isId(elemIdx, outsideElemIdx)] = 0.0;
                }
                continue;
            }

            DimVector faceCenterInside;
            DimVector faceCenterOutside;
            DimVector faceAreaNormal;

            typename std::is_same<Grid, Dune::CpGrid>::type isCpGrid;
            computeFaceProperties(intersection,
                                  elemIdx,
                                  insideFaceIdx,
                                  outsideElemIdx,
                                  outsideFaceIdx,
                                  faceCenterInside,
                                  faceCenterOutside,
                                  faceAreaNormal,
                                  isCpGrid);

            Scalar halfTrans1;
            Scalar halfTrans2;

            computeHalfTrans_(halfTrans1,
                              faceAreaNormal,
                              insideFaceIdx,
                              distanceVector_(faceCenterInside,
                                              intersection.indexInInside(),
                                              elemIdx,
                                              axisCentroids),
                              permeability_[elemIdx]);
            computeHalfTrans_(halfTrans2,
                              faceAreaNormal,
                              outsideFaceIdx,
                              distanceVector_(faceCenterOutside,
                                              intersection.indexInOutside(),
                                              outsideElemIdx,
                                              axisCentroids),
                              permeability_[outsideElemIdx]);

            applyNtg_(halfTrans1, insideFaceIdx, elemIdx, ntg);
            applyNtg_(halfTrans2, outsideFaceIdx, outsideElemIdx, ntg);

            // convert half transmissibilities to full face
            // transmissibilities using the harmonic mean
            Scalar trans;
            if (std::abs(halfTrans1) < 1e-30 || std::abs(halfTrans2) < 1e-30)
                // avoid division by zero
                trans = 0.0;
            else
                trans = 1.0 / (1.0/halfTrans1 + 1.0/halfTrans2);

            // apply the full face transmissibility multipliers
            // for the inside ...
            if(!pinchActive){
                if (insideFaceIdx > 3){// top or bottom
                     auto find_layer = [&cartDims](std::size_t cell){
                        cell /= cartDims[0];
                        auto k = cell / cartDims[1];
                        return k;
                    };
                    int kup = find_layer(insideCartElemIdx);
                    int kdown=find_layer(outsideCartElemIdx);
                    // When a grid is a CpGrid with LGRs, insideCartElemIdx coincides with outsideCartElemIdx
                    // for cells on the leaf with the same parent cell on level zero.
                    assert((kup != kdown) || (insideCartElemIdx == outsideCartElemIdx));
                    if(std::abs(kup -kdown) > 1){
                        trans = 0.0;
                    }
                }
            }

            if (useSmallestMultiplier)
            {
                // Currently PINCH(4) is never queries and hence  PINCH(4) == TOPBOT is assumed
                // and in this branch PINCH(5) == ALL holds
                applyAllZMultipliers_(trans, insideFaceIdx, outsideFaceIdx, insideCartElemIdx,
                                      outsideCartElemIdx, transMult, cartDims,
                                      /* pinchTop= */ false);
            }
            else
            {
                applyMultipliers_(trans, insideFaceIdx, insideCartElemIdx, transMult);
                // ... and outside elements
                applyMultipliers_(trans, outsideFaceIdx, outsideCartElemIdx, transMult);
            }

            // apply the region multipliers (cf. the MULTREGT keyword)
            FaceDir::DirEnum faceDir;
            switch (insideFaceIdx) {
            case 0:
            case 1:
                faceDir = FaceDir::XPlus;
                break;

            case 2:
            case 3:
                faceDir = FaceDir::YPlus;
                break;

            case 4:
            case 5:
                faceDir = FaceDir::ZPlus;
                break;

            default:
                throw std::logic_error("Could not determine a face direction");
            }

            trans *= transMult.getRegionMultiplier(insideCartElemIdx,
                                                   outsideCartElemIdx,
                                                   faceDir);

            trans_[isId(elemIdx, outsideElemIdx)] = trans;

            // update the "thermal half transmissibility" for the intersection
            if (enableEnergy_) {

                Scalar halfDiffusivity1;
                Scalar halfDiffusivity2;

                computeHalfDiffusivity_(halfDiffusivity1,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterInside,
                                                        intersection.indexInInside(),
                                                        elemIdx,
                                                        axisCentroids),
                                        1.0);
                computeHalfDiffusivity_(halfDiffusivity2,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterOutside,
                                                        intersection.indexInOutside(),
                                                        outsideElemIdx,
                                                        axisCentroids),
                                        1.0);
                //TODO Add support for multipliers
                thermalHalfTrans_[directionalIsId(elemIdx, outsideElemIdx)] = halfDiffusivity1;
                thermalHalfTrans_[directionalIsId(outsideElemIdx, elemIdx)] = halfDiffusivity2;
           }

            // update the "diffusive half transmissibility" for the intersection
            if (updateDiffusivity) {

                Scalar halfDiffusivity1;
                Scalar halfDiffusivity2;

                computeHalfDiffusivity_(halfDiffusivity1,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterInside,
                                                        intersection.indexInInside(),
                                                        elemIdx,
                                                        axisCentroids),
                                        porosity_[elemIdx]);
                computeHalfDiffusivity_(halfDiffusivity2,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterOutside,
                                                        intersection.indexInOutside(),
                                                        outsideElemIdx,
                                                        axisCentroids),
                                        porosity_[outsideElemIdx]);

                applyNtg_(halfDiffusivity1, insideFaceIdx, elemIdx, ntg);
                applyNtg_(halfDiffusivity2, outsideFaceIdx, outsideElemIdx, ntg);

                //TODO Add support for multipliers
                Scalar diffusivity;
                if (std::abs(halfDiffusivity1) < 1e-30 || std::abs(halfDiffusivity2) < 1e-30)
                    // avoid division by zero
                    diffusivity = 0.0;
                else
                    diffusivity = 1.0 / (1.0/halfDiffusivity1 + 1.0/halfDiffusivity2);


                diffusivity_[isId(elemIdx, outsideElemIdx)] = diffusivity;
           }

           // update the "dispersivity half transmissibility" for the intersection
            if (updateDispersivity) {

                Scalar halfDispersivity1;
                Scalar halfDispersivity2;

                computeHalfDiffusivity_(halfDispersivity1,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterInside,
                                                        intersection.indexInInside(),
                                                        elemIdx,
                                                        axisCentroids),
                                        dispersion_[elemIdx]);
                computeHalfDiffusivity_(halfDispersivity2,
                                        faceAreaNormal,
                                        distanceVector_(faceCenterOutside,
                                                        intersection.indexInOutside(),
                                                        outsideElemIdx,
                                                        axisCentroids),
                                        dispersion_[outsideElemIdx]);

                applyNtg_(halfDispersivity1, insideFaceIdx, elemIdx, ntg);
                applyNtg_(halfDispersivity2, outsideFaceIdx, outsideElemIdx, ntg);

                //TODO Add support for multipliers
                Scalar dispersivity;
                if (std::abs(halfDispersivity1) < 1e-30 || std::abs(halfDispersivity2) < 1e-30)
                    // avoid division by zero
                    dispersivity = 0.0;
                else
                    dispersivity = 1.0 / (1.0/halfDispersivity1 + 1.0/halfDispersivity2);


                dispersivity_[isId(elemIdx, outsideElemIdx)] = dispersivity;
           }
        }
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

    this->applyEditNncToGridTrans_(globalToLocal);
    this->applyNncToGridTrans_(globalToLocal);
    this->applyEditNncrToGridTrans_(globalToLocal);

    if (applyNncMultregT) {
        this->applyNncMultreg_(globalToLocal);
    }

    // Remove very small non-neighbouring transmissibilities.
    this->removeSmallNonCartesianTransmissibilities_();
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
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
        const std::vector<double>& permxData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp, "PERMX", numElem);

        std::vector<double> permyData;
        if (fp.has_double("PERMY"))
            permyData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMY", numElem);
        else
            permyData = permxData;

        std::vector<double> permzData;
        if (fp.has_double("PERMZ"))
            permzData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMZ", numElem);
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
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
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
        const std::vector<double>& permxData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMX", numElem);

        std::vector<double> permyData;
        if (fp.has_double("PERMY"))
            permyData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMY", numElem);
        else
            permyData = permxData;

        std::vector<double> permzData;
        if (fp.has_double("PERMZ"))
            permzData = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PERMZ", numElem);
        else
            permzData = permxData;

        for (std::size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
            permeability_[dofIdx] = 0.0;
            std::size_t inputDofIdx = map(dofIdx);
            permeability_[dofIdx][0][0] = permxData[inputDofIdx];
            permeability_[dofIdx][1][1] = permyData[inputDofIdx];
            permeability_[dofIdx][2][2] = permzData[inputDofIdx];
        }

        // for now we don't care about non-diagonal entries

    }
    else
        throw std::logic_error("Can't read the intrinsic permeability from the ecl state. "
                               "(The PERM{X,Y,Z} keywords are missing)");
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractPorosity_()
{
    // read the intrinsic porosity from the eclState. Note that all arrays
    // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
    // simulation grid might remove a few elements. (e.g. because it is distributed
    // over several processes.)
    unsigned numElem = gridView_.size(/*codim=*/0);
    const auto& fp = eclState_.fieldProps();
    if (fp.has_double("PORO")) {
        porosity_ = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"PORO", numElem);
    }
    else
        throw std::logic_error("Can't read the porosityfrom the ecl state. "
                               "(The PORO keywords are missing)");
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
extractDispersion_()
{
    if (!enableDispersivity_) {
        throw std::runtime_error("Dispersion disabled at compile time, but the deck "
                                 "contains the DISPERC keyword.");
    }
    unsigned numElem = gridView_.size(/*codim=*/0);
    const auto& fp = eclState_.fieldProps();
    dispersion_ = this-> lookUpData_.assignFieldPropsDoubleOnLeaf(fp,"DISPERC", numElem);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
removeSmallNonCartesianTransmissibilities_()
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    for (auto&& trans: trans_) {
        if (trans.second < transmissibilityThreshold_) {
            const auto& id = trans.first;
            const auto& elements = isIdReverse(id);
            int gc1 = std::min(cartMapper_.cartesianIndex(elements.first), cartMapper_.cartesianIndex(elements.second));
            int gc2 = std::max(cartMapper_.cartesianIndex(elements.first), cartMapper_.cartesianIndex(elements.second));

            // only adjust the NNCs
            // When LGRs, all neighbors in the LGR are cartesian neighbours on the level grid representing the LGR.
            // When elements on the leaf grid view have the same parent cell, gc1 and gc2 coincide.
            if (gc2 - gc1 == 1 || gc2 - gc1 == cartDims[0] || gc2 - gc1 == cartDims[0]*cartDims[1] || gc2 - gc1 == 0)
                continue;

            //remove transmissibilities less than the threshold (by default 1e-6 in the deck's unit system)
            trans.second = 0.0;
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper, Scalar>::
applyAllZMultipliers_(Scalar& trans,
                      unsigned insideFaceIdx,
                      unsigned outsideFaceIdx,
                      unsigned insideCartElemIdx,
                      unsigned outsideCartElemIdx,
                      const TransMult& transMult,
                      const std::array<int, dimWorld>& cartDims,
                      bool pinchTop)
{
    if(grid_.maxLevel()> 0) {
                OPM_THROW(std::invalid_argument, "MULTZ not support with LGRS, yet.");
    }
    if (insideFaceIdx > 3) { // top or or bottom
        assert(insideFaceIdx==5); // as insideCartElemIdx < outsideCartElemIdx holds for the Z column
        // For CpGrid with LGRs, insideCartElemIdx == outsideCartElemIdx when cells on the leaf have the same parent cell on level zero.
        assert(outsideCartElemIdx >= insideCartElemIdx);
        unsigned lastCartElemIdx;
        if (outsideCartElemIdx == insideCartElemIdx) {
            lastCartElemIdx = outsideCartElemIdx;
        }
        else {
            lastCartElemIdx = outsideCartElemIdx - cartDims[0]*cartDims[1];
        }
        // Last multiplier using (Z+)*(Z-)
        Scalar mult = transMult.getMultiplier(lastCartElemIdx , FaceDir::ZPlus) *
            transMult.getMultiplier(outsideCartElemIdx , FaceDir::ZMinus);

        if ( !pinchTop )
        {
            // pick the smallest multiplier using (Z+)*(Z-) while looking down
            // the pillar until reaching the other end of the connection
            for(auto cartElemIdx = insideCartElemIdx; cartElemIdx < lastCartElemIdx;)
            {
                auto multiplier = transMult.getMultiplier(cartElemIdx, FaceDir::ZPlus);
                cartElemIdx += cartDims[0]*cartDims[1];
                multiplier *= transMult.getMultiplier(cartElemIdx, FaceDir::ZMinus);
                mult = std::min(mult, multiplier);
            }
        }

        trans *= mult;
    }
    else
    {
        applyMultipliers_(trans, insideFaceIdx, insideCartElemIdx, transMult);
        applyMultipliers_(trans, outsideFaceIdx, outsideCartElemIdx, transMult);
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
updateFromEclState_(bool global)
{
    const FieldPropsManager* fp =
        (global) ? &(eclState_.fieldProps()) :
        &(eclState_.globalFieldProps());

    std::array<bool,3> is_tran {fp->tran_active("TRANX"),
                                fp->tran_active("TRANY"),
                                fp->tran_active("TRANZ")};

    if( !(is_tran[0] ||is_tran[1] || is_tran[2]) )
    {
        // Skip unneeded expensive traversals
        return;
    }

    std::array<std::string, 3> keywords {"TRANX", "TRANY", "TRANZ"};
    std::array<std::vector<double>,3> trans = createTransmissibilityArrays_(is_tran);
    auto key = keywords.begin();
    auto perform = is_tran.begin();

    for (auto it = trans.begin(); it != trans.end(); ++it, ++key, ++perform)
    {
        if(*perform) {
            if(grid_.maxLevel()>0) {
                OPM_THROW(std::invalid_argument, "Calculations on TRANX/TRANY/TRANZ arrays are not support with LGRS, yet.");
            }
            fp->apply_tran(*key, *it);
        }
    }

    resetTransmissibilityFromArrays_(is_tran, trans);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
std::array<std::vector<double>,3>
EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
createTransmissibilityArrays_(const std::array<bool,3>& is_tran)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    auto numElem = gridView_.size(/*codim=*/0);
    std::array<std::vector<double>,3> trans =
        { std::vector<double>(is_tran[0] ? numElem : 0, 0),
          std::vector<double>(is_tran[1] ? numElem : 0, 0),
          std::vector<double>(is_tran[2] ? numElem : 0, 0)};

    // compute the transmissibilities for all intersections
    for (const auto& elem : elements(gridView_)) {
        for (const auto& intersection : intersections(gridView_, elem)) {
            // store intersection, this might be costly
            if (!intersection.neighbor())
                continue; // intersection is on the domain boundary

            // In the EclState TRANX[c1] is transmissibility in X+
            // direction. Ordering of compressed (c1,c2) and cartesian index
            // (gc1, gc2) is coherent (c1 < c2 <=> gc1 < gc2). This also
            // holds for the global grid. While distributing changes the
            // order of the local indices, the transmissibilities are still
            // stored at the cell with the lower global cartesian index as
            // the fieldprops are communicated by the grid.
            unsigned c1 = elemMapper.index(intersection.inside());
            unsigned c2 = elemMapper.index(intersection.outside());
            int gc1 = cartMapper_.cartesianIndex(c1);
            int gc2 = cartMapper_.cartesianIndex(c2);

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            auto isID = isId(c1, c2);

            if (gc2 - gc1 == 1 && cartDims[0] > 1) {
                if (is_tran[0])
                    // set simulator internal transmissibilities to values from inputTranx
                     trans[0][c1] = trans_[isID];
            }
            else if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                if (is_tran[1])
                    // set simulator internal transmissibilities to values from inputTrany
                     trans[1][c1] = trans_[isID];
            }
            else if (gc2 - gc1 == cartDims[0]*cartDims[1]) {
                if (is_tran[2])
                    // set simulator internal transmissibilities to values from inputTranz
                     trans[2][c1] = trans_[isID];
            }
            //else.. We don't support modification of NNC at the moment.
        }
    }

    return trans;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
resetTransmissibilityFromArrays_(const std::array<bool,3>& is_tran,
                                 const std::array<std::vector<double>,3>& trans)
{
    const auto& cartDims = cartMapper_.cartesianDimensions();
    ElementMapper elemMapper(gridView_, Dune::mcmgElementLayout());

    // compute the transmissibilities for all intersections
    for (const auto& elem : elements(gridView_)) {
        for (const auto& intersection : intersections(gridView_, elem)) {
            if (!intersection.neighbor())
                continue; // intersection is on the domain boundary

            // In the EclState TRANX[c1] is transmissibility in X+
            // direction. Ordering of compressed (c1,c2) and cartesian index
            // (gc1, gc2) is coherent (c1 < c2 <=> gc1 < gc2). This also
            // holds for the global grid. While distributing changes the
            // order of the local indices, the transmissibilities are still
            // stored at the cell with the lower global cartesian index as
            // the fieldprops are communicated by the grid.
            unsigned c1 = elemMapper.index(intersection.inside());
            unsigned c2 = elemMapper.index(intersection.outside());
            int gc1 = cartMapper_.cartesianIndex(c1);
            int gc2 = cartMapper_.cartesianIndex(c2);
            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            auto isID = isId(c1, c2);

            if (gc2 - gc1 == 1 && cartDims[0] > 1) {
                if (is_tran[0])
                    // set simulator internal transmissibilities to values from inputTranx
                    trans_[isID] = trans[0][c1];
            }
            else if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                if (is_tran[1])
                    // set simulator internal transmissibilities to values from inputTrany
                    trans_[isID] = trans[1][c1];
            }
            else if (gc2 - gc1 == cartDims[0]*cartDims[1]) {
                if (is_tran[2])
                    // set simulator internal transmissibilities to values from inputTranz
                    trans_[isID] = trans[2][c1];
            }
            //else.. We don't support modification of NNC at the moment.
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeFaceProperties(const Intersection& intersection,
                      const int,
                      const int,
                      const int,
                      const int,
                      DimVector& faceCenterInside,
                      DimVector& faceCenterOutside,
                      DimVector& faceAreaNormal,
                      /*isCpGrid=*/std::false_type) const
{
    // default implementation for DUNE grids
    const auto& geometry = intersection.geometry();
    faceCenterInside = geometry.center();
    faceCenterOutside = faceCenterInside;

    faceAreaNormal = intersection.centerUnitOuterNormal();
    faceAreaNormal *= geometry.volume();
}


template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
template<class Intersection>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeFaceProperties(const Intersection& intersection,
                      const int insideElemIdx,
                      const int insideFaceIdx,
                      const int outsideElemIdx,
                      const int outsideFaceIdx,
                      DimVector& faceCenterInside,
                      DimVector& faceCenterOutside,
                      DimVector& faceAreaNormal,
                      /*isCpGrid=*/std::true_type) const
{
    int faceIdx = intersection.id();
    faceCenterInside = grid_.faceCenterEcl(insideElemIdx, insideFaceIdx);
    faceCenterOutside = grid_.faceCenterEcl(outsideElemIdx, outsideFaceIdx);
    faceAreaNormal = grid_.faceAreaNormalEcl(faceIdx);
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void
EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
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

        if (low > high)
            std::swap(low, high);

        if (low == -1 && high == -1)
            // Silently discard as it is not between active cells
            continue;

        if (low == -1 || high == -1) {
            // Discard the NNC if it is between active cell and inactive cell
            std::ostringstream sstr;
            sstr << "NNC between active and inactive cells ("
                 << low << " -> " << high << ") with globalcell is (" << c1 << "->" << c2 <<")";
            OpmLog::warning(sstr.str());
            continue;
        }

        {
            auto candidate = trans_.find(isId(low, high));
            if (candidate != trans_.end()) {
                // NNC is represented by the grid and might be a neighboring connection
                // In this case the transmissibilty is added to the value already
                // set or computed.
                candidate->second += nncEntry.trans;
            }
        }
        // if (enableEnergy_) {
        //     auto candidate = thermalHalfTrans_.find(directionalIsId(low, high));
        //     if (candidate != trans_.end()) {
        //         // NNC is represented by the grid and might be a neighboring connection
        //         // In this case the transmissibilty is added to the value already
        //         // set or computed.
        //         candidate->second += nncEntry.transEnergy1;
        //     }
        //     auto candidate = thermalHalfTrans_.find(directionalIsId(high, low));
        //     if (candidate != trans_.end()) {
        //         // NNC is represented by the grid and might be a neighboring connection
        //         // In this case the transmissibilty is added to the value already
        //         // set or computed.
        //         candidate->second += nncEntry.transEnergy2;
        //     }
        // }
        // if (enableDiffusivity_) {
        //     auto candidate = diffusivity_.find(isId(low, high));
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
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal)
{
    const auto& input = eclState_.getInputNNC();
    applyEditNncToGridTransHelper_(globalToLocal, "EDITNNC",
                                   input.edit(),
                                   [&input](const NNCdata& nnc){
                                       return input.edit_location(nnc);},
                                   // Multiply transmissibility with EDITNNC value
                                   [](double& trans, const double& rhs){ trans *= rhs;});
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncrToGridTrans_(const std::unordered_map<std::size_t,int>& globalToLocal)
{
    const auto& input = eclState_.getInputNNC();
    applyEditNncToGridTransHelper_(globalToLocal, "EDITNNCR",
                                   input.editr(),
                                   [&input](const NNCdata& nnc){
                                       return input.editr_location(nnc);},
                                   // Replace Transmissibility with EDITNNCR value
                                   [](double& trans, const double& rhs){ trans = rhs;});
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyEditNncToGridTransHelper_(const std::unordered_map<std::size_t,int>& globalToLocal,
                               const std::string& keyword,
                               const std::vector<NNCdata>& nncs,
                               const std::function<KeywordLocation(const NNCdata&)>& getLocation,
                               const std::function<void(double&, const double&)>& apply)
{
    if (nncs.empty())
        return;
    const auto& cartDims = cartMapper_.cartesianDimensions();

    auto format_ijk = [&cartDims](std::size_t cell) -> std::string {
        auto i = cell % cartDims[0]; cell /= cartDims[0];
        auto j = cell % cartDims[1];
        auto k = cell / cartDims[1];

        return fmt::format("({},{},{})", i + 1,j + 1,k + 1);
    };

    auto print_warning = [&format_ijk, &getLocation, &keyword] (const NNCdata& nnc) {
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
            print_warning(*nnc);
            ++nnc;
            warning_count++;
            continue;
        }

        auto low = lowIt->second, high = highIt->second;

        if (low > high)
            std::swap(low, high);

        auto candidate = trans_.find(isId(low, high));
        if (candidate == trans_.end()) {
            print_warning(*nnc);
            ++nnc;
            warning_count++;
        }
        else {
            // NNC exists
            while (nnc!= end && c1==nnc->cell1 && c2==nnc->cell2) {
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
EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
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
    for (const auto& nncList : { &NNC::input, &NNC::editr }) {
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

            auto candidate = this->trans_.find(isId(low, high));
            if (candidate != this->trans_.end()) {
                candidate->second *= transMult.getRegionMultiplierNNC(c1, c2);
            }
        }
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeHalfTrans_(Scalar& halfTrans,
                  const DimVector& areaNormal,
                  int faceIdx, // in the reference element that contains the intersection
                  const DimVector& distance,
                  const DimMatrix& perm) const
{
    assert(faceIdx >= 0);
    unsigned dimIdx = faceIdx/2;
    assert(dimIdx < dimWorld);
    halfTrans = perm[dimIdx][dimIdx];

    Scalar val = 0;
    for (unsigned i = 0; i < areaNormal.size(); ++i)
        val += areaNormal[i]*distance[i];

    halfTrans *= std::abs(val);
    halfTrans /= distance.two_norm2();
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
computeHalfDiffusivity_(Scalar& halfDiff,
                        const DimVector& areaNormal,
                        const DimVector& distance,
                        const Scalar& poro) const
{
    halfDiff = poro;
    Scalar val = 0;
    for (unsigned i = 0; i < areaNormal.size(); ++i)
        val += areaNormal[i]*distance[i];

    halfDiff *= std::abs(val);
    halfDiff /= distance.two_norm2();
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
typename EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::DimVector
EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
distanceVector_(const DimVector& center,
                int faceIdx, // in the reference element that contains the intersection
                unsigned elemIdx,
                const std::array<std::vector<DimVector>, dimWorld>& axisCentroids) const
{
    assert(faceIdx >= 0);
    unsigned dimIdx = faceIdx/2;
    assert(dimIdx < dimWorld);
    DimVector x = center;
    x -= axisCentroids[dimIdx][elemIdx];

    return x;
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyMultipliers_(Scalar& trans,
                  unsigned faceIdx,
                  unsigned cartElemIdx,
                  const TransMult& transMult) const
{
    // apply multiplyer for the transmissibility of the face. (the
    // face index is the index of the reference-element face which
    // contains the intersection of interest.)
    switch (faceIdx) {
    case 0: // left
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::XMinus);
        break;
    case 1: // right
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::XPlus);
        break;

    case 2: // front
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::YMinus);
        break;
    case 3: // back
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::YPlus);
        break;

    case 4: // bottom
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::ZMinus);
        break;
    case 5: // top
        trans *= transMult.getMultiplier(cartElemIdx, FaceDir::ZPlus);
        break;
    }
}

template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
void EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>::
applyNtg_(Scalar& trans,
          unsigned faceIdx,
          unsigned elemIdx,
          const std::vector<double>& ntg) const
{
    // apply multiplyer for the transmissibility of the face. (the
    // face index is the index of the reference-element face which
    // contains the intersection of interest.)
    switch (faceIdx) {
    case 0: // left
        trans *= ntg[elemIdx];
        break;
    case 1: // right
        trans *= ntg[elemIdx];
        break;

    case 2: // front
        trans *= ntg[elemIdx];
        break;
    case 3: // back
        trans *= ntg[elemIdx];
        break;

        // NTG does not apply to top and bottom faces
    }
}


} // namespace Opm
#endif
