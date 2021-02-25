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
/*!
 * \file
 *
 * \copydoc Opm::EclTransmissibility
 */
#ifndef EWOMS_ECL_TRANSMISSIBILITY_HH
#define EWOMS_ECL_TRANSMISSIBILITY_HH

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/common/multiphasebaseproperties.hh>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>


#include <opm/grid/CpGrid.hpp>

#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <fmt/format.h>
#include <array>
#include <vector>
#include <unordered_map>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the transmissibilites for grid faces according to the
 *        Eclipse Technical Description.
 */
template <class TypeTag>
class EclTransmissibility
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Intersection = typename GridView::Intersection;

    static const bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    static const bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();

    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    static const unsigned elemIdxShift = 32; // bits

public:

    EclTransmissibility(const Vanguard& vanguard)
        : vanguard_(vanguard)
    {
        const Opm::UnitSystem& unitSystem = vanguard_.eclState().getDeckUnitSystem();
        transmissibilityThreshold_  = unitSystem.parse("Transmissibility").getSIScaling() * 1e-6;
    }

    /*!
     * \brief Actually compute the transmissibilty over a face as a pre-compute step.
     *
     * This code actually uses the direction specific "centroids" of
     * each element. These "centroids" are _not_ the identical
     * barycenter of the element, but the middle of the centers of the
     * faces of the logical Cartesian cells, i.e., the centers of the
     * faces of the reference elements. We do things this way because
     * the barycenter of the element can be located outside of the
     * element for sufficiently "ugly" (i.e., thin and "non-flat")
     * elements which in turn leads to quite wrong
     * permeabilities. This approach is probably not always correct
     * either but at least it seems to be much better.
     */
    void finishInit()
    { update(true); }


    /*!
     * \brief Compute all transmissibilities
     *
     * \param global If true, update is called on all processes
     * Also, this updates the "thermal half transmissibilities" if energy is enabled.
     */
    void update(bool global)
    {
        const auto& gridView = vanguard_.gridView();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& eclState = vanguard_.eclState();
        const auto& cartDims = cartMapper.cartesianDimensions();
        auto& transMult = eclState.getTransMult();
        const auto& comm = vanguard_.gridView().comm();
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

        // get the ntg values, the ntg values are modified for the cells merged with minpv
        const std::vector<double>& ntg = eclState.fieldProps().get_double("NTG");
        const bool updateDiffusivity = eclState.getSimulationConfig().isDiffusive() && enableDiffusion;
        unsigned numElements = elemMapper.size();

        extractPermeability_();

        // calculate the axis specific centroids of all elements
        std::array<std::vector<DimVector>, dimWorld> axisCentroids;

        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            axisCentroids[dimIdx].resize(numElements);

        const std::vector<double>& centroids = vanguard_.cellCentroids();

        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        size_t centroidIdx = 0;
        for (; elemIt != elemEndIt; ++elemIt, ++centroidIdx) {
            const auto& elem = *elemIt;
            unsigned elemIdx = elemMapper.index(elem);

            // compute the axis specific "centroids" used for the transmissibilities. for
            // consistency with the flow simulator, we use the element centers as
            // computed by opm-parser's Opm::EclipseGrid class for all axes.
            std::array<double, 3> centroid;
            if (vanguard_.gridView().comm().rank() == 0) {
                const auto& eclGrid = eclState.getInputGrid();
                unsigned cartesianCellIdx = cartMapper.cartesianIndex(elemIdx);
                centroid = eclGrid.getCellCenter(cartesianCellIdx);
            } else
                std::copy(centroids.begin() + centroidIdx * dimWorld,
                          centroids.begin() + (centroidIdx + 1) * dimWorld,
                          centroid.begin());

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
        if (enableEnergy) {
            thermalHalfTrans_->clear();
            thermalHalfTrans_->reserve(numElements*6*1.05);

            thermalHalfTransBoundary_.clear();
        }

        // if diffusion is enabled, let's do the same for the "diffusivity"
        if (updateDiffusivity) {
            diffusivity_->clear();
            diffusivity_->reserve(numElements*3*1.05);
            extractPorosity_();
        }

        // The MULTZ needs special case if the option is ALL
        // Then the smallest multiplier is applied.
        // Default is to apply the top and bottom multiplier
        bool useSmallestMultiplier;
        if (comm.rank() == 0) {
            const auto& eclGrid = eclState.getInputGrid();
            useSmallestMultiplier = eclGrid.getMultzOption() == Opm::PinchMode::ModeEnum::ALL;
        }
        if (global && comm.size() > 1) {
            comm.broadcast(&useSmallestMultiplier, 1, 0);
        }

        // compute the transmissibilities for all intersections
        elemIt = gridView.template begin</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            unsigned elemIdx = elemMapper.index(elem);

            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
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
                    transBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] = transBoundaryIs;

                    // for boundary intersections we also need to compute the thermal
                    // half transmissibilities
                    if (enableEnergy) {
                        const auto& n = intersection.centerUnitOuterNormal();
                        const auto& inPos = elem.geometry().center();
                        const auto& outPos = intersection.geometry().center();
                        const auto& d = outPos - inPos;

                        // eWoms expects fluxes to be area specific, i.e. we must *not*
                        // the transmissibility with the face area here
                        Scalar thermalHalfTrans = std::abs(n*d)/(d*d);

                        thermalHalfTransBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] =
                            thermalHalfTrans;
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

                // update the "thermal half transmissibility" for the intersection
                if (enableEnergy) {
                    const auto& n = intersection.centerUnitOuterNormal();
                    Scalar A = intersection.geometry().volume();

                    const auto& inPos = elem.geometry().center();
                    const auto& outPos = intersection.geometry().center();
                    const auto& d = outPos - inPos;

                    (*thermalHalfTrans_)[directionalIsId_(elemIdx, outsideElemIdx)] =
                        A * (n*d)/(d*d);
                }

                unsigned insideCartElemIdx = cartMapper.cartesianIndex(elemIdx);
                unsigned outsideCartElemIdx = cartMapper.cartesianIndex(outsideElemIdx);

                // we only need to calculate a face's transmissibility
                // once...
                if (insideCartElemIdx > outsideCartElemIdx)
                    continue;

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                int insideFaceIdx  = intersection.indexInInside();
                int outsideFaceIdx = intersection.indexInOutside();

                if (insideFaceIdx == -1) {
                    // NNC. Set zero transmissibility, as it will be
                    // *added to* by applyNncToGridTrans_() later.
                    assert(outsideFaceIdx == -1);
                    trans_[isId_(elemIdx, outsideElemIdx)] = 0.0;
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
                Opm::FaceDir::DirEnum faceDir;
                switch (insideFaceIdx) {
                case 0:
                case 1:
                    faceDir = Opm::FaceDir::XPlus;
                    break;

                case 2:
                case 3:
                    faceDir = Opm::FaceDir::YPlus;
                    break;

                case 4:
                case 5:
                    faceDir = Opm::FaceDir::ZPlus;
                    break;

                default:
                    throw std::logic_error("Could not determine a face direction");
                }

                trans *= transMult.getRegionMultiplier(insideCartElemIdx,
                                                       outsideCartElemIdx,
                                                       faceDir);

                trans_[isId_(elemIdx, outsideElemIdx)] = trans;

                // update the "thermal half transmissibility" for the intersection
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


                    (*diffusivity_)[isId_(elemIdx, outsideElemIdx)] = diffusivity;
               }
            }
        }

        // potentially overwrite and/or modify  transmissibilities based on input from deck
        updateFromEclState_(global);

        // Create mapping from global to local index
        const size_t cartesianSize = cartMapper.cartesianSize();
        // reserve memory
        std::vector<int> globalToLocal(cartesianSize, -1);

        // loop over all elements (global grid) and store Cartesian index
        elemIt = vanguard_.grid().leafGridView().template begin<0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = elemMapper.index(*elemIt);
            int cartElemIdx = vanguard_.cartesianIndexMapper().cartesianIndex(elemIdx);
            globalToLocal[cartElemIdx] = elemIdx;
        }
        applyEditNncToGridTrans_(globalToLocal);
        applyNncToGridTrans_(globalToLocal);

        //remove very small non-neighbouring transmissibilities
        removeSmallNonCartesianTransmissibilities_();
    }

    /*!
     * \brief Return the permeability for an element.
     */
    const DimMatrix& permeability(unsigned elemIdx) const
    { return permeability_[elemIdx]; }

    /*!
     * \brief Return the transmissibility for the intersection between two elements.
     */
    Scalar transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
    { return trans_.at(isId_(elemIdx1, elemIdx2)); }

    /*!
     * \brief Return the transmissibility for a given boundary segment.
     */
    Scalar transmissibilityBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
    { return transBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx)); }

    /*!
     * \brief Return the thermal "half transmissibility" for the intersection between two
     *        elements.
     *
     * The "half transmissibility" features all sub-expressions of the "thermal
     * transmissibility" which can be precomputed, i.e. they are not dependent on the
     * current solution:
     *
     * H_t = A * (n*d)/(d*d);
     *
     * where A is the area of the intersection between the inside and outside elements, n
     * is the outer unit normal, and d is the distance between the center of the inside
     * cell and the center of the intersection.
     */
    Scalar thermalHalfTrans(unsigned insideElemIdx, unsigned outsideElemIdx) const
    { return thermalHalfTrans_->at(directionalIsId_(insideElemIdx, outsideElemIdx)); }

    Scalar thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const
    { return thermalHalfTransBoundary_.at(std::make_pair(insideElemIdx, boundaryFaceIdx)); }

    /*!
     * \brief Return the diffusivity for the intersection between two elements.
     */
    Scalar diffusivity(unsigned elemIdx1, unsigned elemIdx2) const
    {
        if(diffusivity_->empty())
            return 0.0;

        return diffusivity_->at(isId_(elemIdx1, elemIdx2));

    }

private:

    void removeSmallNonCartesianTransmissibilities_()
    {
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        for (auto&& trans: trans_) {
            if (trans.second < transmissibilityThreshold_) {
                const auto& id = trans.first;
                const auto& elements = isIdReverse_(id);
                int gc1 = std::min(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
                int gc2 = std::max(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));

                // only adjust the NNCs
                if (gc2 - gc1 == 1 || gc2 - gc1 == cartDims[0] || gc2 - gc1 == cartDims[0]*cartDims[1])
                    continue;

                //remove transmissibilities less than the threshold (by default 1e-6 in the deck's unit system)
                trans.second = 0.0;
            }
        }
    }

    /// \brief Apply the Multipliers for the case PINCH(4)==TOPBOT
    ///
    /// \param pinchTop Whether PINCH(5) is TOP, otherwise ALL is assumed.
    void applyAllZMultipliers_(Scalar& trans,
                               unsigned insideFaceIdx,
                               unsigned outsideFaceIdx,
                               unsigned insideCartElemIdx,
                               unsigned outsideCartElemIdx,
                               const Opm::TransMult& transMult,
                               const std::array<int, dimWorld>& cartDims,
                               bool pinchTop)
    {
        if (insideFaceIdx > 3) { // top or or bottom
            assert(insideFaceIdx==5); // as insideCartElemIdx < outsideCartElemIdx holds for the Z column
            assert(outsideCartElemIdx > insideCartElemIdx);
            auto lastCartElemIdx = outsideCartElemIdx - cartDims[0]*cartDims[1];
            // Last multiplier
            Scalar mult = transMult.getMultiplier(lastCartElemIdx , Opm::FaceDir::ZPlus);

            if ( !pinchTop )
            {
                // pick the smallest multiplier for Z+ while looking down the pillar until reaching the other end of the connection
                // While Z- is not all used here.
                for(auto cartElemIdx = insideCartElemIdx; cartElemIdx < lastCartElemIdx;
                    cartElemIdx += cartDims[0]*cartDims[1])
                {
                    mult = std::min(mult, transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZPlus));
                }
            }

            trans *= mult;
            applyMultipliers_(trans, outsideFaceIdx, outsideCartElemIdx, transMult);
        }
        else
        {
            applyMultipliers_(trans, insideFaceIdx, insideCartElemIdx, transMult);
            applyMultipliers_(trans, outsideFaceIdx, outsideCartElemIdx, transMult);
        }
    }

    void updateFromEclState_(bool global)
    {
        const FieldPropsManager* fp =
            (global) ? &(vanguard_.eclState().fieldProps()) :
            &(vanguard_.eclState().globalFieldProps());

        std::array<bool,3> is_tran {fp->tran_active("TRANX"),
                                    fp->tran_active("TRANY"),
                                    fp->tran_active("TRANZ")};

        if( !(is_tran[0] ||is_tran[1] || is_tran[2]) )
        {
            // Skip unneeded expensive traversals
            return;
        }

        std::array<std::string, 3> keywords {"TRANX", "TRANY", "TRANZ"};
        std::array<std::vector<double>,3> trans = createTransmissibilityArrays(is_tran);
        auto key = keywords.begin();
        auto perform = is_tran.begin();

        for (auto it = trans.begin(); it != trans.end(); ++it, ++key, ++perform)
        {
            if(perform)
                fp->apply_tran(*key, *it);
        }

        resetTransmissibilityFromArrays(is_tran, trans);
    }

    /// \brief overwrites calculated transmissibilities
    ///
    /// \param is_tran Whether TRAN{XYZ} have been modified.
    /// \param trans Arrays with modified transmissibilities TRAN{XYZ}
    void resetTransmissibilityFromArrays(const std::array<bool,3>& is_tran,
                                         const std::array<std::vector<double>,3>& trans)
    {
        const auto& gridView = vanguard_.gridView();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

        // compute the transmissibilities for all intersections
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;
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
                int gc1 = cartMapper.cartesianIndex(c1);
                int gc2 = cartMapper.cartesianIndex(c2);
                if (gc1 > gc2)
                    continue; // we only need to handle each connection once, thank you.

                auto isId = isId_(c1, c2);

                if (gc2 - gc1 == 1 && cartDims[0] > 1) {
                    if (is_tran[0])
                        // set simulator internal transmissibilities to values from inputTranx
                        trans_[isId] = trans[0][c1];
                }
                else if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                    if (is_tran[1])
                        // set simulator internal transmissibilities to values from inputTrany
                        trans_[isId] = trans[1][c1];
                }
                else if (gc2 - gc1 == cartDims[0]*cartDims[1]) {
                    if (is_tran[2])
                        // set simulator internal transmissibilities to values from inputTranz
                        trans_[isId] = trans[2][c1];
                }
                //else.. We don't support modification of NNC at the moment.
            }
        }
    }

    /// \brief Creates TRANS{XYZ} arrays for modification by FieldProps data
    ///
    /// \param is_tran Whether TRAN{XYZ} will be modified. If entry is false the array will be empty
    /// \returns an array of vector (TRANX, TRANY, TRANZ}
    std::array<std::vector<double>,3>
    createTransmissibilityArrays(const std::array<bool,3>& is_tran)
    {
        const auto& gridView = vanguard_.gridView();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

        auto numElem = vanguard_.gridView().size(/*codim=*/0);
        std::array<std::vector<double>,3> trans =
            { std::vector<double>(is_tran[0] ? numElem : 0, 0),
              std::vector<double>(is_tran[1] ? numElem : 0, 0),
              std::vector<double>(is_tran[2] ? numElem : 0, 0)};

        // compute the transmissibilities for all intersections
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;
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
                int gc1 = cartMapper.cartesianIndex(c1);
                int gc2 = cartMapper.cartesianIndex(c2);

                if (gc1 > gc2)
                    continue; // we only need to handle each connection once, thank you.

                auto isId = isId_(c1, c2);

                if (gc2 - gc1 == 1 && cartDims[0] > 1) {
                    if (is_tran[0])
                        // set simulator internal transmissibilities to values from inputTranx
                         trans[0][c1] = trans_[isId];
                }
                else if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                    if (is_tran[1])
                        // set simulator internal transmissibilities to values from inputTrany
                         trans[1][c1] = trans_[isId];
                }
                else if (gc2 - gc1 == cartDims[0]*cartDims[1]) {
                    if (is_tran[2])
                        // set simulator internal transmissibilities to values from inputTranz
                         trans[2][c1] = trans_[isId];
                }
                //else.. We don't support modification of NNC at the moment.
            }
        }

        return trans;
    }


    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
                               const int insideElemIdx OPM_UNUSED,
                               const int insideFaceIdx OPM_UNUSED,
                               const int outsideElemIdx OPM_UNUSED,
                               const int outsideFaceIdx OPM_UNUSED,
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

    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
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
        faceCenterInside = vanguard_.grid().faceCenterEcl(insideElemIdx, insideFaceIdx);
        faceCenterOutside = vanguard_.grid().faceCenterEcl(outsideElemIdx, outsideFaceIdx);
        faceAreaNormal = vanguard_.grid().faceAreaNormalEcl(faceIdx);
    }

    /*
     * \brief Applies additional transmissibilities specified via NNC keyword.
     *
     * Applies only those NNC that are actually represented by the grid. These may
     * NNCs due to faults or NNCs that are actually neighbours. In both cases that
     * specified transmissibilities (scaled by EDITNNC) will be added to the already
     * existing models.
     *
     * \param cartesianToCompressed Vector containing the compressed index (or -1 for inactive
     *                              cells) as the element at the cartesian index.
     * \return Two vector of NNCs (scaled by EDITNNC). The first one are the NNCs that have been applied
     *         and the second the NNCs not resembled by faces of the grid. NNCs specified for
     *         inactive cells are omitted in these vectors.
     */
    std::tuple<std::vector<Opm::NNCdata>, std::vector<Opm::NNCdata> >
    applyNncToGridTrans_(const std::vector<int>& cartesianToCompressed)
    {
        // First scale NNCs with EDITNNC.
        std::vector<Opm::NNCdata> unprocessedNnc;
        std::vector<Opm::NNCdata> processedNnc;
        const auto& nnc_input = vanguard_.eclState().getInputNNC().input();
        if (nnc_input.empty())
            return make_tuple(processedNnc, unprocessedNnc);

        for (const auto& nncEntry : nnc_input) {
            auto c1 = nncEntry.cell1;
            auto c2 = nncEntry.cell2;
            auto low = cartesianToCompressed[c1];
            auto high = cartesianToCompressed[c2];

            if (low > high)
                std::swap(low, high);

            if (low == -1 && high == -1)
                // Silently discard as it is not between active cells
                continue;

            if (low == -1 || high == -1) {
                // Discard the NNC if it is between active cell and inactive cell
                std::ostringstream sstr;
                sstr << "NNC between active and inactive cells ("
                     << low << " -> " << high << ")";
                Opm::OpmLog::warning(sstr.str());
                continue;
            }

            auto candidate = trans_.find(isId_(low, high));

            if (candidate == trans_.end())
                // This NNC is not resembled by the grid. Save it for later
                // processing with local cell values
                unprocessedNnc.push_back(nncEntry);
            else {
                // NNC is represented by the grid and might be a neighboring connection
                // In this case the transmissibilty is added to the value already
                // set or computed.
                candidate->second += nncEntry.trans;
                processedNnc.push_back(nncEntry);
            }
        }
        return make_tuple(processedNnc, unprocessedNnc);
    }


    /// \brief Multiplies the grid transmissibilities according to EDITNNC.
    void applyEditNncToGridTrans_(const std::vector<int>& globalToLocal)
    {
        const auto& nnc_input = vanguard_.eclState().getInputNNC();
        const auto& editNnc = nnc_input.edit();
        if (editNnc.empty())
            return;
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();

        auto format_ijk = [&cartDims](std::size_t cell) -> std::string {
            auto i = cell % cartDims[0]; cell /= cartDims[0];
            auto j = cell % cartDims[1];
            auto k = cell / cartDims[1];

            return fmt::format("({},{},{})", i + 1,j + 1,k + 1);
        };

        auto make_warning = [&format_ijk] (const Opm::KeywordLocation& location, const Opm::NNCdata& nnc) -> std::string {
            return fmt::format("Problem with EDITNNC keyword\n"
                               "In {} line {} \n"
                               "No NNC defined for connection {} -> {}", location.filename, location.lineno, format_ijk(nnc.cell1), format_ijk(nnc.cell2));

        };

        // editNnc is supposed to only reference non-neighboring connections and not
        // neighboring connections. Use all entries for scaling if there is an NNC.
        // variable nnc incremented in loop body.
        auto nnc = editNnc.begin();
        auto end = editNnc.end();
        std::size_t warning_count = 0;
        while (nnc != end) {
            auto c1 = nnc->cell1;
            auto c2 = nnc->cell2;
            auto low = globalToLocal[c1];
            auto high = globalToLocal[c2];
            if (low > high)
                std::swap(low, high);

            auto candidate = trans_.find(isId_(low, high));
            if (candidate == trans_.end()) {
                const auto& location = nnc_input.edit_location( *nnc );
                auto warning = make_warning(location, *nnc);
                Opm::OpmLog::warning("EDITNNC", warning);
                ++nnc;
                warning_count++;
            }
            else {
                // NNC exists
                while (nnc!= end && c1==nnc->cell1 && c2==nnc->cell2) {
                    candidate->second *= nnc->trans;
                    ++nnc;
                }
            }
        }

        if (warning_count > 0) {
            auto warning = fmt::format("Problems with EDITNNC keyword\n"
                                       "A total of {} connections not defined in grid", warning_count);
            Opm::OpmLog::warning(warning);
        }
    }

    void extractPermeability_()
    {
        unsigned numElem = vanguard_.gridView().size(/*codim=*/0);
        permeability_.resize(numElem);

        // read the intrinsic permeabilities from the eclState. Note that all arrays
        // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
        // simulation grid might remove a few elements. (e.g. because it is distributed
        // over several processes.)
        const auto& fp = vanguard_.eclState().fieldProps();
        if (fp.has_double("PERMX")) {
            const std::vector<double>& permxData = fp.get_double("PERMX");

            std::vector<double> permyData;
            if (fp.has_double("PERMY"))
                permyData = fp.get_double("PERMY");
            else
                permyData = permxData;

            std::vector<double> permzData;
            if (fp.has_double("PERMZ"))
                permzData = fp.get_double("PERMZ");
            else
                permzData = permxData;

            for (size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
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

    void extractPorosity_()
    {
        // read the intrinsic porosity from the eclState. Note that all arrays
        // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
        // simulation grid might remove a few elements. (e.g. because it is distributed
        // over several processes.)
        const auto& fp = vanguard_.eclState().fieldProps();
        if (fp.has_double("PORO")) {
            porosity_ = fp.get_double("PORO");
        }
        else
            throw std::logic_error("Can't read the porosityfrom the ecl state. "
                                   "(The PORO keywords are missing)");
    }

    std::uint64_t isId_(std::uint32_t elemIdx1, std::uint32_t elemIdx2) const
    {
        std::uint32_t elemAIdx = std::min(elemIdx1, elemIdx2);
        std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx<<elemIdxShift) + elemAIdx;
    }

    std::pair<std::uint32_t, std::uint32_t> isIdReverse_(const std::uint64_t& id) const
    {
        // Assigning an unsigned integer to a narrower type discards the most significant bits.
        // See "The C programming language", section A.6.2.
        // NOTE that the ordering of element A and B may have changed
        std::uint32_t elemAIdx = id;
        std::uint32_t elemBIdx = (id - elemAIdx) >> elemIdxShift;

        return std::make_pair(elemAIdx, elemBIdx);
    }

    std::uint64_t directionalIsId_(std::uint32_t elemIdx1, std::uint32_t elemIdx2) const
    {
        return (std::uint64_t(elemIdx1)<<elemIdxShift) + elemIdx2;
    }

    void computeHalfTrans_(Scalar& halfTrans,
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

    void computeHalfDiffusivity_(Scalar& halfDiff,
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

    DimVector distanceVector_(const DimVector& center,
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

    void applyMultipliers_(Scalar& trans,
                           unsigned faceIdx,
                           unsigned cartElemIdx,
                           const Opm::TransMult& transMult) const
    {
        // apply multiplyer for the transmissibility of the face. (the
        // face index is the index of the reference-element face which
        // contains the intersection of interest.)
        switch (faceIdx) {
        case 0: // left
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::XMinus);
            break;
        case 1: // right
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::XPlus);
            break;

        case 2: // front
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::YMinus);
            break;
        case 3: // back
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::YPlus);
            break;

        case 4: // bottom
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZMinus);
            break;
        case 5: // top
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZPlus);
            break;
        }
    }

    void applyNtg_(Scalar& trans,
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

    const Vanguard& vanguard_;
    Scalar transmissibilityThreshold_;
    std::vector<DimMatrix> permeability_;
    std::vector<Scalar> porosity_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
    std::map<std::pair<unsigned, unsigned>, Scalar> transBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> thermalHalfTransBoundary_;
    Opm::ConditionalStorage<enableEnergy,
                            std::unordered_map<std::uint64_t, Scalar> > thermalHalfTrans_;
    Opm::ConditionalStorage<enableDiffusion,
                            std::unordered_map<std::uint64_t, Scalar> > diffusivity_;
};

} // namespace Opm

#endif
