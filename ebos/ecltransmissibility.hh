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
 * \copydoc Ewoms::EclTransmissibility
 */
#ifndef EWOMS_ECL_TRANSMISSIBILITY_HH
#define EWOMS_ECL_TRANSMISSIBILITY_HH

#include <ewoms/common/propertysystem.hh>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/TransMult.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <array>
#include <vector>
#include <unordered_map>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Vanguard);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementMapper);
NEW_PROP_TAG(EnableEnergy);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the transmissibilites for grid faces according to the
 *        Eclipse Technical Description.
 */
template <class TypeTag>
class EclTransmissibility
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GridView::Intersection Intersection;

    static const bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);

    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    EclTransmissibility(const Vanguard& vanguard)
        : vanguard_(vanguard)
    {}

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
    { update(); }


    /*!
     * \brief Compute all transmissibilities
     *
     * Also, this updates the "thermal half transmissibilities" if energy is enabled.
     */
    void update()
    {
        const auto& gridView = vanguard_.gridView();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& eclState = vanguard_.eclState();
        const auto& eclGrid = eclState.getInputGrid();
        auto& transMult = eclState.getTransMult();
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
#else
        ElementMapper elemMapper(gridView);
#endif

        // get the ntg values, the ntg values are modified for the cells merged with minpv
        std::vector<double> ntg;
        minPvFillNtg_(ntg);

        unsigned numElements = elemMapper.size();

        extractPermeability_();

        // calculate the axis specific centroids of all elements
        std::array<std::vector<DimVector>, dimWorld> axisCentroids;

        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            axisCentroids[dimIdx].resize(numElements);

        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            unsigned elemIdx = elemMapper.index(elem);

            // compute the axis specific "centroids" used for the transmissibilities. for
            // consistency with the flow simulator, we use the element centers as
            // computed by opm-parser's Opm::EclipseGrid class for all axes.
            unsigned cartesianCellIdx = cartMapper.cartesianIndex(elemIdx);
            const auto& centroid = eclGrid.getCellCenter(cartesianCellIdx);
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

                if (!intersection.neighbor())
                    // elements can be on process boundaries, i.e. they are not on the
                    // domain boundary yet they don't have neighbors.
                    continue;

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

                // we only need to calculate a face's transmissibility
                // once...
                if (elemIdx > outsideElemIdx)
                    continue;

                unsigned insideCartElemIdx = cartMapper.cartesianIndex(elemIdx);
                unsigned outsideCartElemIdx = cartMapper.cartesianIndex(outsideElemIdx);

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                unsigned insideFaceIdx  = intersection.indexInInside();
                unsigned outsideFaceIdx = intersection.indexInOutside();

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

                applyNtg_(halfTrans1, insideFaceIdx, insideCartElemIdx, ntg);
                applyNtg_(halfTrans2, outsideFaceIdx, outsideCartElemIdx, ntg);

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
                applyMultipliers_(trans, insideFaceIdx, insideCartElemIdx, transMult);
                // ... and outside elements
                applyMultipliers_(trans, outsideFaceIdx, outsideCartElemIdx, transMult);

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
            }
        }
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

private:
    template <class Intersection>
    void computeFaceProperties( const Intersection& intersection,
                                const int insideElemIdx,
                                const int insideFaceIdx,
                                const int outsideElemIdx,
                                const int outsideFaceIdx,
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
    void computeFaceProperties( const Intersection& intersection,
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
        faceCenterInside = vanguard_.grid().faceCenterEcl(insideElemIdx,insideFaceIdx);
        faceCenterOutside = vanguard_.grid().faceCenterEcl(outsideElemIdx,outsideFaceIdx);
        faceAreaNormal = vanguard_.grid().faceAreaNormalEcl(faceIdx);
    }

    void extractPermeability_()
    {
        const auto& props = vanguard_.eclState().get3DProperties();

        unsigned numElem = vanguard_.gridView().size(/*codim=*/0);
        permeability_.resize(numElem);

        // read the intrinsic permeabilities from the eclState. Note that all arrays
        // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
        // simulation grid might remove a few elements. (e.g. because it is distributed
        // over several processes.)
        if (props.hasDeckDoubleGridProperty("PERMX")) {
            const std::vector<double>& permxData =
                props.getDoubleGridProperty("PERMX").getData();
            std::vector<double> permyData(permxData);
            if (props.hasDeckDoubleGridProperty("PERMY"))
                permyData = props.getDoubleGridProperty("PERMY").getData();
            std::vector<double> permzData(permxData);
            if (props.hasDeckDoubleGridProperty("PERMZ"))
                permzData = props.getDoubleGridProperty("PERMZ").getData();

            for (size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
                unsigned cartesianElemIdx = vanguard_.cartesianIndex(dofIdx);
                permeability_[dofIdx] = 0.0;
                permeability_[dofIdx][0][0] = permxData[cartesianElemIdx];
                permeability_[dofIdx][1][1] = permyData[cartesianElemIdx];
                permeability_[dofIdx][2][2] = permzData[cartesianElemIdx];
            }

            // for now we don't care about non-diagonal entries
        }
        else
            throw std::logic_error("Can't read the intrinsic permeability from the ecl state. "
                                   "(The PERM{X,Y,Z} keywords are missing)");
    }

    std::uint64_t isId_(unsigned elemIdx1, unsigned elemIdx2) const
    {
        static const unsigned elemIdxShift = 32; // bits

        unsigned elemAIdx = std::min(elemIdx1, elemIdx2);
        std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx<<elemIdxShift) + elemAIdx;
    }

    std::uint64_t directionalIsId_(unsigned elemIdx1, unsigned elemIdx2) const
    {
        static const unsigned elemIdxShift = 32; // bits

        return (std::uint64_t(elemIdx1)<<elemIdxShift) + elemIdx2;
    }

    void computeHalfTrans_(Scalar& halfTrans,
                           const DimVector& areaNormal,
                           unsigned faceIdx, // in the reference element that contains the intersection
                           const DimVector& distance,
                           const DimMatrix& perm) const
    {
        unsigned dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        halfTrans = perm[dimIdx][dimIdx];

        Scalar val = 0;
        for (unsigned i = 0; i < areaNormal.size(); ++i)
            val += areaNormal[i]*distance[i];

        halfTrans *= std::abs(val);
        halfTrans /= distance.two_norm2();
    }

    DimVector distanceVector_(const DimVector& center,
                              unsigned faceIdx, // in the reference element that contains the intersection
                              unsigned elemIdx,
                              const std::array<std::vector<DimVector>, dimWorld>& axisCentroids) const
    {
        unsigned dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        DimVector x = center;
        x -= axisCentroids[dimIdx][elemIdx];

        return x;
    }

    void applyMultipliers_(Scalar& trans, unsigned faceIdx, unsigned cartElemIdx,
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

    void applyNtg_(Scalar& trans, unsigned faceIdx, unsigned cartElemIdx,
                   const std::vector<double>& ntg) const
    {
        // apply multiplyer for the transmissibility of the face. (the
        // face index is the index of the reference-element face which
        // contains the intersection of interest.)
        switch (faceIdx) {
        case 0: // left
            trans *= ntg[cartElemIdx];
            break;
        case 1: // right
            trans *= ntg[cartElemIdx];
            break;

        case 2: // front
            trans *= ntg[cartElemIdx];
            break;
        case 3: // back
            trans *= ntg[cartElemIdx];
            break;

            // NTG does not apply to top and bottom faces
        }
    }

    void minPvFillNtg_(std::vector<double>& averageNtg) const {

        // compute volume weighted arithmetic average of NTG for
        // cells merged as an result of minpv.
        const auto& eclState = vanguard_.eclState();
        const auto& eclGrid = eclState.getInputGrid();
        const auto& porv = eclState.get3DProperties().getDoubleGridProperty("PORV").getData();
        const auto& actnum = eclState.get3DProperties().getIntGridProperty("ACTNUM").getData();
        const std::vector<double>& ntg =
            eclState.get3DProperties().getDoubleGridProperty("NTG").getData();

        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        assert(dimWorld > 1);
        const size_t nxny = cartDims[0] * cartDims[1];

        averageNtg = ntg;

        for (size_t cartesianCellIdx = 0; cartesianCellIdx < ntg.size(); ++cartesianCellIdx)
        {
            // use the original ntg values for the inactive cells
            if (!actnum[cartesianCellIdx])
                continue;

            // Average properties as long as there exist cells above
            // that has pore volume less than the MINPV threshold
            const double cellVolume = eclGrid.getCellVolume(cartesianCellIdx);
            double ntgCellVolume = ntg[cartesianCellIdx] * cellVolume;
            double totalCellVolume = cellVolume;
            int cartesianCellIdxAbove = cartesianCellIdx - nxny;
            while ( cartesianCellIdxAbove >= 0 &&
                 actnum[cartesianCellIdxAbove] > 0 &&
                 porv[cartesianCellIdxAbove] < eclGrid.getMinpvValue() ) {

                // Volume weighted arithmetic average of NTG
                const double cellAboveVolume = eclGrid.getCellVolume(cartesianCellIdxAbove);
                totalCellVolume += cellAboveVolume;
                ntgCellVolume += ntg[cartesianCellIdxAbove]*cellAboveVolume;
                cartesianCellIdxAbove -= nxny;
            }
            averageNtg[cartesianCellIdx] = ntgCellVolume / totalCellVolume;
        }
    }

    const Vanguard& vanguard_;
    std::vector<DimMatrix> permeability_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
    std::map<std::pair<unsigned, unsigned>, Scalar> transBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> thermalHalfTransBoundary_;
    Opm::ConditionalStorage<enableEnergy,
                            std::unordered_map<std::uint64_t, Scalar> > thermalHalfTrans_;
};

} // namespace Ewoms

#endif
