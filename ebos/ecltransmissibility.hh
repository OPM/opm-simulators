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
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/CpGrid.hpp>

#include <array>
#include <vector>
#include <unordered_map>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridManager);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementMapper);
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
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) GridManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GridView::Intersection Intersection;

    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    EclTransmissibility(const GridManager& gridManager)
        : gridManager_(gridManager)
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


    void update()
    {
        const auto& gridView = gridManager_.gridView();
        const auto& cartMapper = gridManager_.cartesianIndexMapper();
        const auto& eclState = gridManager_.eclState();
        const auto& eclGrid = eclState.getInputGrid();
        auto& transMult = eclState.getTransMult();
        ElementMapper elemMapper(gridView);

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

        // compute the transmissibilities for all intersections
        elemIt = gridView.template begin</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;

                // ignore boundary intersections for now (TODO?)

                // continue if no neighbor is present
                if ( ! intersection.neighbor() )
                    continue;

                const auto& inside = intersection.inside();
                const auto& outside = intersection.outside();
                unsigned insideElemIdx = elemMapper.index(inside);
                unsigned outsideElemIdx = elemMapper.index(outside);

                // we only need to calculate a face's transmissibility
                // once...
                if (insideElemIdx > outsideElemIdx)
                    continue;

                unsigned insideCartElemIdx = cartMapper.cartesianIndex(insideElemIdx);
                unsigned outsideCartElemIdx = cartMapper.cartesianIndex(outsideElemIdx);

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                unsigned insideFaceIdx  = intersection.indexInInside();
                unsigned outsideFaceIdx = intersection.indexInOutside();

                DimVector faceCenterInside;
                DimVector faceCenterOutside;
                DimVector faceAreaNormal;

                typename std::is_same< Grid, Dune::CpGrid> :: type isCpGrid;
                computeFaceProperties( intersection, insideElemIdx, insideFaceIdx, outsideElemIdx, outsideFaceIdx,
                                       faceCenterInside, faceCenterOutside, faceAreaNormal,
                                       isCpGrid );

                Scalar halfTrans1;
                Scalar halfTrans2;

                computeHalfTrans_(halfTrans1,
                                  faceAreaNormal,
                                  insideFaceIdx,
                                  distanceVector_(faceCenterInside,
                                                  intersection.indexInInside(),
                                                  insideElemIdx,
                                                  axisCentroids),
                                  permeability_[insideElemIdx]);
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
                    OPM_THROW(std::logic_error, "Could not determine a face direction");
                }

                trans *= transMult.getRegionMultiplier(insideCartElemIdx,
                                                       outsideCartElemIdx,
                                                       faceDir);

                trans_[isId_(insideElemIdx, outsideElemIdx)] = trans;
            }
        }
    }

    const DimMatrix& permeability(unsigned elemIdx) const
    { return permeability_[elemIdx]; }

    Scalar transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
    { return trans_.at(isId_(elemIdx1, elemIdx2)); }

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
                                std::false_type ) const
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
                                std::true_type ) const
    {
        int faceIdx = intersection.id();
        faceCenterInside = gridManager_.grid().faceCenterEcl(insideElemIdx,insideFaceIdx);
        faceCenterOutside = gridManager_.grid().faceCenterEcl(outsideElemIdx,outsideFaceIdx);
        faceAreaNormal = gridManager_.grid().faceAreaNormalEcl(faceIdx);
    }

    void extractPermeability_()
    {
        const auto& props = gridManager_.eclState().get3DProperties();

        unsigned numElem = gridManager_.gridView().size(/*codim=*/0);
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
                unsigned cartesianElemIdx = gridManager_.cartesianIndex(dofIdx);
                permeability_[dofIdx] = 0.0;
                permeability_[dofIdx][0][0] = permxData[cartesianElemIdx];
                permeability_[dofIdx][1][1] = permyData[cartesianElemIdx];
                permeability_[dofIdx][2][2] = permzData[cartesianElemIdx];
            }

            // for now we don't care about non-diagonal entries
        }
        else
            OPM_THROW(std::logic_error,
                      "Can't read the intrinsic permeability from the ecl state. "
                      "(The PERM{X,Y,Z} keywords are missing)");
    }

    std::uint64_t isId_(unsigned elemIdx1, unsigned elemIdx2) const
    {
        static const unsigned elemIdxShift = 32; // bits

        unsigned elemAIdx = std::min(elemIdx1, elemIdx2);
        std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx<<elemIdxShift) + elemAIdx;
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
        const auto& eclState = gridManager_.eclState();
        const auto& eclGrid = eclState.getInputGrid();
        const auto& porv = eclState.get3DProperties().getDoubleGridProperty("PORV").getData();
        const auto& actnum = eclState.get3DProperties().getIntGridProperty("ACTNUM").getData();
        const std::vector<double>& ntg =
            eclState.get3DProperties().getDoubleGridProperty("NTG").getData();

        const auto& cartMapper = gridManager_.cartesianIndexMapper();
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




    const GridManager& gridManager_;
    std::vector<DimMatrix> permeability_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
};

} // namespace Ewoms

#endif
