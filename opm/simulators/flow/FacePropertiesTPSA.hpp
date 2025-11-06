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
#ifndef TPSA_FACE_PROPERTIES_HPP
#define TPSA_FACE_PROPERTIES_HPP

#include <dune/common/fvector.hh>

#include <opm/grid/LookUpData.hh>

#include <cstdint>
#include <map>
#include <vector>
#include <utility>


namespace Opm {

class EclipseState;


template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
class FacePropertiesTPSA {

    enum { dimWorld = GridView::dimensionworld };

public:
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

    FacePropertiesTPSA(const EclipseState& eclState,
                       const GridView& gridView,
                       const CartesianIndexMapper& cartMapper,
                       const Grid& grid,
                       std::function<std::array<double, dimWorld>(int)> centroids);

    void finishInit();

    void update();

    Scalar weightAverage(unsigned elemIdx1, unsigned elemIdx2) const;
    Scalar weightAverageBoundary(unsigned elemIdx1, unsigned boundaryFaceIdx) const;
    Scalar weightProduct(unsigned elemIdx1, unsigned elemIdx2) const;
    Scalar weightProductBoundary(unsigned elemIdx1, unsigned boundaryFaceIdx) const;
    Scalar normalDistance(unsigned elemIdx1, unsigned elemIdx2) const;
    Scalar normalDistanceBoundary(unsigned elemIdx1, unsigned boundaryFaceIdx) const;
    DimVector cellFaceNormal(unsigned elemIdx1, unsigned elemIdx2);
    const DimVector& cellFaceNormalBoundary(unsigned elemIdx1, unsigned boundaryFaceIdx) const;

    /*!
    * \brief Return shear modulus of an element
    */
    const Scalar shearModulus(unsigned elemIdx) const
    { return sModulus_[elemIdx]; }


protected:
    struct FaceInfo
    {
        DimVector faceCenter;
        int faceIdx;
        unsigned elemIdx;
        unsigned cartElemIdx;
    };

    template <class Intersection>
    void computeCellProperties(const Intersection& intersection,
                               FaceInfo& inside,
                               FaceInfo& outside,
                               DimVector& faceNormal,
                               /*isCpGrid=*/std::false_type) const;

    template <class Intersection>
    void computeCellProperties(const Intersection& intersection,
                               FaceInfo& inside,
                               FaceInfo& outside,
                               DimVector& faceNormal,
                               /*isCpGrid=*/std::true_type) const;

    Scalar computeDistance_(const DimVector& distVec, const DimVector& faceNormal);

    Scalar computeWeight_(const Scalar distance, const Scalar smod);

    DimVector distanceVector_(const DimVector& faceCenter, const unsigned& cellIdx) const;

    void extractSModulus_();

    std::vector<Scalar> sModulus_;
    std::unordered_map<std::uint64_t, Scalar> weightsAvg_;
    std::unordered_map<std::uint64_t, Scalar> weightsProd_;
    std::unordered_map<std::uint64_t, Scalar> distance_;
    std::unordered_map<std::uint64_t, DimVector> faceNormal_;

    std::map<std::pair<unsigned, unsigned>, Scalar> weightsAvgBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> weightsProdBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> distanceBoundary_;
    std::map<std::pair<unsigned, unsigned>, DimVector> faceNormalBoundary_;

    const EclipseState& eclState_;
    const GridView& gridView_;
    const CartesianIndexMapper& cartMapper_;
    const Grid& grid_;
    std::function<std::array<double, dimWorld>(int)> centroids_;
    std::vector<std::array<double, dimWorld>> centroids_cache_;

    const LookUpData<Grid, GridView> lookUpData_;
    const LookUpCartesianData<Grid, GridView> lookUpCartesianData_;

};  // FacePropertiesTPSA

namespace details {
    std::uint64_t isIdTPSA(std::uint32_t elemIdx1, std::uint32_t elemIdx2);
}  // namespace details

}  // namespace Opm

#endif