/*
  Copyright (C) 2014 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::EclTransmissibility
 */
#ifndef EWOMS_ECL_TRANSMISSIBILITY_HH
#define EWOMS_ECL_TRANSMISSIBILITY_HH

#include "eclgridmanager.hh"

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <array>
#include <vector>
#include <unordered_map>

namespace Ewoms {
/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the transmissibilites for grid faces according to the
 *        Eclipse Technical Description.
 */
template <class TypeTag>
class EclTransmissibility
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GridView::Intersection Intersection;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    EclTransmissibility(const Simulator& simulator)
        : simulator_(simulator)
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
    {
        const auto& elementMapper = simulator_.model().elementMapper();
        const auto& gridView = simulator_.gridView();
        const auto& problem = simulator_.problem();

        int numElements = elementMapper.size();

        // this code assumes that the DOFs are the elements. (i.e., an
        // ECFV spatial discretization with TPFA). if you try to use
        // it with something else, you're currently out of luck,
        // sorry!
        assert(simulator_.model().numGridDof() == numElements);

        // calculate the axis specific centroids of all elements
        std::array<std::vector<DimVector>, dimWorld> axisCentroids;

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            axisCentroids[dimIdx].resize(numElements);

        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
            int elemIdx = elementMapper.index(elemIt);
#else
            int elemIdx = elementMapper.map(elemIt);
#endif

            // get the geometry of the current element
            const auto& geom = elemIt->geometry();

            // compute the axis specific "centroids" used for the
            // transmissibilities
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                DimVector x0Local(0.5);
                DimVector x1Local(0.5);

                x0Local[dimIdx] = 0.0;
                x1Local[dimIdx] = 1.0;

                DimVector x = geom.global(x0Local);
                x += geom.global(x1Local);
                x /= 2;

                axisCentroids[dimIdx][elemIdx] = x;
            }
        }

        Opm::EclipseStateConstPtr eclState = simulator_.gridManager().eclState();
        const std::vector<double>& multx =
            eclState->getDoubleGridProperty("MULTX")->getData();
        const std::vector<double>& multy =
            eclState->getDoubleGridProperty("MULTY")->getData();
        const std::vector<double>& multz =
            eclState->getDoubleGridProperty("MULTZ")->getData();
        const std::vector<double>& multxMinus =
            eclState->getDoubleGridProperty("MULTX-")->getData();
        const std::vector<double>& multyMinus =
            eclState->getDoubleGridProperty("MULTY-")->getData();
        const std::vector<double>& multzMinus =
            eclState->getDoubleGridProperty("MULTZ-")->getData();

        const std::vector<double>& ntg =
            eclState->getDoubleGridProperty("NTG")->getData();

        // reserving some space in the hashmap upfront saves quite a
        // bit of time because resizes are costly for hashmaps and
        // there would be quite a few of them if we would not have a
        // rough idea of how large the final map will be (the rough
        // idea is a conforming Cartesian grid)
        trans_.reserve(numElements*3*1.05);

        // compute the transmissibilities for all intersections
        elemIt = gridView.template begin</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            auto isIt = elemIt->ileafbegin();
            const auto& isEndIt = elemIt->ileafend();
            for (; isIt != isEndIt; ++ isIt) {
                // ignore boundary intersections for now (TODO?)
                if (isIt->boundary())
                    continue;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
                int insideElemIdx = elementMapper.index(*isIt->inside());
                int outsideElemIdx = elementMapper.index(*isIt->outside());
#else
                int insideElemIdx = elementMapper.map(*isIt->inside());
                int outsideElemIdx = elementMapper.map(*isIt->outside());
#endif

                // we only need to calculate a face's transmissibility
                // once...
                if (insideElemIdx > outsideElemIdx)
                    continue;

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                int insideFaceIdx = isIt->indexInInside();
                int outsideFaceIdx = isIt->indexInOutside();

                Scalar halfTrans1;
                Scalar halfTrans2;

                computeHalfTrans_(halfTrans1,
                                  *isIt,
                                  insideFaceIdx,
                                  distanceVector_(*isIt,
                                                  isIt->indexInInside(),
                                                  insideElemIdx,
                                                  axisCentroids),
                                  problem.intrinsicPermeability(insideElemIdx));
                computeHalfTrans_(halfTrans2,
                                  *isIt,
                                  outsideFaceIdx,
                                  distanceVector_(*isIt,
                                                  isIt->indexInOutside(),
                                                  outsideElemIdx,
                                                  axisCentroids),
                                  problem.intrinsicPermeability(outsideElemIdx));

                applyNtg_(halfTrans1, insideFaceIdx, insideElemIdx, ntg);
                applyNtg_(halfTrans2, outsideFaceIdx, outsideElemIdx, ntg);

                // convert half transmissibilities to full face
                // transmissibilities using the harmonic mean
                Scalar trans = 1.0 / (1.0/halfTrans1 + 1.0/halfTrans2);

                // apply the full face transmissibility multipliers
                // for the inside ...
                applyMultipliers_(trans, insideFaceIdx, insideElemIdx,
                                  multx, multxMinus,
                                  multy, multyMinus,
                                  multz, multzMinus);
                // ... and outside elements
                applyMultipliers_(trans, outsideFaceIdx, outsideElemIdx,
                                  multx, multxMinus,
                                  multy, multyMinus,
                                  multz, multzMinus);

                trans_[isId_(insideElemIdx, outsideElemIdx)] = trans;
            }
        }
    }

    Scalar transmissibility(int elemIdx1, int elemIdx2) const
    { return trans_.at(isId_(elemIdx1, elemIdx2)); }

private:
    std::uint64_t isId_(int elemIdx1, int elemIdx2) const
    {
        static const int elemIdxShift = 32; // bits

        int elemAIdx = std::min(elemIdx1, elemIdx2);
        std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx<<elemIdxShift) + elemAIdx;
    }

    void computeHalfTrans_(Scalar& halfTrans,
                           const Intersection& is,
                           int faceIdx, // in the reference element that contains the intersection
                           const DimVector& distance,
                           const DimMatrix& perm) const
    {
        int dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        halfTrans = perm[dimIdx][dimIdx];
        halfTrans *= is.geometry().volume();
        halfTrans *= std::abs<Scalar>(is.centerUnitOuterNormal()*distance);
        halfTrans /= distance*distance;
    }

    DimVector distanceVector_(const Intersection& is,
                              int faceIdx, // in the reference element that contains the intersection
                              int elemIdx,
                              const std::array<std::vector<DimVector>, dimWorld>& axisCentroids) const
    {
        int dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        DimVector x = is.geometry().center();
        x -= axisCentroids[dimIdx][elemIdx];

        return x;
    }

    void applyMultipliers_(Scalar &trans, int faceIdx, int elemIdx,
                           const std::vector<Scalar>& multx,
                           const std::vector<Scalar>& multxMinus,
                           const std::vector<Scalar>& multy,
                           const std::vector<Scalar>& multyMinus,
                           const std::vector<Scalar>& multz,
                           const std::vector<Scalar>& multzMinus) const
    {
        // apply multiplyer for the transmissibility of the face. (the
        // face index is the index of the reference-element face which
        // contains the intersection of interest.)
        switch (faceIdx) {
        case 0: // left
            trans *= multxMinus[elemIdx];
            break;
        case 1: // right
            trans *= multx[elemIdx];
            break;

        case 2: // front
            trans *= multyMinus[elemIdx];
            break;
        case 3: // back
            trans *= multy[elemIdx];
            break;

        case 4: // bottom
            trans *= multzMinus[elemIdx];
            break;
        case 5: // top
            trans *= multz[elemIdx];
            break;
        }
    }

    void applyNtg_(Scalar &trans, int faceIdx, int elemIdx,
                   const std::vector<Scalar>& ntg) const
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

    const Simulator& simulator_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
};

} // namespace Ewoms

#endif
