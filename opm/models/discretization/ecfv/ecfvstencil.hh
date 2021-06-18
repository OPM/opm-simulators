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
 * \copydoc Opm::EcfvStencil
 */
#ifndef EWOMS_ECFV_STENCIL_HH
#define EWOMS_ECFV_STENCIL_HH

#include <opm/models/utils/quadraturegeometries.hh>

#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/geometry/type.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Opm {
/*!
 * \ingroup EcfvDiscretization
 *
 * \brief Represents the stencil (finite volume geometry) of a single
 *        element in the ECFV discretization.
 *
 * The ECFV discretization is a element centered finite volume
 * approach. This means that each element corresponds to a control
 * volume.
 */
template <class Scalar,
          class GridView,
          bool needFaceIntegrationPos = true,
          bool needFaceNormal = true>
class EcfvStencil
{
    enum { dimWorld = GridView::dimensionworld };

    using CoordScalar = typename GridView::ctype;
    using Intersection = typename GridView::Intersection;
    using Element = typename GridView::template Codim<0>::Entity;

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using WorldVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    using Entity = Element       ;
    using Mapper = ElementMapper ;

    using LocalGeometry = typename Element::Geometry;

    /*!
     * \brief Represents a sub-control volume.
     *
     * For element centered finite volumes, this is equivalent to the
     * element, in the vertex centered finite volume approach, this
     * corresponds to the intersection of a finite volume and the
     * grid element.
     */
    class SubControlVolume
    {
    public:
        // default construct an uninitialized object.
        // this is only here because std::vector needs it...
        SubControlVolume()
        {}

        SubControlVolume(const Element& element)
            : element_(element)
        { update(); }

        void update(const Element& element)
        { element_ = element; }

        void update()
        {
            const auto& geometry = element_.geometry();
            centerPos_ = geometry.center();
            volume_ = geometry.volume();
        }

        /*!
         * \brief The global position associated with the sub-control volume
         */
        const GlobalPosition& globalPos() const
        { return centerPos_; }

        /*!
         * \brief The center of the sub-control volume
         */
        const GlobalPosition& center() const
        { return centerPos_; }

        /*!
         * \brief The volume [m^3] occupied by the sub-control volume
         */
        Scalar volume() const
        { return volume_; }

        /*!
         * \brief The geometry of the sub-control volume.
         */
        const LocalGeometry geometry() const
        { return element_.geometry(); }

        /*!
         * \brief Geometry of the sub-control volume relative to parent.
         */
        const LocalGeometry localGeometry() const
        { return element_.geometryInFather(); }

    private:
        GlobalPosition centerPos_;
        Scalar volume_;
        Element element_;
    };

    /*!
     * \brief Represents a face of a sub-control volume.
     */
    template <bool needIntegrationPos, bool needNormal>
    class EcfvSubControlVolumeFace
    {
    public:
        EcfvSubControlVolumeFace()
        {}

        EcfvSubControlVolumeFace(const Intersection& intersection, unsigned localNeighborIdx)
        {
            exteriorIdx_ = static_cast<unsigned short>(localNeighborIdx);

            if (needNormal)
                (*normal_) = intersection.centerUnitOuterNormal();

            const auto& geometry = intersection.geometry();
            if (needIntegrationPos)
                (*integrationPos_) = geometry.center();
            area_ = geometry.volume();
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's interior.
         */
        unsigned short interiorIndex() const
        {
            // The local index of the control volume in the interior
            // of a face of the stencil in the element centered finite
            // volume discretization is always the "central"
            // element. In this implementation this element always has
            // index 0....
            return 0;
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's outside.
         */
        unsigned short exteriorIndex() const
        { return exteriorIdx_; }

        /*!
         * \brief Returns the global position of the face's
         *        integration point.
         */
        const GlobalPosition& integrationPos() const
        { return *integrationPos_; }

        /*!
         * \brief Returns the outer unit normal at the face's
         *        integration point.
         */
        const WorldVector& normal() const
        { return *normal_; }

        /*!
         * \brief Returns the area [m^2] of the face
         */
        Scalar area() const
        { return area_; }

    private:
        ConditionalStorage<needIntegrationPos, GlobalPosition> integrationPos_;
        ConditionalStorage<needNormal, WorldVector> normal_;
        Scalar area_;

        unsigned short exteriorIdx_;
    };

    using SubControlVolumeFace = EcfvSubControlVolumeFace<needFaceIntegrationPos, needFaceNormal>;
    using BoundaryFace = EcfvSubControlVolumeFace</*needFaceIntegrationPos=*/true, needFaceNormal>;

    EcfvStencil(const GridView& gridView, const Mapper& mapper)
        : gridView_(gridView)
        , elementMapper_(mapper)
    {
        // try to ensure that the mapper passed indeed maps elements
        assert(int(gridView.size(/*codim=*/0)) == int(elementMapper_.size()));
    }

    void updateTopology(const Element& element)
    {
        auto isIt = gridView_.ibegin(element);
        const auto& endIsIt = gridView_.iend(element);

        // add the "center" element of the stencil
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/element/*)*/);
        elements_.clear();
        elements_.emplace_back(element);

        interiorFaces_.clear();
        boundaryFaces_.clear();

        for (; isIt != endIsIt; ++isIt) {
            const auto& intersection = *isIt;
            // if the current intersection has a neighbor, add a
            // degree of freedom and an internal face, else add a
            // boundary face
            if (intersection.neighbor()) {
                elements_.emplace_back( intersection.outside() );
                subControlVolumes_.emplace_back(/*SubControlVolume(*/elements_.back()/*)*/);
                interiorFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, subControlVolumes_.size() - 1/*)*/);
            }
            else {
                boundaryFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, - 10000/*)*/);
            }
        }
    }

    void updatePrimaryTopology(const Element& element)
    {
        // add the "center" element of the stencil
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/element/*)*/);
        elements_.clear();
        elements_.emplace_back(element);
    }

    void update(const Element& element)
    {
        updateTopology(element);
    }

    void updateCenterGradients()
    {
        assert(false); // not yet implemented
    }

    /*!
     * \brief Return the element to which the stencil refers.
     */
    const Element& element() const
    { return element( 0 ); }

    /*!
     * \brief Return the grid view of the element to which the stencil
     *        refers.
     */
    const GridView& gridView() const
    { return *gridView_; }

    /*!
     * \brief Returns the number of degrees of freedom which the
     *        current element interacts with.
     */
    size_t numDof() const
    { return subControlVolumes_.size(); }

    /*!
     * \brief Returns the number of degrees of freedom which are contained
     *        by within the current element.
     *
     * Primary DOFs are always expected to have a lower index than
     * "secondary" DOFs.
     *
     * For element centered finite elements, this is only the central DOF.
     */
    size_t numPrimaryDof() const
    { return 1; }

    /*!
     * \brief Return the global space index given the index of a degree of
     *        freedom.
     */
    unsigned globalSpaceIndex(unsigned dofIdx) const
    {
        assert(dofIdx < numDof());

        return static_cast<unsigned>(elementMapper_.index(element(dofIdx)));
    }

    /*!
     * \brief Return partition type of a given degree of freedom
     */
    Dune::PartitionType partitionType(unsigned dofIdx) const
    { return elements_[dofIdx]->partitionType(); }

    /*!
     * \brief Return the element given the index of a degree of
     *        freedom.
     *
     * If no degree of freedom index is passed, the element which was
     * passed to the update() method is returned...
     */
    const Element& element(unsigned dofIdx) const
    {
        assert(dofIdx < numDof());

        return elements_[dofIdx];
    }

    /*!
     * \brief Return the entity given the index of a degree of
     *        freedom.
     */
    const Entity& entity(unsigned dofIdx) const
    {
        return element( dofIdx );
    }

    /*!
     * \brief Returns the sub-control volume object belonging to a
     *        given degree of freedom.
     */
    const SubControlVolume& subControlVolume(unsigned dofIdx) const
    { return subControlVolumes_[dofIdx]; }

    /*!
     * \brief Returns the number of interior faces of the stencil.
     */
    size_t numInteriorFaces() const
    { return interiorFaces_.size(); }

    /*!
     * \brief Returns the face object belonging to a given face index
     *        in the interior of the domain.
     */
    const SubControlVolumeFace& interiorFace(unsigned faceIdx) const
    { return interiorFaces_[faceIdx]; }

    /*!
     * \brief Returns the number of boundary faces of the stencil.
     */
    size_t numBoundaryFaces() const
    { return boundaryFaces_.size(); }

    /*!
     * \brief Returns the boundary face object belonging to a given
     *        boundary face index.
     */
    const BoundaryFace& boundaryFace(unsigned bfIdx) const
    { return boundaryFaces_[bfIdx]; }

protected:
    const GridView&       gridView_;
    const ElementMapper&  elementMapper_;

    std::vector<Element> elements_;
    std::vector<SubControlVolume>      subControlVolumes_;
    std::vector<SubControlVolumeFace>  interiorFaces_;
    std::vector<BoundaryFace>  boundaryFaces_;
};

} // namespace Opm


#endif

