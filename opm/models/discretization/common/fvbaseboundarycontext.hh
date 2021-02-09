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
 * \copydoc Opm::FvBaseBoundaryContext
 */
#ifndef EWOMS_FV_BASE_BOUNDARY_CONTEXT_HH
#define EWOMS_FV_BASE_BOUNDARY_CONTEXT_HH

#include "fvbaseproperties.hh"

#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Represents all quantities which available on boundary segments
 */
template<class TypeTag>
class FvBaseBoundaryContext
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using GradientCalculator = GetPropType<TypeTag, Properties::GradientCalculator>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using IntersectionIterator = typename GridView::IntersectionIterator;
    using Intersection = typename GridView::Intersection;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Vector = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor.
     */
    explicit FvBaseBoundaryContext(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
        , intersectionIt_(gridView().ibegin(element()))
    { }

    void increment()
    {
        const auto& iend = gridView().iend(element());

        if(intersectionIt_ == iend)
          return;

        ++intersectionIt_;
        // iterate to the next boundary intersection
        while (intersectionIt_ != iend && intersectionIt_->neighbor()) {
            ++intersectionIt_;
        }
    }

    /*!
     * \copydoc Opm::ElementContext::problem()
     */
    const Problem& problem() const
    { return elemCtx_.problem(); }

    /*!
     * \copydoc Opm::ElementContext::model()
     */
    const Model& model() const
    { return elemCtx_.model(); }

    /*!
     * \copydoc Opm::ElementContext::gridView()
     */
    const GridView& gridView() const
    { return elemCtx_.gridView(); }

    /*!
     * \copydoc Opm::ElementContext::element()
     */
    const Element& element() const
    { return elemCtx_.element(); }

    /*!
     * \brief Returns a reference to the element context object.
     */
    const ElementContext& elementContext() const
    { return elemCtx_; }

    /*!
     * \brief Returns a reference to the current gradient calculator.
     */
    const GradientCalculator& gradientCalculator() const
    { return elemCtx_.gradientCalculator(); }

    /*!
     * \copydoc Opm::ElementContext::numDof()
     */
    size_t numDof(unsigned timeIdx) const
    { return elemCtx_.numDof(timeIdx); }

    /*!
     * \copydoc Opm::ElementContext::numPrimaryDof()
     */
    size_t numPrimaryDof(unsigned timeIdx) const
    { return elemCtx_.numPrimaryDof(timeIdx); }

    /*!
     * \copydoc Opm::ElementContext::numInteriorFaces()
     */
    size_t numInteriorFaces(unsigned timeIdx) const
    { return elemCtx_.numInteriorFaces(timeIdx); }

    /*!
     * \brief Return the number of boundary segments of the current element
     */
    size_t numBoundaryFaces(unsigned timeIdx) const
    { return elemCtx_.stencil(timeIdx).numBoundaryFaces(); }

    /*!
     * \copydoc Opm::ElementContext::stencil()
     */
    const Stencil& stencil(unsigned timeIdx) const
    { return elemCtx_.stencil(timeIdx); }

    /*!
     * \brief Returns the outer unit normal of the boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    Vector normal(unsigned boundaryFaceIdx, unsigned timeIdx) const
    {
        auto tmp = stencil(timeIdx).boundaryFace[boundaryFaceIdx].normal;
        tmp /= tmp.two_norm();
        return tmp;
    }

    /*!
     * \brief Returns the area [m^2] of a given boudary segment.
     */
    Scalar boundarySegmentArea(unsigned boundaryFaceIdx, unsigned timeIdx) const
    { return elemCtx_.stencil(timeIdx).boundaryFace(boundaryFaceIdx).area(); }

    /*!
     * \brief Return the position of a local entity in global coordinates.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const GlobalPosition& pos(unsigned boundaryFaceIdx, unsigned timeIdx) const
    { return stencil(timeIdx).boundaryFace(boundaryFaceIdx).integrationPos(); }

    /*!
     * \brief Return the position of a control volume's center in global coordinates.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const GlobalPosition& cvCenter(unsigned boundaryFaceIdx, unsigned timeIdx) const
    {
        unsigned scvIdx = stencil(timeIdx).boundaryFace(boundaryFaceIdx).interiorIndex();
        return stencil(timeIdx).subControlVolume(scvIdx).globalPos();
    }

    /*!
     * \brief Return the local sub-control volume index upon which the linearization is
     *        currently focused.
     */
    unsigned focusDofIndex() const
    { return elemCtx_.focusDofIndex(); }

    /*!
     * \brief Return the local sub-control volume index of the
     *        interior of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    unsigned interiorScvIndex(unsigned boundaryFaceIdx, unsigned timeIdx) const
    { return stencil(timeIdx).boundaryFace(boundaryFaceIdx).interiorIndex(); }

    /*!
     * \brief Return the global space index of the sub-control volume
     *        at the interior of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    unsigned globalSpaceIndex(unsigned boundaryFaceIdx, unsigned timeIdx) const
    { return elemCtx_.globalSpaceIndex(interiorScvIndex(boundaryFaceIdx, timeIdx), timeIdx); }

    /*!
     * \brief Return the intensive quantities for the finite volume in the
     *        interiour of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const IntensiveQuantities& intensiveQuantities(unsigned boundaryFaceIdx, unsigned timeIdx) const
    {
        unsigned interiorScvIdx = this->interiorScvIndex(boundaryFaceIdx, timeIdx);
        return elemCtx_.intensiveQuantities(interiorScvIdx, timeIdx);
    }

    /*!
     * \brief Return the extensive quantities for a given boundary face.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const ExtensiveQuantities& extensiveQuantities(unsigned boundaryFaceIdx, unsigned timeIdx) const
    { return elemCtx_.boundaryExtensiveQuantities(boundaryFaceIdx, timeIdx); }

    /*!
     * \brief Return the intersection for the neumann segment
     *
     * TODO/HACK: The intersection should take a local index as an
     * argument. since that's not supported efficiently by the DUNE
     * grid interface, we just ignore the index argument here!
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     */
    const Intersection intersection(unsigned boundaryFaceIdx OPM_UNUSED) const
    { return *intersectionIt_; }

    /*!
     * \brief Return the intersection for the neumann segment
     *
     * TODO/HACK: the intersection iterator can basically be
     * considered as an index which is manipulated externally, but
     * context classes should not store any indices. it is done this
     * way for performance reasons
     */
    IntersectionIterator& intersectionIt()
    { return intersectionIt_; }

protected:
    const ElementContext& elemCtx_;
    IntersectionIterator intersectionIt_;
};

} // namespace Opm

#endif
