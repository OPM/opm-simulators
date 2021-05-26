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
 * \copydoc Opm::FvBaseGradientCalculator
 */
#ifndef EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH
#define EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH

#include "fvbaseproperties.hh"

#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Opm {
template<class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief This class calculates gradients of arbitrary quantities at
 *        flux integration points using the two-point approximation scheme
 */
template<class TypeTag>
class FvBaseGradientCalculator
{
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // maximum number of flux approximation points. to calculate this,
    // we assume that the geometry with the most pointsq is a cube.
    enum { maxFap = 2 << dim };

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using EvalDimVector = Dune::FieldVector<Evaluation, dimWorld>;

public:
    /*!
     * \brief Register all run-time parameters for the gradient calculator
     *        of the base class of the discretization.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Precomputes the common values to calculate gradients and values of
     *        quantities at every interior flux approximation point.
     *
     * \param elemCtx The current execution context
     * \param timeIdx The index used by the time discretization.
     */
    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext& elemCtx OPM_UNUSED, unsigned timeIdx OPM_UNUSED)
    { /* noting to do */ }

    /*!
     * \brief Calculates the value of an arbitrary scalar quantity at any interior flux
     *        approximation point.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point in the current
     *               element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    auto calculateScalarValue(const ElementContext& elemCtx,
                              unsigned fapIdx,
                              const QuantityCallback& quantityCallback) const
        -> typename std::remove_reference<decltype(quantityCallback.operator()(0))>::type
    {
        using RawReturnType = decltype(quantityCallback.operator()(0));
        using ReturnType = typename std::remove_const<typename std::remove_reference<RawReturnType>::type>::type;

        Scalar interiorDistance;
        Scalar exteriorDistance;
        computeDistances_(interiorDistance, exteriorDistance, elemCtx, fapIdx);

        const auto& face = elemCtx.stencil(/*timeIdx=*/0).interiorFace(fapIdx);
        auto i = face.interiorIndex();
        auto j = face.exteriorIndex();
        auto focusDofIdx = elemCtx.focusDofIndex();

        // use the average weighted by distance...
        ReturnType value;
        if (i == focusDofIdx)
            value = quantityCallback(i)*interiorDistance;
        else
            value = getValue(quantityCallback(i))*interiorDistance;

        if (j == focusDofIdx)
            value += quantityCallback(j)*exteriorDistance;
        else
            value += getValue(quantityCallback(j))*exteriorDistance;

        value /= interiorDistance + exteriorDistance;

        return value;
    }

    /*!
     * \brief Calculates the value of an arbitrary vectorial quantity at any interior flux
     *        approximation point.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point in the current
     *               element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    auto calculateVectorValue(const ElementContext& elemCtx,
                              unsigned fapIdx,
                              const QuantityCallback& quantityCallback) const
        -> typename std::remove_reference<decltype(quantityCallback.operator()(0))>::type
    {
        using RawReturnType = decltype(quantityCallback.operator()(0));
        using ReturnType = typename std::remove_const<typename std::remove_reference<RawReturnType>::type>::type;

        Scalar interiorDistance;
        Scalar exteriorDistance;
        computeDistances_(interiorDistance, exteriorDistance, elemCtx, fapIdx);

        const auto& face = elemCtx.stencil(/*timeIdx=*/0).interiorFace(fapIdx);
        auto i = face.interiorIndex();
        auto j = face.exteriorIndex();
        auto focusDofIdx = elemCtx.focusDofIndex();

        // use the average weighted by distance...
        ReturnType value;
        if (i == focusDofIdx) {
            value = quantityCallback(i);
            for (int k = 0; k < value.size(); ++k)
                value[k] *= interiorDistance;
        }
        else {
            const auto& dofVal = getValue(quantityCallback(i));
            for (int k = 0; k < dofVal.size(); ++k)
                value[k] = getValue(dofVal[k])*interiorDistance;
        }

        if (j == focusDofIdx) {
            const auto& dofVal = quantityCallback(j);
            for (int k = 0; k < dofVal.size(); ++k)
                value[k] += dofVal[k]*exteriorDistance;
        }
        else {
            const auto& dofVal = quantityCallback(j);
            for (int k = 0; k < dofVal.size(); ++k)
                value[k] += getValue(dofVal[k])*exteriorDistance;
        }

        Scalar totDistance = interiorDistance + exteriorDistance;
        for (int k = 0; k < value.size(); ++k)
            value[k] /= totDistance;

        return value;
    }

    /*!
     * \brief Calculates the gradient of an arbitrary quantity at any
     *        flux approximation point.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity given the index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateGradient(EvalDimVector& quantityGrad,
                           const ElementContext& elemCtx,
                           unsigned fapIdx,
                           const QuantityCallback& quantityCallback) const
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.interiorFace(fapIdx);

        auto i = face.interiorIndex();
        auto j = face.exteriorIndex();
        auto focusIdx = elemCtx.focusDofIndex();

        const auto& interiorPos = stencil.subControlVolume(i).globalPos();
        const auto& exteriorPos = stencil.subControlVolume(j).globalPos();

        Evaluation deltay;
        if (i == focusIdx) {
            deltay =
                getValue(quantityCallback(j))
                - quantityCallback(i);
        }
        else if (j == focusIdx) {
            deltay =
                quantityCallback(j)
                - getValue(quantityCallback(i));
        }
        else
            deltay =
                getValue(quantityCallback(j))
                - getValue(quantityCallback(i));

        Scalar distSquared = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
            distSquared += tmp*tmp;
        }

        // divide the gradient by the squared distance between the centers of the
        // sub-control volumes: the gradient is the normalized directional vector between
        // the two centers times the ratio of the difference of the values and their
        // distance, i.e., d/abs(d) * delta y / abs(d) = d*delta y / abs(d)^2.
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
            quantityGrad[dimIdx] = deltay*(tmp/distSquared);
        }
    }

    /*!
     * \brief Calculates the value of an arbitrary quantity at any
     *        flux approximation point on the grid boundary.
     *
     * Boundary values are always calculated using the two-point
     * approximation.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity given the index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    auto calculateBoundaryValue(const ElementContext& elemCtx OPM_UNUSED,
                                unsigned fapIdx OPM_UNUSED,
                                const QuantityCallback& quantityCallback)
        -> decltype(quantityCallback.boundaryValue())
    { return quantityCallback.boundaryValue(); }

    /*!
     * \brief Calculates the gradient of an arbitrary quantity at any
     *        flux approximation point on the boundary.
     *
     * Boundary gradients are always calculated using the two-point
     * approximation.
     *
     * \param elemCtx The current execution context
     * \param faceIdx The local index of the flux approximation point
     *                in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateBoundaryGradient(EvalDimVector& quantityGrad,
                                   const ElementContext& elemCtx,
                                   unsigned faceIdx,
                                   const QuantityCallback& quantityCallback) const
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.boundaryFace(faceIdx);

        Evaluation deltay;
        if (face.interiorIndex() == elemCtx.focusDofIndex())
            deltay = quantityCallback.boundaryValue() - quantityCallback(face.interiorIndex());
        else
            deltay =
                getValue(quantityCallback.boundaryValue())
                - getValue(quantityCallback(face.interiorIndex()));

        const auto& boundaryFacePos = face.integrationPos();
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).center();

        Scalar distSquared = 0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = boundaryFacePos[dimIdx] - interiorPos[dimIdx];
            distSquared += tmp*tmp;
        }

        // divide the gradient by the squared distance between the center of the
        // sub-control and the center of the boundary face: the gradient is the
        // normalized directional vector between the two centers times the ratio of the
        // difference of the values and their distance, i.e., d/abs(d) * deltay / abs(d)
        // = d*deltay / abs(d)^2.
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = boundaryFacePos[dimIdx] - interiorPos[dimIdx];
            quantityGrad[dimIdx] = deltay*(tmp/distSquared);
        }
    }

private:
    void computeDistances_(Scalar& interiorDistance,
                           Scalar& exteriorDistance,
                           const ElementContext& elemCtx,
                           unsigned fapIdx) const
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.interiorFace(fapIdx);

        // calculate the distances of the position of the interior and of the exterior
        // finite volume to the position of the integration point.
        const auto& normal = face.normal();
        auto i = face.interiorIndex();
        auto j = face.exteriorIndex();
        const auto& interiorPos = stencil.subControlVolume(i).globalPos();
        const auto& exteriorPos = stencil.subControlVolume(j).globalPos();
        const auto& integrationPos = face.integrationPos();

        interiorDistance = 0.0;
        exteriorDistance = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            interiorDistance +=
                (interiorPos[dimIdx] - integrationPos[dimIdx])
                * normal[dimIdx];

            exteriorDistance +=
                (exteriorPos[dimIdx] - integrationPos[dimIdx])
                * normal[dimIdx];
        }

        interiorDistance = std::sqrt(std::abs(interiorDistance));
        exteriorDistance = std::sqrt(std::abs(exteriorDistance));
    }
};
} // namespace Opm

#endif
