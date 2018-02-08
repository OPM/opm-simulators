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
 * \copydoc Ewoms::EclDummyGradientCalculator
 */
#ifndef EWOMS_ECL_DUMMY_GRADIENT_CALCULATOR_HH
#define EWOMS_ECL_DUMMY_GRADIENT_CALCULATOR_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

namespace Ewoms {
/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This is a "dummy" gradient calculator which does not do anything.
 *
 * The ECL blackoil simulator does not need any gradients: Volume fluxes are calculated
 * via pressure differences instead of pressure gradients (i.e., transmissibilities
 * instead of permeabilities), and an energy equation and molecular diffusion are not
 * supported.
 */
template<class TypeTag>
class EclDummyGradientCalculator
{

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    static void registerParameters()
    { }

    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext& elemCtx OPM_UNUSED, unsigned timeIdx OPM_UNUSED)
    { }

    template <class QuantityCallback, class QuantityType = Scalar>
    QuantityType calculateValue(const ElementContext& elemCtx OPM_UNUSED,
                                unsigned fapIdx OPM_UNUSED,
                                const QuantityCallback& quantityCallback OPM_UNUSED) const
    {
        throw std::logic_error("Generic values are not supported by the ECL black-oil simulator");
    }

    template <class QuantityCallback>
    void calculateGradient(DimVector& quantityGrad OPM_UNUSED,
                           const ElementContext& elemCtx OPM_UNUSED,
                           unsigned fapIdx OPM_UNUSED,
                           const QuantityCallback& quantityCallback OPM_UNUSED) const
    {
        throw std::logic_error("Generic gradients are not supported by the ECL black-oil simulator");
    }

    template <class QuantityCallback>
    Scalar calculateBoundaryValue(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned fapIdx OPM_UNUSED,
                                  const QuantityCallback& quantityCallback OPM_UNUSED)
    {
        throw std::logic_error("Generic boundary values are not supported by the ECL black-oil simulator");
    }

    template <class QuantityCallback>
    void calculateBoundaryGradient(DimVector& quantityGrad OPM_UNUSED,
                                   const ElementContext& elemCtx OPM_UNUSED,
                                   unsigned fapIdx OPM_UNUSED,
                                   const QuantityCallback& quantityCallback OPM_UNUSED) const
    {
        throw std::logic_error("Generic boundary gradients are not supported by the ECL black-oil simulator");
    }
};
} // namespace Ewoms

#endif
