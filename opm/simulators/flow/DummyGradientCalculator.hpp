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
 * \copydoc Opm::DummyGradientCalculator
 */
#ifndef OPM_DUMMY_GRADIENT_CALCULATOR_HPP
#define OPM_DUMMY_GRADIENT_CALCULATOR_HPP

#include <dune/common/fvector.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <stdexcept>

namespace Opm {
/*!
 * \ingroup BlackOilSimulator
 *
 * \brief This is a "dummy" gradient calculator which does not do anything.
 *
 * The blackoil simulator does not need any gradients: Volume fluxes are calculated
 * via pressure differences instead of pressure gradients (i.e., transmissibilities
 * instead of permeabilities), and an energy equation and molecular diffusion are not
 * supported.
 */
template<class TypeTag>
class DummyGradientCalculator
{
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    enum { dimWorld = GridView::dimensionworld };

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    static void registerParameters()
    { }

    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext&, unsigned)
    { }

    template <class QuantityCallback, class QuantityType = Scalar>
    QuantityType calculateValue(const ElementContext&,
                                unsigned,
                                const QuantityCallback&) const
    {
        throw std::logic_error("Generic values are not supported by the black-oil simulator");
    }

    template <class QuantityCallback>
    void calculateGradient(DimVector&,
                           const ElementContext&,
                           unsigned,
                           const QuantityCallback&) const
    {
        throw std::logic_error("Generic gradients are not supported by the black-oil simulator");
    }

    template <class QuantityCallback>
    Scalar calculateBoundaryValue(const ElementContext&,
                                  unsigned,
                                  const QuantityCallback&)
    {
        throw std::logic_error("Generic boundary values are not supported by the black-oil simulator");
    }

    template <class QuantityCallback>
    void calculateBoundaryGradient(DimVector&,
                                   const ElementContext&,
                                   unsigned,
                                   const QuantityCallback&) const
    {
        throw std::logic_error("Generic boundary gradients are not supported by the black-oil simulator");
    }
};
} // namespace Opm

#endif // OPM_DUMMY_GRADIENT_CALCULATOR_HPP
