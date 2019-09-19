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
 * \copydoc Opm::RichardsRateVector
 */
#ifndef EWOMS_RICHARDS_RATE_VECTOR_HH
#define EWOMS_RICHARDS_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "richardsintensivequantities.hh"

namespace Opm {

/*!
 * \ingroup RichardsModel
 *
 * \brief Implements a vector representing mass, molar or volumetric rates.
 *
 * This class is basically a Dune::FieldVector which can be set using either mass, molar
 * or volumetric rates.
 */
template <class TypeTag>
class RichardsRateVector
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Evaluation),
                               GET_PROP_VALUE(TypeTag, NumEq)>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) EnergyModule;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidCompIdx = GET_PROP_VALUE(TypeTag, LiquidComponentIndex) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Evaluation, numEq> ParentType;

public:
    RichardsRateVector() : ParentType()
    { Opm::Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(Scalar)
     */
    RichardsRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(const
     * ImmiscibleRateVector& )
     */
    RichardsRateVector(const RichardsRateVector& value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::setMassRate
     */
    void setMassRate(const ParentType& value)
    { ParentType::operator=(value); }

    /*!
     * \copydoc ImmiscibleRateVector::setMolarRate
     */
    void setMolarRate(const ParentType& value)
    {
        // convert to mass rates
        ParentType::operator[](contiEqIdx) =
            value[contiEqIdx]*FluidSystem::molarMass(liquidCompIdx);
    }

    /*!
     * \copydoc ImmiscibleRateVector::setEnthalpyRate
     */
    template <class RhsEval>
    void setEnthalpyRate(const RhsEval& rate)
    { EnergyModule::setEnthalpyRate(*this, rate); }

    /*!
     * \copydoc ImmiscibleRateVector::setVolumetricRate
     */
    template <class FluidState, class RhsEval>
    void setVolumetricRate(const FluidState& fluidState, unsigned phaseIdx, const RhsEval& volume)
    {
       (*this)[contiEqIdx] =
            fluidState.density(phaseIdx)
            * fluidState.massFraction(phaseIdx, liquidCompIdx)
            * volume;

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volume);
    }


    /*!
     * \brief Assignment operator from a scalar or a function evaluation
     */
    template <class RhsEval>
    RichardsRateVector& operator=(const RhsEval& value)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = value;
        return *this;
    }

    /*!
     * \brief Assignment operator from another rate vector
     */
    RichardsRateVector& operator=(const RichardsRateVector& other)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = other[i];
        return *this;
    }
};

} // namespace Opm

#endif
