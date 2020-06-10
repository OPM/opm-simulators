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
 * \copydoc Opm::FlashRateVector
 */
#ifndef EWOMS_FLASH_RATE_VECTOR_HH
#define EWOMS_FLASH_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/models/common/energymodule.hh>
#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/common/Valgrind.hpp>

#include "flashintensivequantities.hh"

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \copydoc Opm::ImmiscibleRateVector
 */
template <class TypeTag>
class FlashRateVector
    : public Dune::FieldVector<GetPropType<TypeTag, Properties::Evaluation>,
                               getPropValue<TypeTag, Properties::NumEq>()>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };

    using ParentType = Dune::FieldVector<Evaluation, numEq>;
    using EnergyModule = Opm::EnergyModule<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>;

public:
    FlashRateVector() : ParentType()
    { Opm::Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(Scalar)
     */
    FlashRateVector(const Evaluation& value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(const
     * ImmiscibleRateVector& )
     */
    FlashRateVector(const FlashRateVector& value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::setMassRate
     */
    void setMassRate(const ParentType& value)
    {
        // convert to molar rates
        ParentType molarRate(value);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            molarRate[conti0EqIdx + compIdx] /= FluidSystem::molarMass(compIdx);

        setMolarRate(molarRate);
    }

    /*!
     * \copydoc ImmiscibleRateVector::setMolarRate
     */
    void setMolarRate(const ParentType& value)
    { ParentType::operator=(value); }

    /*!
     * \copydoc ImmiscibleRateVector::setEnthalpyRate
     */
    void setEnthalpyRate(const Evaluation& rate)
    { EnergyModule::setEnthalpyRate(*this, rate); }

    /*!
     * \copydoc ImmiscibleRateVector::setVolumetricRate
     */
    template <class FluidState, class RhsEval>
    void setVolumetricRate(const FluidState& fluidState, unsigned phaseIdx, const RhsEval& volume)
    {
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] =
                fluidState.density(phaseIdx, compIdx)
                * fluidState.moleFraction(phaseIdx, compIdx)
                * volume;

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volume);
    }

    /*!
     * \brief Assignment operator from a scalar or a function evaluation
     */
    template <class RhsEval>
    FlashRateVector& operator=(const RhsEval& value)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = value;
        return *this;
    }

    /*!
     * \brief Assignment operator from another rate vector
     */
    FlashRateVector& operator=(const FlashRateVector& other)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = other[i];
        return *this;
    }
};

} // namespace Opm

#endif
