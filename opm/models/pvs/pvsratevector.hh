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
 * \copydoc Opm::PvsRateVector
 */
#ifndef EWOMS_PVS_RATE_VECTOR_HH
#define EWOMS_PVS_RATE_VECTOR_HH

#include "pvsindices.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup PvsModel
 *
 * \brief Implements a vector representing molar, mass or volumetric rates.
 *
 * This class is basically a Dune::FieldVector which can be set using either mass, molar
 * or volumetric rates.
 */
template <class TypeTag>
class PvsRateVector
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Evaluation),
                               GET_PROP_VALUE(TypeTag, NumEq)>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Dune::FieldVector<Evaluation, numEq> ParentType;
    typedef Opm::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    PvsRateVector() : ParentType()
    { Opm::Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(Scalar)
     */
    PvsRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(const ImmiscibleRateVector&)
     */
    PvsRateVector(const PvsRateVector& value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::setMassRate
     */
    void setMassRate(const ParentType& value)
    {
        // convert to molar rates
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            (*this)[conti0EqIdx + compIdx] = value[conti0EqIdx + compIdx];
            (*this)[conti0EqIdx + compIdx] /= FluidSystem::molarMass(compIdx);
        }
    }

    /*!
     * \copydoc ImmiscibleRateVector::setMolarRate
     */
    void setMolarRate(const ParentType& value)
    { ParentType::operator=(value); }

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
    PvsRateVector& operator=(const RhsEval& value)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = value;
        return *this;
    }

    /*!
     * \brief Assignment operator from another rate vector
     */
    PvsRateVector& operator=(const PvsRateVector& other)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = other[i];
        return *this;
    }
};

} // namespace Opm

#endif
