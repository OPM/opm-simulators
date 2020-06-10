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
 * \copydoc Opm::ImmiscibleRateVector
 */
#ifndef EWOMS_IMMISCIBLE_RATE_VECTOR_HH
#define EWOMS_IMMISCIBLE_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "immiscibleintensivequantities.hh"

namespace Opm {
/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Implements a vector representing rates of conserved quantities.
 *
 * This class is basically a Dune::FieldVector which can be set using
 * either mass, molar or volumetric rates.
 */
template <class TypeTag>
class ImmiscibleRateVector
    : public Dune::FieldVector<GetPropType<TypeTag, Properties::Evaluation>,
                               getPropValue<TypeTag, Properties::NumEq>()>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    using ParentType = Dune::FieldVector<Evaluation, numEq>;
    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;

public:
    /*!
     * \brief Default constructor
     */
    ImmiscibleRateVector() : ParentType()
    { Opm::Valgrind::SetUndefined(*this); }

    /*!
     * \brief Constructor with assignment from scalar
     *
     * \param value The scalar value to which all entries of the vector will be set.
     */
    ImmiscibleRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \brief Copy constructor
     *
     * \param value The rate vector that will be duplicated.
     */
    ImmiscibleRateVector(const ImmiscibleRateVector& value)
        : ParentType(value)
    {}

    /*!
     * \brief Set a mass rate of the conservation quantities.
     *
     * Enthalpy is _not_ taken into account seperately here. This
     * means that it must be set to the desired value in the
     * parameter.
     *
     * \param value The mass rate in \f$[kg/(m^2\,s)]\f$ (unit for areal fluxes)
     */
    void setMassRate(const ParentType& value)
    { ParentType::operator=(value); }

    /*!
     * \brief Set a molar rate of the conservation quantities.
     *
     * Enthalpy is _not_ taken into account seperately here. This
     * means that it must be set to the desired value in the
     * parameter.
     *
     * \param value The new molar rate in \f$[mol/(m^2\,s)]\f$
     */
    void setMolarRate(const ParentType& value)
    {
        // convert to mass rates
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            ParentType::operator[](conti0EqIdx + compIdx) =
                value[conti0EqIdx + compIdx]*FluidSystem::molarMass(compIdx);
    }

    /*!
     * \brief Set an enthalpy rate [J/As] where \f$A \in \{m^2, m^3\}\f$
     *
     * If the energy equation is not enabled, this method is a no-op.
     *
     * \param rate The enthalpy rate in \f$[J/(m^2\,s)]\f$
     */
    template <class RhsEval>
    void setEnthalpyRate(const RhsEval& rate)
    { EnergyModule::setEnthalpyRate(*this, rate); }

    /*!
     * \brief Set a volumetric rate of a phase.
     *
     * The enthalpy transported into the domain is taken into account
     * by this method.
     *
     * \param fluidState The thermodynamic state of the fluids which
     *                   should be considered. The density and the
     *                   composition of the considered phase must be
     *                   specified before calling this method.
     * \param phaseIdx The index of the fluid phase for which the
     *                 given amount of volume should be specified.
     * \param volume The volumetric rate of the fluid phase in
     *               \f$[m^3/(m^2\,s)]\f$ (unit for areal fluxes)
     */
    template <class FluidState, class RhsEval>
    void setVolumetricRate(const FluidState& fluidState, unsigned phaseIdx, const RhsEval& volume)
    {
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] =
                fluidState.density(phaseIdx, compIdx)
                * fluidState.massFraction(phaseIdx, compIdx)
                * volume;

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volume);
    }

    /*!
     * \brief Assignment operator from a scalar or a function evaluation
     */
    template <class RhsEval>
    ImmiscibleRateVector& operator=(const RhsEval& value)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = value;
        return *this;
    }

    /*!
     * \brief Assignment operator from another rate vector
     */
    ImmiscibleRateVector& operator=(const ImmiscibleRateVector& other)
    {
        for (unsigned i=0; i < this->size(); ++i)
            (*this)[i] = other[i];
        return *this;
    }
};

} // namespace Opm

#endif
