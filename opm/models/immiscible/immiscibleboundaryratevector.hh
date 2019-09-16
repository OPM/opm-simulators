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
 * \copydoc Opm::ImmiscibleBoundaryRateVector
 */
#ifndef EWOMS_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "immiscibleintensivequantities.hh"

namespace Opm {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Implements a boundary vector for the fully implicit
 *        multi-phase model which assumes immiscibility.
 */
template <class TypeTag>
class ImmiscibleBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Opm::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    ImmiscibleBoundaryRateVector()
        : ParentType()
    {}

    /*!
     * \brief Constructor that assigns all entries to a scalar value.
     *
     * \param value The scalar value to which all components of the
     *              boundary rate vector will be set.
     */
    ImmiscibleBoundaryRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \brief Copy constructor
     *
     * \param value The boundary rate vector to be duplicated.
     */
    ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector& value) = default;

    ImmiscibleBoundaryRateVector& operator=(const ImmiscibleBoundaryRateVector& value) = default;

    /*!
     * \brief Specify a free-flow boundary
     *
     * \param context The execution context for which the boundary rate should
     *                be specified.
     * \param bfIdx The local space index of the boundary segment.
     * \param timeIdx The index used by the time discretization.
     * \param fluidState The repesentation of the thermodynamic state
     *                   of the system on the integration point of the
     *                   boundary segment.
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context, unsigned bfIdx, unsigned timeIdx, const FluidState& fluidState)
    {
        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);
        unsigned focusDofIdx = context.focusDofIndex();
        unsigned interiorDofIdx = context.interiorScvIndex(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = Evaluation(0.0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto& pBoundary = fluidState.pressure(phaseIdx);
            const Evaluation& pInside = insideIntQuants.fluidState().pressure(phaseIdx);

            // mass conservation
            Evaluation density;
            if  (pBoundary > pInside) {
                if (focusDofIdx == interiorDofIdx)
                    density = fluidState.density(phaseIdx);
                else
                    density = Opm::getValue(fluidState.density(phaseIdx));
            }
            else if (focusDofIdx == interiorDofIdx)
                density = insideIntQuants.fluidState().density(phaseIdx);
            else
                density = Opm::getValue(insideIntQuants.fluidState().density(phaseIdx));

            Opm::Valgrind::CheckDefined(density);
            Opm::Valgrind::CheckDefined(extQuants.volumeFlux(phaseIdx));

            (*this)[conti0EqIdx + phaseIdx] += extQuants.volumeFlux(phaseIdx)*density;

            // energy conservation
            if (enableEnergy) {
                Evaluation specificEnthalpy;
                if (pBoundary > pInside) {
                    if (focusDofIdx == interiorDofIdx)
                        specificEnthalpy = fluidState.enthalpy(phaseIdx);
                    else
                        specificEnthalpy = Opm::getValue(fluidState.enthalpy(phaseIdx));
                }
                else if (focusDofIdx == interiorDofIdx)
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);
                else
                    specificEnthalpy = Opm::getValue(insideIntQuants.fluidState().enthalpy(phaseIdx));

                Evaluation enthalpyRate = density*extQuants.volumeFlux(phaseIdx)*specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        // thermal conduction
        EnergyModule::addToEnthalpyRate(*this, EnergyModule::thermalConductionRate(extQuants));

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i)
            Opm::Valgrind::CheckDefined((*this)[i]);
        Opm::Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \brief Specify an inflow boundary
     *
     * An inflow boundary condition is basically a free flow boundary
     * condition that is not prevented from specifying a flow out of
     * the domain.
     *
     * \copydetails setFreeFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::min(0.0, val);
        }
    }

    /*!
     * \brief Specify an outflow boundary
     *
     * An outflow boundary condition is basically a free flow boundary
     * condition that is not prevented from specifying a flow into
     * the domain.
     *
     * \copydetails setFreeFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context& context,
                    unsigned bfIdx,
                    unsigned timeIdx,
                    const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::max(0.0, val);
        }
    }

    /*!
     * \brief Specify a no-flow boundary for all conserved quantities.
     */
    void setNoFlow()
    { (*this) = Evaluation(0.0); }
};

} // namespace Opm

#endif
