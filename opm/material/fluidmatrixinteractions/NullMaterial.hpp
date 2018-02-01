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
 * \copydoc Opm::NullMaterial
 */
#ifndef OPM_NULL_MATERIAL_HPP
#define OPM_NULL_MATERIAL_HPP

#include "NullMaterialParams.hpp"

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements a dummy linear saturation-capillary pressure
 *        relation which just disables capillary pressure.
 */
template <class TraitsT>
class NullMaterial : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef NullMaterialParams<TraitsT> Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const unsigned numPhases = Traits::numPhases;

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = (numPhases == 2);

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = (numPhases == 2);

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    //!
    //! In this law, the relative permeabilities are saturation
    //! dependent, even if capillary pressure is always zero
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    /*!
     * \brief Returns constant 0 for all phases.
     *
     * \param values Container for the return values
     * \param params Parameters
     * \param state The fluid state
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT& values,
                                   const Params& /*params*/,
                                   const FluidState& /*fluidState*/)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            values[phaseIdx] = 0.0;
    }

    /*!
     * \brief The inverse of the capillary pressure
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT& /*values*/,
                            const Params& /*params*/,
                            const FluidState& /*fluidState*/)
    { throw std::logic_error("Not defined: NullMaterial::saturations()"); }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT& values,
                                       const Params& /*params*/,
                                       const FluidState& fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Evaluation& S =
                Opm::decay<Evaluation>(fluidState.saturation(phaseIdx));
            values[phaseIdx] = Opm::max(Opm::min(S, 1.0), 0.0);
        }
    }

    /*!
     * \brief The difference between the pressures of the non-wetting and wetting phase.
     */
    template <class FluidState, class Evaluation  = typename FluidState::Scalar>
    static typename std::enable_if<(numPhases > 1), Evaluation>::type
    pcnw(const Params& /*params*/, const FluidState& /*fluidState*/)
    { return 0; }

    template <class Evaluation>
    static typename std::enable_if<numPhases == 2, Evaluation>::type
    twoPhaseSatPcnw(const Params& /*params*/, const Evaluation& /*Sw*/)
    { return 0; }

    /*!
     * \brief Calculate wetting phase saturation given that the rest
     *        of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Scalar Sw(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not defined: Sw()"); }

    template <class Evaluation>
    static typename std::enable_if<numPhases == 2, Evaluation>::type
    twoPhaseSatSw(const Params& /*params*/, const Evaluation& /*pcnw*/)
    { throw std::logic_error("Not defined: twoPhaseSatSw()"); }

    /*!
     * \brief Calculate non-wetting phase saturation given that the
     *        rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Scalar Sn(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not defined: Sn()"); }

    template <class Evaluation>
    static typename std::enable_if<numPhases == 2, Evaluation>::type
    twoPhaseSatSn(const Params& /*params*/, const Evaluation& /*pcnw*/)
    { throw std::logic_error("Not defined: twoPhaseSatSn()"); }

    /*!
     * \brief Calculate gas phase saturation given that the rest of
     *        the fluid state has been initialized
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if< (numPhases > 2), Evaluation>::type
    Sg(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not defined: Sg()"); }

    /*!
     * \brief The relative permability of the wetting phase
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if<(numPhases > 1), Evaluation>::type
    krw(const Params& /*params*/, const FluidState& fluidState)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fluidState.saturation(Traits::wettingPhaseIdx));

        return Opm::max(0.0, Opm::min(1.0, Sw));
    }

    template <class Evaluation>
    static typename std::enable_if<numPhases == 2, Evaluation>::type
    twoPhaseSatKrw(const Params& /*params*/, const Evaluation& Sw)
    { return Opm::max(0.0, Opm::min(1.0, Sw)); }

    /*!
     * \brief The relative permability of the liquid non-wetting phase
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if<(numPhases > 1), Evaluation>::type
    krn(const Params& /*params*/, const FluidState& fluidState)
    {
        const Evaluation& Sn =
            Opm::decay<Evaluation>(fluidState.saturation(Traits::nonWettingPhaseIdx));

        return Opm::max(0.0, Opm::min(1.0, Sn));
    }

    template <class Evaluation>
    static typename std::enable_if<numPhases == 2, Evaluation>::type
    twoPhaseSatKrn(const Params& /*params*/, const Evaluation& Sw)
    {
        return Opm::max(0.0, Opm::min(1.0, 1.0 - Opm::decay<Evaluation>(Sw)));
    }

    /*!
     * \brief The relative permability of the gas phase
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if< (numPhases > 2), Evaluation>::type
    krg(const Params& /*params*/, const FluidState& fluidState)
    {
        const Evaluation& Sg =
            Opm::decay<Evaluation>(fluidState.saturation(Traits::gasPhaseIdx));

        return Opm::max(0.0, Opm::min(1.0, Sg));
    }

    /*!
     * \brief The difference between the pressures of the gas and the non-wetting phase.
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), Evaluation>::type
    pcgn(const Params& /*params*/, const FluidState& /*fluidState*/)
    { return 0.0; }
};
} // namespace Opm

#endif
