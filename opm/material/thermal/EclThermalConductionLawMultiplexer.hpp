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
 * \copydoc Opm::EclThermalConductionLawMultiplexer
 */
#ifndef OPM_ECL_THERMAL_CONDUCTION_LAW_MULTIPLEXER_HPP
#define OPM_ECL_THERMAL_CONDUCTION_LAW_MULTIPLEXER_HPP

#include "EclThermalConductionLawMultiplexerParams.hpp"

#include "EclThconrLaw.hpp"
#include "EclThcLaw.hpp"
#include "NullThermalConductionLaw.hpp"

#include <opm/material/densead/Math.hpp>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements the total thermal conductivity and rock enthalpy relations used by ECL.
 */
template <class ScalarT,
          class FluidSystem,
          class ParamsT = EclThermalConductionLawMultiplexerParams<ScalarT>>
class EclThermalConductionLawMultiplexer
{
    enum { numPhases = FluidSystem::numPhases };

    typedef Opm::EclThconrLaw<ScalarT, FluidSystem, typename ParamsT::ThconrLawParams> ThconrLaw;
    typedef Opm::EclThcLaw<ScalarT, typename ParamsT::ThcLawParams> ThcLaw;
    typedef Opm::NullThermalConductionLaw<ScalarT> NullLaw;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, compute the volumetric internal energy of the rock [W/m^3].
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params,
                                       const FluidState& fluidState)
    {
        switch (params.thermalConductionApproach()) {
        case Params::thconrApproach:
            // relevant ECL keywords: THCONR and THCONSF
            return ThconrLaw::thermalConductivity(params.template getRealParams<Params::thconrApproach>(), fluidState);

        case Params::thcApproach:
            // relevant ECL keywords: THCROCK, THCOIL, THCGAS and THCWATER
            return ThcLaw::thermalConductivity(params.template getRealParams<Params::thcApproach>(), fluidState);

        case Params::nullApproach:
            // relevant ECL keywords: none or none recognized
            return NullLaw::thermalConductivity(0, fluidState);

        default:
            throw std::logic_error("Invalid thermal conductivity approach: "+std::to_string(int(params.thermalConductionApproach())));
        }
    }
};
} // namespace Opm

#endif
