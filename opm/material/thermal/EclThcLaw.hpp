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
 * \copydoc Opm::EclThcLaw
 */
#ifndef OPM_ECL_THC_LAW_HPP
#define OPM_ECL_THC_LAW_HPP

#include "EclThcLawParams.hpp"

#include <opm/material/densead/Math.hpp>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements the total thermal conductivity and rock enthalpy relations used by ECL.
 *
 * This is the thermal conduction law based on the THCROCK, THCOIL, THCGAS and THCWATER
 * keywords.
 */
template <class ScalarT,
          class ParamsT = EclThcLawParams<ScalarT> >
class EclThcLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the total thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
      */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params,
                                          const FluidState&)
    {
        // The thermal conductivity approach based on the THC* keywords.

        // let's assume that the porosity of the rock at standard condition is meant
        Scalar poro = params.porosity();

        // IMO this approach is very questionable because the total thermal conductivity
        // should at least depend on the current solution's phase saturation. Since ECL
        // is king, let's follow their lead and throw ourselfs down the cliff of obvious
        // incorrectness.
        //
        // TODO: also follow their fine leadership in the twophase case.
        Scalar numPhases = 3.0;
        Scalar thconAvg =
            poro*(params.thcoil() + params.thcgas() + params.thcwater()) / numPhases
            + (1.0 - poro)*params.thcrock();

        return thconAvg;
    }
};
} // namespace Opm

#endif
