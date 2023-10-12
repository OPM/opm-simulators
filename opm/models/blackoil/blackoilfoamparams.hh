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
 * \brief Contains the parameters to extend the black-oil model to include the effects of foam.
 */
#ifndef EWOMS_BLACK_OIL_FOAM_PARAMS_HH
#define EWOMS_BLACK_OIL_FOAM_PARAMS_HH

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#include <vector>

namespace Opm {

//! \brief Struct holding the parameters for the BlackoilFoamModule class.
template<class Scalar>
struct BlackOilFoamParams {
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    /*!
     * \brief Specify the number of saturation regions.
     */
    void setNumSatRegions(unsigned numRegions)
    {
        foamCoefficients_.resize(numRegions);
        foamRockDensity_.resize(numRegions);
        foamAllowDesorption_.resize(numRegions);
        adsorbedFoamTable_.resize(numRegions);
    }

    // a struct containing constants to calculate change to relative permeability,
    // based on model (1-9) in Table 1 of
    // Kun Ma, Guangwei Ren, Khalid Mateen, Danielle Morel, and Philippe Cordelier:
    // "Modeling techniques for foam flow in porous media", SPE Journal, 20(03):453–470, jun 2015.
    // The constants are provided by various deck keywords as shown in the comments below.
    struct FoamCoefficients {
        Scalar fm_min = 1e-20;   // FOAMFSC
        Scalar fm_mob = 1.0;     // FOAMFRM

        Scalar fm_surf = 1.0;    // FOAMFSC
        Scalar ep_surf = 1.0;    // FOAMFSC

        Scalar fm_oil = 1.0;     // FOAMFSO
        Scalar fl_oil = 0.0;     // FOAMFSO
        Scalar ep_oil = 0.0;     // FOAMFSO

        Scalar fm_cap = 1.0;     // FOAMFCN
        Scalar ep_cap = 0.0;     // FOAMFCN

        Scalar fm_dry = 1.0;     // FOAMFSW
        Scalar ep_dry = 0.0;     // FOAMFSW
    };

    std::vector<Scalar> foamRockDensity_;
    std::vector<bool> foamAllowDesorption_;
    std::vector<FoamCoefficients> foamCoefficients_;
    std::vector<TabulatedFunction> adsorbedFoamTable_;
    std::vector<TabulatedFunction> gasMobilityMultiplierTable_;
    Opm::Phase transport_phase_;
};

} // namespace Opm

#endif
