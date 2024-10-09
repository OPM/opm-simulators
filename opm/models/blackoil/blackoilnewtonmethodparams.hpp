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
 * \copydoc Opm::BlackOilNewtonMethod
 */
#ifndef EWOMS_BLACK_OIL_NEWTON_METHOD_PARAMETERS_HH
#define EWOMS_BLACK_OIL_NEWTON_METHOD_PARAMETERS_HH

namespace Opm::Parameters {

template<class Scalar>
struct DpMaxRel { static constexpr Scalar value = 0.3; };

template<class Scalar>
struct DsMax { static constexpr Scalar value = 0.2; };

template<class Scalar>
struct PriVarOscilationThreshold { static constexpr Scalar value = 1e-5; };

struct ProjectSaturations { static constexpr bool value = false; };

template<class Scalar>
struct MaxTemperatureChange { static constexpr Scalar value = 5.0; }; // Kelvin

template<class Scalar>
struct TemperatureMax { static constexpr Scalar value = 1e9; }; // Kelvin

template<class Scalar>
struct TemperatureMin { static constexpr Scalar value = 0.0; }; // Kelvin

template<class Scalar>
struct PressureMax { static constexpr Scalar value = 1e99; };

template<class Scalar>
struct PressureMin { static constexpr Scalar value = -1e99; };

template<class Scalar>
struct MaximumWaterSaturation { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct WaterOnlyThreshold { static constexpr Scalar value = 1.0; };

} // namespace Opm::Parameters

namespace Opm {

/*!
 * \brief Struct holding the parameters for BlackoilNewtonMethod.
 */
template<class Scalar>
struct BlackoilNewtonParams
{
    //! \brief Registers the parameters in parameter system.
    static void registerParameters();

    //! \brief Reads the parameter values from the parameter system.
    void read();

    Scalar priVarOscilationThreshold_;
    Scalar waterSaturationMax_;
    Scalar waterOnlyThreshold_;

    Scalar dpMaxRel_;
    Scalar dsMax_;
    bool projectSaturations_;
    Scalar maxTempChange_;
    Scalar tempMax_;
    Scalar tempMin_;
    Scalar pressMax_;
    Scalar pressMin_;
};

} // namespace Opm

#endif
