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
 * \copydoc Opm::VtkBlackOilModule
 */
#ifndef OPM_VTK_BLACK_OIL_PARAMS_HPP
#define OPM_VTK_BLACK_OIL_PARAMS_HPP

namespace Opm::Parameters {

// set default values for what quantities to output
struct VtkWriteGasDissolutionFactor { static constexpr bool value = false; };
struct VtkWriteOilVaporizationFactor { static constexpr bool value = false; };
struct VtkWriteOilFormationVolumeFactor { static constexpr bool value = false; };
struct VtkWriteGasFormationVolumeFactor { static constexpr bool value = false; };
struct VtkWriteWaterFormationVolumeFactor { static constexpr bool value = false; };
struct VtkWriteOilSaturationPressure { static constexpr bool value = false; };
struct VtkWriteGasSaturationPressure { static constexpr bool value = false; };
struct VtkWriteSaturationRatios { static constexpr bool value = false; };
struct VtkWriteSaturatedOilGasDissolutionFactor { static constexpr bool value = false; };
struct VtkWriteSaturatedGasOilVaporizationFactor { static constexpr bool value = false; };
struct VtkWritePrimaryVarsMeaning { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

/*!
 * \brief Struct holding the parameters for VtkBlackoilOutputModule.
 */
struct VtkBlackoilParams
{
    //! \brief Registers the parameters in parameter system.
    static void registerParameters();

    //! \brief Reads the parameter values from the parameter system.
    void read();

    bool gasDissolutionFactorOutput_;
    bool oilVaporizationFactorOutput_;
    bool oilFormationVolumeFactorOutput_;
    bool gasFormationVolumeFactorOutput_;
    bool waterFormationVolumeFactorOutput_;
    bool oilSaturationPressureOutput_;
    bool gasSaturationPressureOutput_;
    bool saturatedOilGasDissolutionFactorOutput_;
    bool saturatedGasOilVaporizationFactorOutput_;
    bool saturationRatiosOutput_;
    bool primaryVarsMeaningOutput_;
};

} // namespace Opm

#endif // OPM_VTK_BLACK_OIL_PARAMS_HPP
