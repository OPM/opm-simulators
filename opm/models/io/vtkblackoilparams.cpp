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

#include <config.h>
#include <opm/models/io/vtkblackoilparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void VtkBlackoilParams::registerParameters()
{
    Parameters::Register<Parameters::VtkWriteGasDissolutionFactor>
        ("Include the gas dissolution factor (R_s) of the observed oil "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteOilVaporizationFactor>
        ("Include the oil vaporization factor (R_v) of the observed gas "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteOilFormationVolumeFactor>
        ("Include the oil formation volume factor (B_o) in the VTK output files");
    Parameters::Register<Parameters::VtkWriteGasFormationVolumeFactor>
        ("Include the gas formation volume factor (B_g) in the "
         "VTK output files");
    Parameters::Register<Parameters::VtkWriteWaterFormationVolumeFactor>
        ("Include the water formation volume factor (B_w) in the "
         "VTK output files");
    Parameters::Register<Parameters::VtkWriteOilSaturationPressure>
        ("Include the saturation pressure of oil (p_o,sat) in the "
         "VTK output files");
    Parameters::Register<Parameters::VtkWriteGasSaturationPressure>
        ("Include the saturation pressure of gas (p_g,sat) in the "
         "VTK output files");
    Parameters::Register<Parameters::VtkWriteSaturatedOilGasDissolutionFactor>
        ("Include the gas dissolution factor (R_s,sat) of gas saturated "
         "oil in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSaturatedGasOilVaporizationFactor>
        ("Include the oil vaporization factor (R_v,sat) of oil saturated "
         "gas in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSaturationRatios>
        ("Write the ratio of the actually and maximum dissolved component of "
         "the mixtures");
    Parameters::Register<Parameters::VtkWritePrimaryVarsMeaning>
        ("Include how the primary variables should be interpreted to the "
         "VTK output files");
}

void VtkBlackoilParams::read()
{
    gasDissolutionFactorOutput_ = Parameters::Get<Parameters::VtkWriteGasDissolutionFactor>();
    oilVaporizationFactorOutput_ = Parameters::Get<Parameters::VtkWriteOilVaporizationFactor>();
    oilFormationVolumeFactorOutput_ = Parameters::Get<Parameters::VtkWriteOilFormationVolumeFactor>();
    gasFormationVolumeFactorOutput_ = Parameters::Get<Parameters::VtkWriteGasFormationVolumeFactor>();
    waterFormationVolumeFactorOutput_ = Parameters::Get<Parameters::VtkWriteWaterFormationVolumeFactor>();
    oilSaturationPressureOutput_ = Parameters::Get<Parameters::VtkWriteOilSaturationPressure>();
    gasSaturationPressureOutput_ = Parameters::Get<Parameters::VtkWriteGasSaturationPressure>();
    saturatedOilGasDissolutionFactorOutput_ = Parameters::Get<Parameters::VtkWriteSaturatedOilGasDissolutionFactor>();
    saturatedGasOilVaporizationFactorOutput_ = Parameters::Get<Parameters::VtkWriteSaturatedGasOilVaporizationFactor>();
    saturationRatiosOutput_ = Parameters::Get<Parameters::VtkWriteSaturationRatios>();
    primaryVarsMeaningOutput_ = Parameters::Get<Parameters::VtkWritePrimaryVarsMeaning>();
}

} // namespace Opm
