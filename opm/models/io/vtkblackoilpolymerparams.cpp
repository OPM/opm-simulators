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
#include <opm/models/io/vtkblackoilpolymerparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void VtkBlackoilPolymerParams::registerParameters()
{
    Parameters::Register<Parameters::VtkWritePolymerConcentration>
        ("Include the concentration of the polymer component in the water phase "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWritePolymerDeadPoreVolume>
        ("Include the fraction of the \"dead\" pore volume "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWritePolymerRockDensity>
        ("Include the amount of already adsorbed polymer component"
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWritePolymerAdsorption>
        ("Include the adsorption rate of the polymer component"
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWritePolymerViscosityCorrection>
        ("Include the viscosity correction of the polymer component "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteWaterViscosityCorrection>
        ("Include the viscosity correction of the water component "
         "due to polymers in the VTK output files");
}

void VtkBlackoilPolymerParams::read()
{
    polymerConcentrationOutput_ = Parameters::Get<Parameters::VtkWritePolymerConcentration>();
    polymerDeadPoreVolumeOutput_ = Parameters::Get<Parameters::VtkWritePolymerDeadPoreVolume>();
    polymerRockDensityOutput_ = Parameters::Get<Parameters::VtkWritePolymerRockDensity>();
    polymerAdsorptionOutput_ = Parameters::Get<Parameters::VtkWritePolymerAdsorption>();
    polymerViscosityCorrectionOutput_ = Parameters::Get<Parameters::VtkWritePolymerViscosityCorrection>();
    waterViscosityCorrectionOutput_ = Parameters::Get<Parameters::VtkWritePolymerViscosityCorrection>();
}

} // namespace Opm
