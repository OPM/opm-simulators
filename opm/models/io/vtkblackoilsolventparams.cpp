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
#include <opm/models/io/vtkblackoilsolventparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void VtkBlackOilSolventParams::registerParameters()
{
    Parameters::Register<Parameters::VtkWriteSolventSaturation>
        ("Include the \"saturation\" of the solvent component "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSolventRsw>
        ("Include the \"dissolved volume in water\" of the solvent component "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSolventDensity>
        ("Include the \"density\" of the solvent component "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSolventViscosity>
        ("Include the \"viscosity\" of the solvent component "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSolventMobility>
        ("Include the \"mobility\" of the solvent component "
         "in the VTK output files");
}

void VtkBlackOilSolventParams::read()
{
    solventSaturationOutput_ = Parameters::Get<Parameters::VtkWriteSolventSaturation>();
    solventRswOutput_ = Parameters::Get<Parameters::VtkWriteSolventRsw>();
    solventDensityOutput_ = Parameters::Get<Parameters::VtkWriteSolventDensity>();
    solventViscosityOutput_ = Parameters::Get<Parameters::VtkWriteSolventViscosity>();
    solventMobilityOutput_ = Parameters::Get<Parameters::VtkWriteSolventMobility>();
}

} // namespace Opm
