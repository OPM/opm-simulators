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
#include <opm/models/io/vtkblackoilbioeffectsparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void VtkBlackOilBioeffectsParams::registerParameters(const bool isMICP)
{
    Parameters::Register<Parameters::VtkWriteMicrobialConcentration>
        ("Include the concentration of the microbial component in the water phase "
         "in the VTK output files");
    Parameters::Register<Parameters::VtkWriteBiofilmVolumeFraction>
        ("Include the biofilm volume fraction in the VTK output files");
    if (isMICP) {
        Parameters::Register<Parameters::VtkWriteOxygenConcentration>
            ("Include the concentration of the oxygen component in the water phase "
            "in the VTK output files");
        Parameters::Register<Parameters::VtkWriteUreaConcentration>
            ("Include the concentration of the urea component in the water phase "
            "in the VTK output files");
        Parameters::Register<Parameters::VtkWriteCalciteVolumeFraction>
            ("Include the calcite volume fraction in the VTK output files");
    }
}

void VtkBlackOilBioeffectsParams::read(const bool isMICP)
{
    microbialConcentrationOutput_ = Parameters::Get<Parameters::VtkWriteMicrobialConcentration>();
    biofilmVolumeFractionOutput_ = Parameters::Get<Parameters::VtkWriteBiofilmVolumeFraction>();
    if (isMICP) {
        oxygenConcentrationOutput_ = Parameters::Get<Parameters::VtkWriteOxygenConcentration>();
        ureaConcentrationOutput_ = Parameters::Get<Parameters::VtkWriteUreaConcentration>();
        calciteVolumeFractionOutput_ = Parameters::Get<Parameters::VtkWriteCalciteVolumeFraction>();
    }
}

} // namespace Opm
