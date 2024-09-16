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
#include <opm/models/io/vtkmultiphaseparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void VtkMultiPhaseParams::registerParameters()
{
    Parameters::Register<Parameters::VtkWriteExtrusionFactor>
        ("Include the extrusion factor of the degrees of freedom into the VTK output files");
    Parameters::Register<Parameters::VtkWritePressures>
        ("Include the phase pressures in the VTK output files");
    Parameters::Register<Parameters::VtkWriteDensities>
        ("Include the phase densities in the VTK output files");
    Parameters::Register<Parameters::VtkWriteSaturations>
        ("Include the phase saturations in the VTK output files");
    Parameters::Register<Parameters::VtkWriteMobilities>
        ("Include the phase mobilities in the VTK output files");
    Parameters::Register<Parameters::VtkWriteRelativePermeabilities>
        ("Include the phase relative permeabilities in the VTK output files");
    Parameters::Register<Parameters::VtkWriteViscosities>
        ("Include component phase viscosities in the VTK output files");
    Parameters::Register<Parameters::VtkWriteAverageMolarMasses>
        ("Include the average phase mass in the VTK output files");
    Parameters::Register<Parameters::VtkWritePorosity>
        ("Include the porosity in the VTK output files");
    Parameters::Register<Parameters::VtkWriteIntrinsicPermeabilities>
        ("Include the intrinsic permeability in the VTK output files");
    Parameters::Register<Parameters::VtkWriteFilterVelocities>
        ("Include in the filter velocities of the phases the VTK output files");
    Parameters::Register<Parameters::VtkWritePotentialGradients>
        ("Include the phase pressure potential gradients in the VTK output files");
}

void VtkMultiPhaseParams::read()
{
    extrusionFactorOutput_ = Parameters::Get<Parameters::VtkWriteExtrusionFactor>();
    pressureOutput_ = Parameters::Get<Parameters::VtkWritePressures>();
    densityOutput_ = Parameters::Get<Parameters::VtkWriteDensities>();
    saturationOutput_ = Parameters::Get<Parameters::VtkWriteSaturations>();
    mobilityOutput_ = Parameters::Get<Parameters::VtkWriteMobilities>();
    relativePermeabilityOutput_ = Parameters::Get<Parameters::VtkWriteRelativePermeabilities>();
    viscosityOutput_ = Parameters::Get<Parameters::VtkWriteViscosities>();
    averageMolarMassOutput_ = Parameters::Get<Parameters::VtkWriteAverageMolarMasses>();
    porosityOutput_ = Parameters::Get<Parameters::VtkWritePorosity>();
    intrinsicPermeabilityOutput_ = Parameters::Get<Parameters::VtkWriteIntrinsicPermeabilities>();
    velocityOutput_ = Parameters::Get<Parameters::VtkWriteFilterVelocities>();
    potentialGradientOutput_ = Parameters::Get<Parameters::VtkWritePotentialGradients>();
}

} // namespace Opm
