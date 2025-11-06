// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#include <opm/models/io/vtktpsaparams.hpp>
#include <opm/models/utils/parametersystem.hpp>


namespace Opm {

/*!
* \brief Register runtime parameters
*/
void VtkTpsaParams::registerParameters()
{
    Parameters::Register<Parameters::VtkWriteDisplacement>
        ("Include displacement in VTK output files");
    Parameters::Register<Parameters::VtkWriteRotation>
        ("Include rotation in VTK output files");
    Parameters::Register<Parameters::VtkWriteSolidPressure>
        ("Include solid pressure in VTK output files");
}

void VtkTpsaParams::read()
{
    displacementOutput_ = Parameters::Get<Parameters::VtkWriteDisplacement>();
    rotationOutput_ = Parameters::Get<Parameters::VtkWriteRotation>();
    solidPressureOutput_ = Parameters::Get<Parameters::VtkWriteSolidPressure>();
}

}  // namespace Opm