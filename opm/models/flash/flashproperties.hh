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
 * \ingroup FlashModel
 *
 * \brief Declares the properties required by the compositional
 *        multi-phase model based on flash calculations.
 */
#ifndef EWOMS_FLASH_PROPERTIES_HH
#define EWOMS_FLASH_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>

BEGIN_PROPERTIES

//! Provides the thermodynamic relations
NEW_PROP_TAG(FluidSystem);
//! The type of the flash constraint solver
NEW_PROP_TAG(FlashSolver);
//! The maximum accepted error of the flash solver
NEW_PROP_TAG(FlashTolerance);

//! The thermal conduction law which ought to be used
NEW_PROP_TAG(ThermalConductionLaw);
//! The parameters of the thermal conduction law
NEW_PROP_TAG(ThermalConductionLawParams);

//! Specifies whether energy should be considered as a conservation quantity or not
NEW_PROP_TAG(EnableEnergy);
//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

END_PROPERTIES

#endif
