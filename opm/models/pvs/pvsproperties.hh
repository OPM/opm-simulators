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
 * \ingroup PvsModel
 *
 * \brief Declares the properties required for the compositional
 *        multi-phase primary variable switching model.
 */
#ifndef EWOMS_PVS_PROPERTIES_HH
#define EWOMS_PVS_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/common/diffusionmodule.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkphasepresencemodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>

BEGIN_PROPERTIES

//! Specifies whether energy is considered as a conservation quantity or not
NEW_PROP_TAG(EnableEnergy);
//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! The verbosity of the model (0 -> do not print anything, 2 -> spam stdout a lot)
NEW_PROP_TAG(PvsVerbosity);
//! The basis value for the weight of the pressure primary variable
NEW_PROP_TAG(PvsPressureBaseWeight);
//! The basis value for the weight of the saturation primary variables
NEW_PROP_TAG(PvsSaturationsBaseWeight);
//! The basis value for the weight of the mole fraction primary variables
NEW_PROP_TAG(PvsMoleFractionsBaseWeight);

END_PROPERTIES

#endif
