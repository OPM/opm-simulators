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
 * \ingroup NcpModel
 *
 * \brief Declares the properties required for the NCP compositional
 *        multi-phase model.
 */
#ifndef EWOMS_NCP_PROPERTIES_HH
#define EWOMS_NCP_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>

namespace Opm::Properties {

//! The unmodified weight for the pressure primary variable
template<class TypeTag, class MyTypeTag>
struct NcpPressureBaseWeight { using type = UndefinedProperty; };
//! The weight for the saturation primary variables
template<class TypeTag, class MyTypeTag>
struct NcpSaturationsBaseWeight { using type = UndefinedProperty; };
//! The unmodified weight for the fugacity primary variables
template<class TypeTag, class MyTypeTag>
struct NcpFugacitiesBaseWeight { using type = UndefinedProperty; };

//! The themodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
template<class TypeTag, class MyTypeTag>
struct NcpCompositionFromFugacitiesSolver { using type = UndefinedProperty; };

} // namespace Opm::Properties

#endif
