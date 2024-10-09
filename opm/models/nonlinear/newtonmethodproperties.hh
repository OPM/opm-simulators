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
#ifndef EWOMS_NEWTON_METHOD_PROPERTIES_HH
#define EWOMS_NEWTON_METHOD_PROPERTIES_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Properties {

//! Specifies the type of the actual Newton method
template<class TypeTag, class MyTypeTag>
struct NewtonMethod { using type = UndefinedProperty; };

//! The class which linearizes the non-linear system of equations
template<class TypeTag, class MyTypeTag>
struct Linearizer { using type = UndefinedProperty; };

//! Specifies the type of the class which writes out the Newton convergence
template<class TypeTag, class MyTypeTag>
struct NewtonConvergenceWriter { using type = UndefinedProperty; };

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
template<class TypeTag, class MyTypeTag>
struct ConvergenceWriter { using type = UndefinedProperty; };

} // end namespace Opm::Properties

#endif
