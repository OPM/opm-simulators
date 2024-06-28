// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center
  
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
 *
 * \copydoc Opm::DamarisWriter
 */
#ifndef OPM_DAMARIS_PARAMETERS_HPP
#define OPM_DAMARIS_PARAMETERS_HPP

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutput { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisOutputHdfCollective { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisSaveMeshToHdf { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisSaveToHdf { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisPythonScript { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisPythonParaviewScript { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisSimName { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisDedicatedCores { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisDedicatedNodes { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisSharedMemoryName { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisSharedMemorySizeBytes { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisLogLevel { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisDaskFile { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DamarisLimitVariables { using type = Properties::UndefinedProperty; };

} // namespace Opm::Parameters

#endif // OPM_DAMARIS_PARAMETERS_HPP
