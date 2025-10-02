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

namespace Opm::Parameters {

struct EnableDamarisOutput { static constexpr bool value = false; };
struct DamarisOutputHdfCollective { static constexpr bool value = true; };
struct DamarisSaveMeshToHdf { static constexpr bool value = false; };
struct DamarisSaveToHdf { static constexpr bool value = true; };
struct DamarisPythonScript { static constexpr auto value = ""; };
struct DamarisPythonParaviewScript { static constexpr auto value = ""; };
struct DamarisSimName { static constexpr auto value = ""; };
struct DamarisDedicatedCores { static constexpr int value = 1; };
struct DamarisDedicatedNodes { static constexpr int value = 0; };
struct DamarisSharedMemoryName { static constexpr auto value = "" ; };
struct DamarisSharedMemorySizeBytes { static constexpr long value = 536870912; }; // 512 MB
struct DamarisLogLevel { static constexpr auto value = "info"; };
struct DamarisDaskFile { static constexpr auto value = ""; };
struct DamarisLimitVariables  { static constexpr auto value = ""; };

} // namespace Opm::Parameters

namespace Opm {

//! \brief Register damaris runtime parameters.
void registerDamarisParameters();

}

#endif // OPM_DAMARIS_PARAMETERS_HPP
