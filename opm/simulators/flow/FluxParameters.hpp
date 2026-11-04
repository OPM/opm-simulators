/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_FLUX_PARAMETERS_HPP
#define OPM_FLUX_PARAMETERS_HPP

namespace Opm::Parameters {

// In a USEFLUX (sector) run, trigger the ACTIONX whose conditions read
// quantities the sector cannot compute itself when the parent run triggered
// them, as recorded in the FLUX file, instead of evaluating their conditions.
struct FluxReplayActions { static constexpr bool value = true; };

// In a USEFLUX (sector) run, evaluate FPR as FPRFLUX in UDQ definitions and
// ACTIONX conditions.
struct FluxFprAsFprflux { static constexpr bool value = false; };

} // namespace Opm::Parameters

#endif // OPM_FLUX_PARAMETERS_HPP
