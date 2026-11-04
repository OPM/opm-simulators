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

#ifndef OPM_INJECTED_PHASE_HPP_INCLUDED
#define OPM_INJECTED_PHASE_HPP_INCLUDED

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>

#include <optional>

namespace Opm {

/// Canonical phase index (IndexTraits::waterPhaseIdx, oilPhaseIdx or
/// gasPhaseIdx) of the single phase injected by an injector of the given
/// type.  Returns nullopt for MULTI, which has no single injected phase.
template <class IndexTraits>
std::optional<unsigned> injectedCanonicalPhaseIdx(const InjectorType type)
{
    switch (type) {
    case InjectorType::WATER: return IndexTraits::waterPhaseIdx;
    case InjectorType::OIL:   return IndexTraits::oilPhaseIdx;
    case InjectorType::GAS:   return IndexTraits::gasPhaseIdx;
    case InjectorType::MULTI: return std::nullopt;
    }
    return std::nullopt;
}

/// Active phase index of the phase injected by an injector of the given
/// type, i.e. the position of that phase in per-phase well and connection
/// rate arrays such as PerfData::phase_rates.  Returns nullopt for MULTI
/// injectors and if the injected phase is not active in the run.
template <class IndexTraits>
std::optional<int> injectedActivePhaseIdx(const PhaseUsageInfo<IndexTraits>& pu,
                                          const InjectorType type)
{
    const auto canonical = injectedCanonicalPhaseIdx<IndexTraits>(type);
    if (!canonical.has_value() || !pu.phaseIsActive(*canonical)) {
        return std::nullopt;
    }
    return pu.canonicalToActivePhaseIdx(*canonical);
}

} // namespace Opm

#endif // OPM_INJECTED_PHASE_HPP_INCLUDED
