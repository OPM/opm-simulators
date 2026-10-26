/*
  Copyright 2026, SINTEF Digital

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

#ifndef OPM_EXPLICIT_INITIAL_SATURATIONS_HPP
#define OPM_EXPLICIT_INITIAL_SATURATIONS_HPP

#include <array>
#include <optional>

namespace Opm
{

struct ExplicitInitialSaturations {
    std::array<double, 3> saturations {};
    /// Absolute input closure error for a nontrivial normalization.
    double normalizationDeviation = 0.0;
};

/// Complete explicit saturations in water, oil, gas order. Inactive phases
/// have zero saturation and their input is ignored. At most one active phase
/// may be omitted; its saturation is the remaining pore volume. When all
/// active phases are supplied in a two- or three-phase model, normalize totals
/// within 0.05 of unity. With one phase omitted, a small excess within the
/// input tolerance is normalized with zero remainder; a positive remainder
/// is retained, however small.
ExplicitInitialSaturations
completeExplicitSaturations(const std::array<std::optional<double>, 3>& supplied,
                            const std::array<bool, 3>& active);

} // namespace Opm

#endif // OPM_EXPLICIT_INITIAL_SATURATIONS_HPP
