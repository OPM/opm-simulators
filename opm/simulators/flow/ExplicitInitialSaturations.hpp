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

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

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
inline ExplicitInitialSaturations
completeExplicitSaturations(const std::array<std::optional<double>, 3>& supplied,
                            const std::array<bool, 3>& active)
{
    // Allow closure errors from rounded input, including float32 restart data.
    constexpr double inputTolerance = 1.0e-6;
    constexpr double normalizationTolerance = 0.05;
    constexpr std::array names {"SWAT", "SOIL", "SGAS"};
    ExplicitInitialSaturations result;
    auto& saturations = result.saturations;
    std::size_t activeCount = 0;
    std::size_t suppliedCount = 0;
    std::size_t missingPhase = 0;
    double sum = 0.0;

    for (std::size_t phase = 0; phase < active.size(); ++phase) {
        if (!active[phase]) {
            continue;
        }
        ++activeCount;
        if (!supplied[phase]) {
            missingPhase = phase;
            continue;
        }
        const double saturation = *supplied[phase];
        if (!std::isfinite(saturation) || saturation < 0.0) {
            throw std::invalid_argument(
                fmt::format("{} must be finite and nonnegative", names[phase]));
        }
        saturations[phase] = saturation;
        sum += saturation;
        ++suppliedCount;
    }

    if (activeCount == 0) {
        throw std::invalid_argument("Explicit initialization requires an active phase");
    }
    if (suppliedCount + 1 < activeCount) {
        throw std::invalid_argument("Explicit initialization needs saturations for all but one "
                                    "active phase (SWAT, SOIL, SGAS)");
    }
    if (suppliedCount < activeCount) {
        if (!std::isfinite(sum) || sum > 1.0 + inputTolerance) {
            throw std::invalid_argument(
                fmt::format("The supplied saturations sum to more than 1 + {}", inputTolerance));
        }
        if (sum > 1.0) {
            for (double& saturation : saturations) {
                saturation /= sum;
            }
        } else {
            saturations[missingPhase] = 1.0 - sum;
        }
        return result;
    }

    // All active phases are supplied, so normalize the input instead of
    // assigning the remaining pore volume to an omitted phase.
    const double tolerance = activeCount > 1 ? normalizationTolerance : inputTolerance;
    if (!std::isfinite(sum) || sum < 1.0 - tolerance || sum > 1.0 + tolerance) {
        throw std::invalid_argument(fmt::format(
            "{} must sum to 1 within a tolerance of {}",
            activeCount == 3 ? "SWAT, SOIL and SGAS" : "Saturations of the active phases",
            tolerance));
    }
    const double deviation = std::abs(sum - 1.0);
    if (activeCount > 1 && deviation > inputTolerance) {
        result.normalizationDeviation = deviation;
    }
    for (double& saturation : saturations) {
        saturation /= sum;
    }
    return result;
}

} // namespace Opm

#endif // OPM_EXPLICIT_INITIAL_SATURATIONS_HPP
