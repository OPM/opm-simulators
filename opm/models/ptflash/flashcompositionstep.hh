// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 SINTEF Digital

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
*/
/*!
 * \file
 *
 * \copydoc Opm::applyFlashCompositionStep
 */
#ifndef OPM_FLASH_COMPOSITION_STEP_HH
#define OPM_FLASH_COMPOSITION_STEP_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace Opm {

namespace detail {

//! Reserves the floor for every component and rescales the excess so that the fractions
//! sum to one. The input fractions must sum to one.
template <class Scalar, std::size_t numComponents>
void keepCompositionFloor(std::array<Scalar, numComponents>& z, const Scalar compositionFloor)
{
    Scalar excessSum = 0.0;
    for (auto& zc : z) {
        zc = std::max(zc - compositionFloor, Scalar{0});
        excessSum += zc;
    }
    // Fractions summing to one leave a positive excessSum.
    const Scalar excessTotal = 1 - numComponents * compositionFloor;
    for (auto& zc : z) {
        zc = compositionFloor + excessTotal * zc / excessSum;
    }
}

//! Largest change of the hydrocarbon share h or of an amount h z between two states.
template <class Scalar, std::size_t numComponents>
Scalar largestAmountChange(const Scalar share,
                           const std::array<Scalar, numComponents>& z,
                           const Scalar otherShare,
                           const std::array<Scalar, numComponents>& otherZ)
{
    Scalar change = std::abs(otherShare - share);
    for (std::size_t compIdx = 0; compIdx < numComponents; ++compIdx) {
        change = std::max(std::abs(otherShare * otherZ[compIdx] - share * z[compIdx]), change);
    }
    return change;
}

//! Replaces the old state (z, sw) by the point on the straight line in h and h z towards
//! (newZ, newShare) whose largest change from the old state is maxAmountChange. The line
//! starts from the old state moved onto the bounds, and that offset counts against the
//! limit. The bounds are linear in h and h z, so every point on the line keeps them.
template <class Scalar, std::size_t numComponents>
void shortenStep(std::array<Scalar, numComponents>& z,
                 Scalar& sw,
                 const std::array<Scalar, numComponents>& newZ,
                 const Scalar newShare,
                 const Scalar compositionFloor,
                 const Scalar hydrocarbonFloor,
                 const Scalar maxAmountChange)
{
    const Scalar share = std::max(1 - sw, hydrocarbonFloor);
    std::array<Scalar, numComponents> startZ = z;
    keepCompositionFloor(startZ, compositionFloor);
    const Scalar startShare = std::max(1 - std::clamp(sw, Scalar{0}, Scalar{1}),
                                       hydrocarbonFloor);

    // A start that uses up the limit is kept.
    const Scalar offset = largestAmountChange(share, z, startShare, startZ);
    const Scalar length = largestAmountChange(startShare, startZ, newShare, newZ);
    const Scalar scale = offset < maxAmountChange ? (maxAmountChange - offset) / length
                                                  : Scalar{0};

    const Scalar finalShare = startShare + scale * (newShare - startShare);
    for (std::size_t compIdx = 0; compIdx < numComponents; ++compIdx) {
        const Scalar startAmount = startShare * startZ[compIdx];
        const Scalar newAmount = newShare * newZ[compIdx];
        z[compIdx] = (startAmount + scale * (newAmount - startAmount)) / finalShare;
    }
    // Dividing by a small share amplifies roundoff in the amounts, so the floors and the
    // sum are restored.
    keepCompositionFloor(z, compositionFloor);
    // Sw follows from the share: interpolated separately, the two would disagree where
    // the floor on h applies.
    sw = 1 - finalShare;
}

} // namespace detail

/*!
 * \ingroup FlashModel
 *
 * \brief Applies a Newton step to the overall composition and the water saturation.
 *
 * The step is taken and limited in h and the component amounts h z, where
 * h = max(1 - Sw, hydrocarbonFloor) is the regularized hydrocarbon share of the pore
 * space. Component-storage derivatives scale with h, so Newton steps in the mole
 * fractions can become large in cells containing almost only water. Where h ends at its
 * floor, the hydrocarbon has vanished and z is kept, since the equations no longer
 * determine it.
 *
 * The returned fractions sum to one and are at or above compositionFloor,
 * and Sw stays within [0, 1]. Changes in h and each h z are limited to
 * maxAmountChange. If restoring an input state to these bounds requires a
 * larger change, restoring the bounds takes priority.
 *
 * \param z Overall mole fractions, summing to one, updated in place; zeros are allowed.
 * \param sw Water saturation, updated in place; zero without water.
 * \param dz Newton step of the mole fractions, summing to zero.
 * \param dSw Newton step of the water saturation; zero without water.
 * \param compositionFloor Lower bound of every mole fraction.
 * \param hydrocarbonFloor Lower bound of h.
 * \param maxAmountChange Largest change of h and of every amount.
 *
 * \pre Inputs are finite, 0 <= numComponents * compositionFloor < 1,
 *      0 < hydrocarbonFloor <= 1, and maxAmountChange > 0.
 */
template <class Scalar, std::size_t numComponents>
void applyFlashCompositionStep(std::array<Scalar, numComponents>& z,
                               Scalar& sw,
                               const std::array<Scalar, numComponents>& dz,
                               const Scalar dSw,
                               const Scalar compositionFloor,
                               const Scalar hydrocarbonFloor,
                               const Scalar maxAmountChange)
{
    const Scalar hydrocarbonShare = std::max(1 - sw, hydrocarbonFloor);
    const Scalar dHydrocarbonShare = -dSw;

    // One damping factor limits the linearized changes in h and h z. The changes after
    // clamping Sw and enforcing the floors are checked at the end.
    Scalar maxChange = std::abs(dHydrocarbonShare);
    for (std::size_t compIdx = 0; compIdx < numComponents; ++compIdx) {
        const Scalar dAmount = hydrocarbonShare * dz[compIdx] + z[compIdx] * dHydrocarbonShare;
        maxChange = std::max(std::abs(dAmount), maxChange);
    }
    const Scalar alpha = maxChange > maxAmountChange ? maxAmountChange / maxChange : 1.0;

    const Scalar newSw = std::clamp(sw + alpha * dSw, Scalar{0}, Scalar{1});
    const Scalar newHydrocarbonShare = std::max(1 - newSw, hydrocarbonFloor);

    // Back from amounts to fractions: z moves by hydrocarbonShare/newHydrocarbonShare of
    // the damped step. A vanished hydrocarbon keeps its composition: the equations see z
    // only through the floor of h then, so the step would mostly carry linear-solver error.
    const bool vanished = newHydrocarbonShare <= hydrocarbonFloor;
    const Scalar zStepScale = vanished ? Scalar{0}
                                       : alpha * hydrocarbonShare / newHydrocarbonShare;
    std::array<Scalar, numComponents> newZ{};
    for (std::size_t compIdx = 0; compIdx < numComponents; ++compIdx) {
        newZ[compIdx] = z[compIdx] + zStepScale * dz[compIdx];
    }
    detail::keepCompositionFloor(newZ, compositionFloor);

    // Clamping Sw and enforcing the floors can lengthen the step beyond the limit.
    const Scalar change = detail::largestAmountChange(hydrocarbonShare, z,
                                                      newHydrocarbonShare, newZ);
    if (change <= maxAmountChange) {
        z = newZ;
        sw = newSw;
    }
    else {
        detail::shortenStep(z, sw, newZ, newHydrocarbonShare,
                            compositionFloor, hydrocarbonFloor, maxAmountChange);
    }
}

} // namespace Opm

#endif
