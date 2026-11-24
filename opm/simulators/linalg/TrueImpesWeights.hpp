/*
  Copyright 2026 SINTEF Digital

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

#ifndef OPM_TRUE_IMPES_WEIGHTS_HPP
#define OPM_TRUE_IMPES_WEIGHTS_HPP

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace Opm::Amg::detail {

template<class FluidSystem, class FluidState, class PrimaryVariables, class Weights>
OPM_HOST_DEVICE void analyticWeights(const FluidSystem& system, const FluidState& fs,
                                     const PrimaryVariables& priVars, Weights&& weights)
{
    if (system.phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::waterPhaseIdx));
        weights[activeCompIdx]
            = 1.0 / getValue(fs.invB(FluidSystem::waterPhaseIdx));
    }

    double denominator = 1.0;
    double rs = getValue(fs.Rs());
    double rv = getValue(fs.Rv());

    if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
        rs = 0.0;
    }
    if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
        rv = 0.0;
    }
    if (system.phaseIsActive(FluidSystem::oilPhaseIdx)
        && system.phaseIsActive(FluidSystem::gasPhaseIdx)) {
        denominator = 1 - rs * rv;
    }

    if (system.phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::oilPhaseIdx));
        weights[activeCompIdx] = (1.0 / getValue(fs.invB(FluidSystem::oilPhaseIdx))
                              - rs / getValue(fs.invB(FluidSystem::gasPhaseIdx))) / denominator;
    }
    if (system.phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::gasPhaseIdx));
        weights[activeCompIdx] = (1.0 / getValue(fs.invB(FluidSystem::gasPhaseIdx))
                              - rv / getValue(fs.invB(FluidSystem::oilPhaseIdx))) / denominator;
    }
}

template <int numEq, class Matrix, class Index>
OPM_HOST_DEVICE inline auto& trueImpesJacobianEntry(Matrix& block, Index pv, Index eq)
{
    if constexpr (std::is_array_v<Matrix>) {
        return block[pv * numEq + eq];
    } else {
        return block[pv][eq];
    }
}

// Each row of the transposed storage Jacobian corresponds to a primary variable.
// Use the native matrix on the host and a flat row-major block in the kernel.
template <int numEq, class Storage, class Scale, class Index, class Matrix>
OPM_HOST_DEVICE inline void formTrueImpesStorageJacobian(const Storage& storage,
                                                        Scale storageScale,
                                                        Index pressureVariable,
                                                        Matrix& block)
{
    constexpr double pressureScale = 50e5;
    for (Index eq = 0; eq < numEq; ++eq) {
        for (Index pv = 0; pv < numEq; ++pv) {
            auto& entry = trueImpesJacobianEntry<numEq>(block, pv, eq);
            entry = storage[eq].derivative(pv) / storageScale;
            if (pv == pressureVariable) {
                entry *= pressureScale;
            }
        }
    }
}

// The result may alias the solution or point directly into the output vector.
template <int numEq, class Scalar>
OPM_HOST_DEVICE inline bool normalizeTrueImpesWeights(const Scalar* solution, Scalar* weights)
{
#if OPM_IS_INSIDE_DEVICE_FUNCTION
    using ::fabs;
    using ::isfinite;
    using ::max;
#else
    using std::fabs;
    using std::isfinite;
    using std::max;
#endif
    Scalar maxAbs = 0;
    for (int eq = 0; eq < numEq; ++eq) {
        const Scalar weight = solution[eq];
        if (!isfinite(weight)) {
            return false;
        }
        maxAbs = max(maxAbs, fabs(weight));
    }
    if (maxAbs == 0) {
        return false;
    }
    for (int eq = 0; eq < numEq; ++eq) {
        weights[eq] = solution[eq] / maxAbs;
    }
    return true;
}

} // namespace Opm::Amg::detail

#endif // OPM_TRUE_IMPES_WEIGHTS_HPP
