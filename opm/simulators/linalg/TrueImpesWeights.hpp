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

namespace Opm::Amg {

enum class TrueImpesMethod { Storage, Analytic };

} // namespace Opm::Amg

namespace Opm::Amg::detail {

// Volume coefficients of the water, oil and gas mass balances in the
// sequential pressure equation. Only the values of the PVT quantities enter.
template<class FluidSystem, class FluidState, class PrimaryVariables, class Weights>
OPM_HOST_DEVICE void analyticWeights(const FluidSystem& system, const FluidState& fs,
                                     const PrimaryVariables& priVars, Weights&& weights)
{
    const auto component = [&system](unsigned phase) {
        return system.canonicalToActiveCompIdx(system.solventComponentIndex(phase));
    };
    constexpr auto water = FluidSystem::waterPhaseIdx;
    constexpr auto oil = FluidSystem::oilPhaseIdx;
    constexpr auto gas = FluidSystem::gasPhaseIdx;

    if (system.phaseIsActive(water)) {
        weights[component(water)] = 1 / getValue(fs.invB(water));
    }

    const bool oilActive = system.phaseIsActive(oil);
    const bool gasActive = system.phaseIsActive(gas);
    if (oilActive && gasActive) {
        const auto bo = getValue(fs.invB(oil));
        const auto bg = getValue(fs.invB(gas));
        using Scalar = std::remove_cv_t<decltype(bo)>;
        const auto meaning = priVars.primaryVarsMeaningGas();
        const Scalar rs = meaning == PrimaryVariables::GasMeaning::Rv ? 0 : getValue(fs.Rs());
        const Scalar rv = meaning == PrimaryVariables::GasMeaning::Rs ? 0 : getValue(fs.Rv());
        // Match the value part of the host AD calculation for float and double.
        const Scalar scale = 1.0 / (1.0 - static_cast<double>(rs) * rv);
        weights[component(oil)] = (1 / bo - rs / bg) * scale;
        weights[component(gas)] = (1 / bg - rv / bo) * scale;
    }
    else if (oilActive) {
        weights[component(oil)] = 1 / getValue(fs.invB(oil));
    }
    else if (gasActive) {
        weights[component(gas)] = 1 / getValue(fs.invB(gas));
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
