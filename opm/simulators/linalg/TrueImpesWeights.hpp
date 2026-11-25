/*
  Copyright 2026 Equinor ASA

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

#include <type_traits>

namespace Opm::Amg {

enum class TrueImpesMethod { Storage, Analytic };

} // namespace Opm::Amg

namespace Opm::Amg::detail {

// The sequential residual is a linear combination of the
// mass balance residuals, with coefficients equal to (for
// water, oil, gas):
//    1/bw,
//    (1/bo - rs/bg)/(1-rs*rv)
//    (1/bg - rv/bo)/(1-rs*rv)
// These coefficients must be applied for both the residual and
// Jacobian.
template<class FluidSystem, class FluidState, class PrimaryVariables, class Weights>
OPM_HOST_DEVICE void analyticWeights(const FluidSystem& system, const FluidState& fs,
                                     const PrimaryVariables& priVars, Weights&& bweights)
{
    if (system.phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::waterPhaseIdx));
        bweights[activeCompIdx] = 1 / getValue(fs.invB(FluidSystem::waterPhaseIdx));
    }

    using Scalar = std::decay_t<decltype(getValue(fs.Rs()))>;
    double denominator = 1.0;
    Scalar rs = 0.0;
    Scalar rv = 0.0;
    // The dissolved gas and vaporised oil terms only exist when both
    // oil and gas are active.
    const bool oilAndGas = system.phaseIsActive(FluidSystem::oilPhaseIdx)
        && system.phaseIsActive(FluidSystem::gasPhaseIdx);
    if (oilAndGas) {
        rs = getValue(fs.Rs());
        rv = getValue(fs.Rv());
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
            rs = 0.0;
        }
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
            rv = 0.0;
        }
        denominator = 1.0 - static_cast<double>(rs) * rv;
    }
    // DenseAd divides by a scalar by multiplying by its reciprocal in the value type.
    const Scalar inverseDenominator = 1.0 / denominator;

    if (system.phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::oilPhaseIdx));
        auto oilWeight = 1 / getValue(fs.invB(FluidSystem::oilPhaseIdx));
        if (oilAndGas) {
            oilWeight -= rs / getValue(fs.invB(FluidSystem::gasPhaseIdx));
        }
        bweights[activeCompIdx] = oilWeight * inverseDenominator;
    }
    if (system.phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned activeCompIdx = system.canonicalToActiveCompIdx(
            system.solventComponentIndex(FluidSystem::gasPhaseIdx));
        auto gasWeight = 1 / getValue(fs.invB(FluidSystem::gasPhaseIdx));
        if (oilAndGas) {
            gasWeight -= rv / getValue(fs.invB(FluidSystem::oilPhaseIdx));
        }
        bweights[activeCompIdx] = gasWeight * inverseDenominator;
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
                                                        Scale storage_scale,
                                                        Index pressureVarIndex,
                                                        Matrix& block_transpose)
{
    const double pressure_scale = 50e5;

    // Build the transposed matrix directly to avoid separate transpose step
    for (Index ii = 0; ii < numEq; ++ii) {
        for (Index jj = 0; jj < numEq; ++jj) {
            auto& entry = trueImpesJacobianEntry<numEq>(block_transpose, jj, ii);
            entry = storage[ii].derivative(jj)/storage_scale;
            if (jj == pressureVarIndex) {
                entry *= pressure_scale;
            }
        }
    }
}

} // namespace Opm::Amg::detail

#endif // OPM_TRUE_IMPES_WEIGHTS_HPP
