/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_WELL_PRIMARY_VARIABLES_IMPL_HPP
#define OPM_COMP_WELL_PRIMARY_VARIABLES_IMPL_HPP

#include <flowexperimental/comp/wells/SingleCompWellState.hpp>

namespace Opm {

template <typename FluidSystem, typename Indices>
void
CompWellPrimaryVariables<FluidSystem, Indices>::
update(const SingleWellState& well_state)
{
    value_[QTotal] = well_state.get_total_surface_rate();
    // the mole fractions of the first n-1 component
    std::vector<Scalar> mole_fractions(FluidSystem::numComponents, 0.);
    Scalar sum_mole_fraction = 0.;
    for (int i = 0; i < FluidSystem::numComponents; ++i) {
        mole_fractions[i] = std::max(well_state.total_molar_fractions[i], 1.e-10);
        sum_mole_fraction += mole_fractions[i];
    }
    assert(sum_mole_fraction != 0.);
    for (int i = 0; i < numWellConservationEq - 1; ++i) {
        value_[i + 1] = mole_fractions[i] / sum_mole_fraction;
    }
    value_[Bhp] = well_state.bhp;

    temperature_ = well_state.temperature;

    updateEvaluation();
}

template <typename FluidSystem, typename Indices>
void
CompWellPrimaryVariables<FluidSystem, Indices>::
updateEvaluation()
{
    for (std::size_t idx = 0; idx < numWellEq; ++idx) {
        evaluation_[idx] = 0.;
        evaluation_[idx].setValue(value_[idx]);
        evaluation_[idx].setDerivative(idx + numResEq, 1.);
    }
}

template <typename FluidSystem, typename Indices>
typename CompWellPrimaryVariables<FluidSystem, Indices>::EvalWell
CompWellPrimaryVariables<FluidSystem, Indices>::
getBhp() const
{
    return evaluation_[Bhp];
}

template <typename FluidSystem, typename Indices>
typename CompWellPrimaryVariables<FluidSystem, Indices>::EvalWell
CompWellPrimaryVariables<FluidSystem, Indices>::
getTotalRate() const
{
    return evaluation_[QTotal];
}

template <typename FluidSystem, typename Indices>
typename CompWellPrimaryVariables<FluidSystem, Indices>::EvalWell
CompWellPrimaryVariables<FluidSystem, Indices>::
extendEval(const Eval& in)
{
    EvalWell out = 0.0;
    out.setValue(in.value());
    for(int eq_idx = 0; eq_idx < Indices::numEq;++eq_idx) {
        out.setDerivative(eq_idx, in.derivative(eq_idx));
    }
    return out;
}

template <typename FluidSystem, typename Indices>
typename CompWellPrimaryVariables<FluidSystem, Indices>::Eval
CompWellPrimaryVariables<FluidSystem, Indices>::
restrictEval(const EvalWell& in)
{
    Eval out = 0.0;
    out.setValue(in.value());
    for(int eq_idx = 0; eq_idx < Indices::numEq;++eq_idx) {
        out.setDerivative(eq_idx, in.derivative(eq_idx));
    }
    return out;
}

template <typename FluidSystem, typename Indices>
void
CompWellPrimaryVariables<FluidSystem, Indices>::
updateNewton(const BVectorWell& dwells)
{
    constexpr Scalar damping = 1.0;

    for (unsigned i = 0; i < value_.size(); ++i) {
        value_[i] -= damping * dwells[0][i];
    }
    // TODO: more general indices here
    value_[1] = std::clamp(value_[1], 1.e-10, 1.);
    value_[2] = std::clamp(value_[2], 1.e-10, 1.);
    std::vector<Scalar> mole_fractions(FluidSystem::numComponents, 0.);
    Scalar sum_mole_fraction = 0.;
    for (int i = 0; i < FluidSystem::numComponents-1; ++i) {
        mole_fractions[i] = std::max(value_[i + 1], 1.e-10);
        sum_mole_fraction += mole_fractions[i];
    }
    mole_fractions[FluidSystem::numComponents - 1] = std::max(1.0 - sum_mole_fraction, 1.e-10);
    sum_mole_fraction += mole_fractions[FluidSystem::numComponents - 1];
    assert(sum_mole_fraction != 0.);
    value_[1] = mole_fractions[0] / sum_mole_fraction;
    value_[2] = mole_fractions[1] / sum_mole_fraction;

    updateEvaluation();
}

template <typename FluidSystem, typename Indices>
template <typename T>
T
CompWellPrimaryVariables<FluidSystem, Indices>::
getValue_(int index) const
{
    static_assert(std::is_same_v<T, Scalar> || std::is_same_v<T, EvalWell>, "Unsupported type in CompWellPrimaryVariables::getValue_");

    if constexpr (std::is_same_v<T, Scalar>) {
        return value_[index];
    } else {
        return evaluation_[index];
    }
}

template <typename FluidSystem, typename Indices>
template <typename T>
typename CompWellPrimaryVariables<FluidSystem, Indices>::template FluidState<T>
CompWellPrimaryVariables<FluidSystem, Indices>::
toFluidState() const
{
    static_assert(std::is_same_v<T, Scalar> || std::is_same_v<T, EvalWell>, "Unsupported type in CompWellPrimaryVariables::toFluidState");

    CompositionalFluidState<T, FluidSystem> fluid_state;
    const auto& pressure = getValue_<T>(Bhp);
    std::array<T, FluidSystem::numComponents> total_molar_fractions;
    T sum = 0.;
    for (int i = 0; i < FluidSystem::numComponents - 1; ++i) {
        total_molar_fractions[i] = getValue_<T>(i + 1);
        sum += total_molar_fractions[i];
    }
    total_molar_fractions[FluidSystem::numComponents - 1] = 1.0 - sum;

    for (int i = 0; i < FluidSystem::numComponents; ++i) {
        if constexpr (std::is_same_v<T, EvalWell>) {
            total_molar_fractions[i].setValue(std::max(getValue(total_molar_fractions[i]), 1.e-10));
        } else { // Scalar
            total_molar_fractions[i] = std::max(total_molar_fractions[i], 1.e-10);
        }
        fluid_state.setMoleFraction(i, total_molar_fractions[i]);
    }

    fluid_state.setPressure(FluidSystem::oilPhaseIdx, pressure);
    fluid_state.setPressure(FluidSystem::gasPhaseIdx, pressure);

    fluid_state.setTemperature(temperature_);

    for (int i = 0; i < FluidSystem::numComponents; ++i) {
        fluid_state.setKvalue(i, fluid_state.wilsonK_(i));
    }

    fluid_state.setLvalue(-1.);

    return fluid_state;
}

} // end of namespace Opm

#endif // OPM_COMP_WELL_PRIMARY_VARIABLES_IMPL_HPP
