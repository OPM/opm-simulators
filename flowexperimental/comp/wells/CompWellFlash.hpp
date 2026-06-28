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

#ifndef OPM_COMP_WELL_FLASH_HPP
#define OPM_COMP_WELL_FLASH_HPP

#include <opm/material/Constants.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <type_traits>

namespace Opm {

/// Flash the wellbore fluid at its current (pressure, temperature, overall
/// composition) and fill in the per-phase saturations, compressibility factors
/// and densities for the two-phase (oil/gas) system.
///
/// \tparam FluidSystem  compositional fluid system (oil/gas, two phases)
/// \tparam T            value type of the fluid state: the fluid system's Scalar
///                      for the explicit/old-time-level quantities, or an
///                      Evaluation for the AD path. The dispatch between the
///                      scalar and the AD flash solver is done on this type.
/// \param flash_tolerance  convergence tolerance handed to the flash solver. The
///                      default reproduces the value used in production; the unit
///                      test tightens it so the finite-difference comparison is
///                      not dominated by flash-convergence noise.
///
/// This is a free-function extraction of the wellbore flash that previously
/// lived inline in CompWell, so it can be unit tested (e.g. by a finite
/// difference check of the AD derivatives) without instantiating a full
/// Simulator. The behaviour is intentionally identical to the original.
template <typename FluidSystem, typename T>
void flashWellboreFluidState(CompositionalFluidState<T, FluidSystem>& fluid_state,
                             const typename FluidSystem::Scalar flash_tolerance = 1.e-6)
{
    using Scalar = typename FluidSystem::Scalar;
    using EOSType = CompositionalConfig::EOSType;

    bool single_phase = false;
    if constexpr (std::is_same_v<T, Scalar>) {
        single_phase = PTFlash<Scalar, FluidSystem>::flash_solve_scalar_(fluid_state, "ssi", flash_tolerance, EOSType::PR);
    } else { // Evaluation
        single_phase = PTFlash<Scalar, FluidSystem>::solve(fluid_state, "ssi", flash_tolerance, EOSType::PR);
    }

    constexpr Scalar R = Constants<Scalar>::R;
    typename FluidSystem::template ParameterCache<T> param_cache {EOSType::PR};
    param_cache.updatePhase(fluid_state, FluidSystem::oilPhaseIdx);
    const auto Z_L = (param_cache.molarVolume(FluidSystem::oilPhaseIdx) * fluid_state.pressure(FluidSystem::oilPhaseIdx)) /
                     (R * fluid_state.temperature(FluidSystem::oilPhaseIdx));
    param_cache.updatePhase(fluid_state, FluidSystem::gasPhaseIdx);
    const auto Z_V = (param_cache.molarVolume(FluidSystem::gasPhaseIdx) * fluid_state.pressure(FluidSystem::gasPhaseIdx)) /
                     (R * fluid_state.temperature(FluidSystem::gasPhaseIdx));

    auto L = fluid_state.L();
    if (single_phase) {
        // we check whether the phase label is correct
        if (L > 0.9 && Z_L > 0.8) { // marked as liquid phase while compress factor shows it is gas
            L = 0.;
            fluid_state.setLvalue(L);
        } else if (L < 0.1 && Z_V < 0.5) { // marked as gas phase while compress factor shows it is liquid
            L = 1.;
            fluid_state.setLvalue(L);
        }
    }
    T So = Opm::max((L * Z_L / (L * Z_L + (1 - L) * Z_V)), 0.0);
    T Sg = Opm::max(1 - So, 0.0);
    T sumS = So + Sg;
    So /= sumS;
    Sg /= sumS;

    fluid_state.setSaturation(FluidSystem::oilPhaseIdx, So);
    fluid_state.setSaturation(FluidSystem::gasPhaseIdx, Sg);

    fluid_state.setCompressFactor(FluidSystem::oilPhaseIdx, Z_L);
    fluid_state.setCompressFactor(FluidSystem::gasPhaseIdx, Z_V);

    fluid_state.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fluid_state, param_cache, FluidSystem::oilPhaseIdx));
    fluid_state.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fluid_state, param_cache, FluidSystem::gasPhaseIdx));
}

/// Mass [kg] of each component held in the wellbore control volume, given an
/// already-flashed two-phase fluid state and the wellbore volume.
///
/// mass_c = (y^o_c * rho_o * S_o + y^g_c * rho_g * S_g) * V
///
/// where y^p_c is the mass fraction of component c in phase p. The value type T
/// carries through, so when T is an Evaluation the returned masses carry the AD
/// derivatives with respect to the wellbore primary variables.
template <typename FluidSystem, typename T, typename Scalar>
std::array<T, FluidSystem::numComponents>
wellboreComponentMasses(const CompositionalFluidState<T, FluidSystem>& fluid_state,
                        const Scalar wellbore_volume)
{
    std::array<T, FluidSystem::numComponents> component_masses;

    const auto& so = fluid_state.saturation(FluidSystem::oilPhaseIdx);
    const auto& sg = fluid_state.saturation(FluidSystem::gasPhaseIdx);
    const auto& density_oil = fluid_state.density(FluidSystem::oilPhaseIdx);
    const auto& density_gas = fluid_state.density(FluidSystem::gasPhaseIdx);

    for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
        const auto oil_mass_fraction = fluid_state.massFraction(FluidSystem::oilPhaseIdx, comp_idx);
        const auto gas_mass_fraction = fluid_state.massFraction(FluidSystem::gasPhaseIdx, comp_idx);
        component_masses[comp_idx] = (oil_mass_fraction * density_oil * so +
                                      gas_mass_fraction * density_gas * sg) * wellbore_volume;
    }

    return component_masses;
}

} // end of namespace Opm

#endif // OPM_COMP_WELL_FLASH_HPP
