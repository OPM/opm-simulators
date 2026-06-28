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
/*!
 * \file
 *
 * \brief Finite-difference check of the automatic-differentiation derivatives
 *        that flow through the compositional wellbore flash.
 *
 * The compositional well (flowexperimental/comp/wells) assembles its Jacobian
 * from Evaluation (dense-AD) quantities that depend on the wellbore primary
 * variables (bottom-hole pressure and the overall mole fractions) through a
 * PT flash. A comment in CompWell_impl.hpp flagged the mass-fraction
 * derivatives as suspicious. This test pins those derivatives down: it builds a
 * wellbore fluid state with the pressure and composition as AD variables, runs
 * the extracted flashWellboreFluidState() + wellboreComponentMasses() helpers,
 * and compares the resulting analytical derivatives of
 *   - the per-component masses in the wellbore,
 *   - the per-component mass fractions, and
 *   - the wellbore fluid density,
 * against central finite differences computed by perturbing the scalar primary
 * variables and re-flashing.
 */
#include "config.h"

#define BOOST_TEST_MODULE CompWellJacobian
#include <boost/test/unit_test.hpp>

#include <flowexperimental/comp/wells/CompWellFlash.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <type_traits>

namespace {

using Scalar = double;

// Two-phase (oil/gas) three-component PT-flash fluid system. The generic fluid
// system registers its component data at runtime, see
// registerFluidSystemComponents() below.
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, /*enableWater=*/false>;

constexpr int numComponents = FluidSystem::numComponents; // 3

// Register the fixed CO2/Methane/Decane composition with the generic fluid
// system. The component data is shared static state, so this must run before any
// flash.
void registerFluidSystemComponents()
{
    using CO2 = Opm::SimpleCO2<Scalar>;
    using C1  = Opm::C1<Scalar>;
    using C10 = Opm::C10<Scalar>;
    using CompParam = FluidSystem::ComponentParam;

    FluidSystem::init();
    FluidSystem::addComponent(CompParam{CO2::name(), CO2::molarMass(), CO2::criticalTemperature(),
                                        CO2::criticalPressure(), CO2::criticalVolume(), CO2::acentricFactor()});
    FluidSystem::addComponent(CompParam{C1::name(), C1::molarMass(), C1::criticalTemperature(),
                                        C1::criticalPressure(), C1::criticalVolume(), C1::acentricFactor()});
    FluidSystem::addComponent(CompParam{C10::name(), C10::molarMass(), C10::criticalTemperature(),
                                        C10::criticalPressure(), C10::criticalVolume(), C10::acentricFactor()});
}

// The wellbore primary variables that the component masses depend on are the
// bottom-hole pressure and the first (numComponents - 1) overall mole
// fractions; the last mole fraction is the dependent 1 - sum. (The total rate
// and the constant temperature do not enter the masses, so they are left out.)
// Derivative slot layout: 0 -> pressure, 1 -> z0, 2 -> z1.
constexpr int numDeriv = numComponents; // pressure + (z0, z1)
constexpr int pIdx = 0;
constexpr int z0Idx = 1;
constexpr int z1Idx = 2;

using Evaluation = Opm::DenseAd::Evaluation<Scalar, numDeriv>;

// Quantities derived from a flashed wellbore fluid state, carrying value type T.
template <typename T>
struct WellboreQuantities
{
    std::array<T, numComponents> component_masses{};
    std::array<T, numComponents> mass_fractions{};
    T fluid_density{};
};

// Build the wellbore fluid state the same way CompWellPrimaryVariables::toFluidState
// does: overall mole fractions (clamped to >= 1e-10), pressure on both phases,
// constant temperature, Wilson initial K-values and an unset L flag.
template <typename T>
Opm::CompositionalFluidState<T, FluidSystem>
makeWellboreFluidState(const T& pressure,
                       const std::array<T, numComponents>& z,
                       const Scalar temperature)
{
    Opm::CompositionalFluidState<T, FluidSystem> fs;

    for (int i = 0; i < numComponents; ++i) {
        T zi = z[i];
        if constexpr (std::is_same_v<T, Scalar>) {
            zi = std::max(zi, Scalar(1.e-10));
        } else {
            zi.setValue(std::max(Opm::getValue(z[i]), Scalar(1.e-10)));
        }
        fs.setMoleFraction(i, zi);
    }

    fs.setPressure(FluidSystem::oilPhaseIdx, pressure);
    fs.setPressure(FluidSystem::gasPhaseIdx, pressure);
    fs.setTemperature(T(temperature));

    for (int i = 0; i < numComponents; ++i) {
        fs.setKvalue(i, fs.wilsonK_(i));
    }
    fs.setLvalue(T(-1.));

    return fs;
}

// Flash the wellbore fluid and assemble the quantities under test.
template <typename T>
WellboreQuantities<T>
computeWellboreQuantities(const T& pressure,
                          const std::array<T, numComponents>& z,
                          const Scalar temperature,
                          const Scalar wellbore_volume,
                          const Scalar flash_tolerance)
{
    auto fs = makeWellboreFluidState<T>(pressure, z, temperature);
    Opm::flashWellboreFluidState(fs, flash_tolerance);

    WellboreQuantities<T> q;
    q.component_masses = Opm::wellboreComponentMasses(fs, wellbore_volume);

    T total_mass = 0.;
    for (int c = 0; c < numComponents; ++c) {
        total_mass += q.component_masses[c];
    }
    for (int c = 0; c < numComponents; ++c) {
        q.mass_fractions[c] = q.component_masses[c] / total_mass;
    }

    const auto& so = fs.saturation(FluidSystem::oilPhaseIdx);
    const auto& sg = fs.saturation(FluidSystem::gasPhaseIdx);
    q.fluid_density = fs.density(FluidSystem::oilPhaseIdx) * so
                    + fs.density(FluidSystem::gasPhaseIdx) * sg;

    return q;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(WellboreFlashDerivatives)
{
    registerFluidSystemComponents();

    // A composition/pressure that sits comfortably inside the two-phase region,
    // so the flash is smooth and the central differences are meaningful.
    const Scalar temperature = 300.0;            // K
    const Scalar wellbore_volume = 21.6e-3;      // m^3 (matches CompWell)
    const Scalar p0 = 10.0e5;                    // Pa
    const Scalar z0_0 = 0.5;
    const Scalar z1_0 = 0.3;

    // Tighten the flash tolerance well below the perturbation-induced signal so
    // the finite differences are not swamped by flash-convergence noise.
    const Scalar flash_tolerance = 1.e-8;

    // --- Analytical (AD) quantities at the base point ---------------------
    Evaluation P  = Evaluation::createVariable(p0, pIdx);
    std::array<Evaluation, numComponents> z;
    z[0] = Evaluation::createVariable(z0_0, z0Idx);
    z[1] = Evaluation::createVariable(z1_0, z1Idx);
    z[2] = 1.0 - z[0] - z[1];

    const auto qad = computeWellboreQuantities<Evaluation>(P, z, temperature,
                                                           wellbore_volume, flash_tolerance);

    // Sanity guard: the base point must be genuinely two-phase, otherwise the
    // saturation-dependent derivatives vanish and the test becomes vacuous. A
    // substantial composition sensitivity of the fluid density confirms we are
    // exercising the two-phase flash coupling.
    BOOST_TEST_MESSAGE("base-point fluid density = " << qad.fluid_density.value());
    BOOST_REQUIRE_GT(qad.fluid_density.value(), 0.0);
    BOOST_REQUIRE_GT(std::abs(qad.fluid_density.derivative(z0Idx)), 1.0);

    // Characteristic magnitudes, used to build the absolute tolerance floor so
    // that near-zero derivatives are not held to an unreachable relative bound.
    Scalar mass_scale = 0.0;
    for (int c = 0; c < numComponents; ++c) {
        mass_scale = std::max(mass_scale, std::abs(qad.component_masses[c].value()));
    }
    const Scalar rho_scale = std::abs(qad.fluid_density.value());

    // --- Central finite differences ---------------------------------------
    // Relative step per primary variable (pressure is O(1e6), mole fractions O(1)).
    const std::array<Scalar, numDeriv> base_value{p0, z0_0, z1_0};
    const std::array<Scalar, numDeriv> eps{p0 * 1.e-5, 1.e-4, 1.e-4};

    // The AD and central-difference derivatives are observed to agree to ~1e-7
    // (relative) for the significant derivatives, so 1e-3 is a strong guard with
    // a comfortable margin against platform-dependent flash-convergence noise.
    // The absolute floor lets the genuinely-zero derivatives (e.g. the overall
    // mass fractions are pressure-independent) pass without a relative bound.
    const Scalar rel_tol = 1.e-3;
    const Scalar abs_floor = 1.e-3; // times the quantity scale

    auto quantitiesAt = [&](Scalar p, Scalar zz0, Scalar zz1) {
        const std::array<Scalar, numComponents> zs{zz0, zz1, 1.0 - zz0 - zz1};
        return computeWellboreQuantities<Scalar>(p, zs, temperature,
                                                 wellbore_volume, flash_tolerance);
    };

    auto checkDeriv = [&](Scalar ad_deriv, Scalar fd_deriv, Scalar scale,
                          const std::string& what) {
        const Scalar tol = rel_tol * std::abs(ad_deriv) + abs_floor * scale;
        // Diagnostic, silent at the default log level (use --log_level=message).
        BOOST_TEST_MESSAGE(what << ": AD=" << ad_deriv << " FD=" << fd_deriv
                                << " |diff|=" << std::abs(ad_deriv - fd_deriv) << " tol=" << tol);
        BOOST_CHECK_MESSAGE(std::abs(ad_deriv - fd_deriv) <= tol,
            what << ": AD=" << ad_deriv << " FD=" << fd_deriv
                 << " |diff|=" << std::abs(ad_deriv - fd_deriv) << " tol=" << tol);
    };

    for (int s = 0; s < numDeriv; ++s) {
        const Scalar h = eps[s];
        std::array<Scalar, numDeriv> vp = base_value;
        std::array<Scalar, numDeriv> vm = base_value;
        vp[s] += h;
        vm[s] -= h;

        const auto qp = quantitiesAt(vp[pIdx], vp[z0Idx], vp[z1Idx]);
        const auto qm = quantitiesAt(vm[pIdx], vm[z0Idx], vm[z1Idx]);

        for (int c = 0; c < numComponents; ++c) {
            const Scalar fd_mass = (qp.component_masses[c] - qm.component_masses[c]) / (2.0 * h);
            checkDeriv(qad.component_masses[c].derivative(s), fd_mass, mass_scale,
                       "d(mass[" + std::to_string(c) + "])/dx[" + std::to_string(s) + "]");

            const Scalar fd_mf = (qp.mass_fractions[c] - qm.mass_fractions[c]) / (2.0 * h);
            checkDeriv(qad.mass_fractions[c].derivative(s), fd_mf, 1.0,
                       "d(massfrac[" + std::to_string(c) + "])/dx[" + std::to_string(s) + "]");
        }

        const Scalar fd_rho = (qp.fluid_density - qm.fluid_density) / (2.0 * h);
        checkDeriv(qad.fluid_density.derivative(s), fd_rho, rho_scale,
                   "d(fluid_density)/dx[" + std::to_string(s) + "]");
    }
}
