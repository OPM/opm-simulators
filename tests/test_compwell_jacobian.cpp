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
 * variables and re-flashing. It does the same for the wellbore contents of
 * wellboreContents(), where water fills a volume fraction of the wellbore that
 * is a primary variable as well.
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

// Dimensionless SSHIFT coefficients for CO2, C1 and C10.
constexpr std::array<Scalar, 3> noVolumeShift{0.0, 0.0, 0.0};
constexpr std::array<Scalar, 3> volumeShift{-0.0817, -0.1540, 0.0855};

// Register the fixed CO2/Methane/Decane composition with the generic fluid
// system. The component data is shared static state, so this must run before any
// flash.
void registerFluidSystemComponents(const std::array<Scalar, 3>& sshift = noVolumeShift)
{
    using CO2 = Opm::SimpleCO2<Scalar>;
    using C1  = Opm::C1<Scalar>;
    using C10 = Opm::C10<Scalar>;
    using CompParam = FluidSystem::ComponentParam;

    FluidSystem::init();
    FluidSystem::addComponent(CompParam{CO2::name(), CO2::molarMass(), CO2::criticalTemperature(),
                                        CO2::criticalPressure(), CO2::criticalVolume(), CO2::acentricFactor(),
                                        sshift[0]});
    FluidSystem::addComponent(CompParam{C1::name(), C1::molarMass(), C1::criticalTemperature(),
                                        C1::criticalPressure(), C1::criticalVolume(), C1::acentricFactor(),
                                        sshift[1]});
    FluidSystem::addComponent(CompParam{C10::name(), C10::molarMass(), C10::criticalTemperature(),
                                        C10::criticalPressure(), C10::criticalVolume(), C10::acentricFactor(),
                                        sshift[2]});
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
    T oil_saturation{};
    T gas_saturation{};
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
    q.oil_saturation = so;
    q.gas_saturation = sg;
    q.fluid_density = fs.density(FluidSystem::oilPhaseIdx) * so
                    + fs.density(FluidSystem::gasPhaseIdx) * sg;

    return q;
}

// A state inside the two-phase region, away from phase transitions that would
// make the central-difference comparison unreliable.
constexpr Scalar temperature = 300.0;            // K
constexpr Scalar wellbore_volume = 21.6e-3;      // m^3 (matches CompWell)
constexpr Scalar p0 = 10.0e5;                    // Pa
constexpr Scalar z0_0 = 0.5;
constexpr Scalar z1_0 = 0.3;

// Resolve the flash more accurately than the finite-difference perturbations.
constexpr Scalar flash_tolerance = 1.e-8;

// Wellbore quantities at the test state using the registered components.
WellboreQuantities<Scalar> baseWellboreQuantities()
{
    const std::array<Scalar, numComponents> z{z0_0, z1_0, 1.0 - z0_0 - z1_0};
    return computeWellboreQuantities<Scalar>(p0, z, temperature,
                                             wellbore_volume, flash_tolerance);
}

// Compare wellbore AD derivatives with central differences for the current
// component configuration.
void checkWellboreFlashDerivatives()
{

    // --- Analytical (AD) quantities at the base point ---------------------
    Evaluation P  = Evaluation::createVariable(p0, pIdx);
    std::array<Evaluation, numComponents> z;
    z[0] = Evaluation::createVariable(z0_0, z0Idx);
    z[1] = Evaluation::createVariable(z1_0, z1Idx);
    z[2] = 1.0 - z[0] - z[1];

    const auto qad = computeWellboreQuantities<Evaluation>(P, z, temperature,
                                                           wellbore_volume, flash_tolerance);

    // Require both phases so the derivative check exercises phase-volume coupling.
    BOOST_REQUIRE_GT(qad.oil_saturation.value(), 0.0);
    BOOST_REQUIRE_LT(qad.oil_saturation.value(), 1.0);
    BOOST_REQUIRE_GT(qad.gas_saturation.value(), 0.0);
    BOOST_REQUIRE_LT(qad.gas_saturation.value(), 1.0);
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

// The wellbore contents also depend on the water volume fraction w.
// Derivative slot layout: 0 -> pressure, 1 -> z0, 2 -> z1, 3 -> w.
constexpr int numContentsDeriv = numDeriv + 1;
constexpr int wIdx = numDeriv;
constexpr Scalar w0 = 0.3;

using ContentsEvaluation = Opm::DenseAd::Evaluation<Scalar, numContentsDeriv>;

// A slightly compressible water, standing in for the fluid system's water PVT.
template <typename T>
T
waterDensity(const T& pressure)
{
    return 1000.0 * (1.0 + 4.5e-10 * (pressure - 1.0e5));
}

template <typename T>
Opm::WellboreContents<T, numComponents>
wellboreContentsAt(const T& pressure,
                   const std::array<T, numComponents>& z,
                   const T& water_fraction)
{
    auto fs = makeWellboreFluidState<T>(pressure, z, temperature);
    Opm::flashWellboreFluidState(fs, flash_tolerance);
    return Opm::wellboreContents(fs, water_fraction, waterDensity(pressure), wellbore_volume);
}

std::array<ContentsEvaluation, numComponents>
compositionVariables()
{
    std::array<ContentsEvaluation, numComponents> z;
    z[0] = ContentsEvaluation::createVariable(z0_0, z0Idx);
    z[1] = ContentsEvaluation::createVariable(z1_0, z1Idx);
    z[2] = 1.0 - z[0] - z[1];
    return z;
}

void
checkWellboreContentsDerivatives()
{
    const auto pressure = ContentsEvaluation::createVariable(p0, pIdx);
    const auto water_fraction = ContentsEvaluation::createVariable(w0, wIdx);
    const auto ad = wellboreContentsAt(pressure, compositionVariables(), water_fraction);

    BOOST_REQUIRE_GT(ad.water_mass_fraction.value(), 0.0);
    BOOST_REQUIRE_LT(ad.water_mass_fraction.value(), 1.0);
    // the pressure dependence of the water density must reach the water mass
    BOOST_REQUIRE_GT(ad.water_mass.derivative(pIdx), 0.0);

    Scalar mass_scale = ad.water_mass.value();
    for (const auto& mass : ad.component_masses) {
        mass_scale = std::max(mass_scale, mass.value());
    }
    const Scalar density_scale = ad.density.value();

    const std::array<Scalar, numContentsDeriv> base_value {p0, z0_0, z1_0, w0};
    const std::array<Scalar, numContentsDeriv> eps {p0 * 1.e-5, 1.e-4, 1.e-4, 1.e-4};
    // Scale of each primary variable, so that the absolute floor below is
    // small against the derivatives with respect to pressure as well.
    const std::array<Scalar, numContentsDeriv> variable_scale {p0, 1.0, 1.0, 1.0};

    const auto contentsAt = [](const std::array<Scalar, numContentsDeriv>& v) {
        const std::array<Scalar, numComponents> z {v[z0Idx], v[z1Idx], 1.0 - v[z0Idx] - v[z1Idx]};
        return wellboreContentsAt<Scalar>(v[pIdx], z, v[wIdx]);
    };

    for (int s = 0; s < numContentsDeriv; ++s) {
        auto vp = base_value;
        auto vm = base_value;
        vp[s] += eps[s];
        vm[s] -= eps[s];
        const auto qp = contentsAt(vp);
        const auto qm = contentsAt(vm);

        const auto check = [&](const ContentsEvaluation& quantity,
                               const Scalar plus,
                               const Scalar minus,
                               const Scalar scale,
                               const std::string& what) {
            const Scalar ad_deriv = quantity.derivative(s);
            const Scalar fd_deriv = (plus - minus) / (2.0 * eps[s]);
            const Scalar tol = 1.e-3 * std::abs(ad_deriv) + 1.e-6 * scale / variable_scale[s];
            BOOST_CHECK_MESSAGE(std::abs(ad_deriv - fd_deriv) <= tol,
                                what << "/dx[" << s << "]: AD=" << ad_deriv << " FD=" << fd_deriv
                                     << " |diff|=" << std::abs(ad_deriv - fd_deriv)
                                     << " tol=" << tol);
        };

        for (int c = 0; c < numComponents; ++c) {
            check(ad.component_masses[c],
                  qp.component_masses[c],
                  qm.component_masses[c],
                  mass_scale,
                  "d(mass[" + std::to_string(c) + "])");
            check(ad.mass_fractions[c],
                  qp.mass_fractions[c],
                  qm.mass_fractions[c],
                  1.0,
                  "d(massfrac[" + std::to_string(c) + "])");
        }
        check(ad.water_mass, qp.water_mass, qm.water_mass, mass_scale, "d(water mass)");
        check(ad.water_mass_fraction,
              qp.water_mass_fraction,
              qm.water_mass_fraction,
              1.0,
              "d(water massfrac)");
        check(ad.density, qp.density, qm.density, density_scale, "d(density)");
    }
}

// With water alone in the wellbore every derivative of the flash is multiplied
// by 1 - w = 0, so the contents from a scalar flash match those from the AD
// flash, which CompWell relies on to skip the AD flash.
void
checkWaterFilledWellboreContents()
{
    const auto pressure = ContentsEvaluation::createVariable(p0, pIdx);
    const auto water_fraction = ContentsEvaluation::createVariable(1.0, wIdx);
    const auto z = compositionVariables();
    const auto ad_flash = wellboreContentsAt(pressure, z, water_fraction);

    std::array<Scalar, numComponents> z_value;
    for (int c = 0; c < numComponents; ++c) {
        z_value[c] = z[c].value();
    }
    auto fs = makeWellboreFluidState<Scalar>(p0, z_value, temperature);
    Opm::flashWellboreFluidState(fs, flash_tolerance);
    const auto scalar_flash
        = Opm::wellboreContents(fs, water_fraction, waterDensity(pressure), wellbore_volume);

    const auto checkSame = [](const ContentsEvaluation& expected,
                              const ContentsEvaluation& actual,
                              const std::string& what) {
        const Scalar tol = 1.e-12 * std::max(1.0, std::abs(expected.value()));
        BOOST_CHECK_MESSAGE(std::abs(expected.value() - actual.value()) <= tol,
                            what << ": AD flash " << expected.value() << ", scalar flash "
                                 << actual.value());
        for (int s = 0; s < numContentsDeriv; ++s) {
            const Scalar deriv_tol = 1.e-12 * std::max(1.0, std::abs(expected.derivative(s)));
            BOOST_CHECK_MESSAGE(std::abs(expected.derivative(s) - actual.derivative(s))
                                    <= deriv_tol,
                                what << "/dx[" << s << "]: AD flash " << expected.derivative(s)
                                     << ", scalar flash " << actual.derivative(s));
        }
    };

    for (int c = 0; c < numComponents; ++c) {
        checkSame(ad_flash.component_masses[c],
                  scalar_flash.component_masses[c],
                  "mass[" + std::to_string(c) + "]");
        checkSame(ad_flash.mass_fractions[c],
                  scalar_flash.mass_fractions[c],
                  "massfrac[" + std::to_string(c) + "]");
    }
    checkSame(ad_flash.water_mass, scalar_flash.water_mass, "water mass");
    checkSame(ad_flash.water_mass_fraction, scalar_flash.water_mass_fraction, "water massfrac");
    checkSame(ad_flash.density, scalar_flash.density, "density");
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(WellboreFlashDerivatives)
{
    registerFluidSystemComponents();
    checkWellboreFlashDerivatives();
}

BOOST_AUTO_TEST_CASE(WellboreContentsDerivatives)
{
    registerFluidSystemComponents();
    checkWellboreContentsDerivatives();
}

BOOST_AUTO_TEST_CASE(WaterFilledWellboreContents)
{
    registerFluidSystemComponents();
    checkWaterFilledWellboreContents();
}

BOOST_AUTO_TEST_CASE(WellboreFlashDerivativesWithVolumeShift)
{
    // Exercise the derivatives with SSHIFT applied to densities and saturations.
    registerFluidSystemComponents();
    const auto unshifted = baseWellboreQuantities();

    registerFluidSystemComponents(volumeShift);
    const auto shifted = baseWellboreQuantities();

    // Changing only phase densities must not satisfy the saturation check.
    BOOST_TEST_MESSAGE("oil saturation unshifted = " << unshifted.oil_saturation
                       << ", shifted = " << shifted.oil_saturation);
    BOOST_REQUIRE_GT(std::abs(shifted.oil_saturation - unshifted.oil_saturation), 1.e-6);
    BOOST_REQUIRE_GT(std::abs(shifted.gas_saturation - unshifted.gas_saturation), 1.e-6);

    // Require a measurable density change so ignoring SSHIFT cannot pass the
    // derivative check. At this state, the change is about 0.0146 kg/m3;
    // the 0.001 kg/m3 threshold is well below that change.
    BOOST_TEST_MESSAGE("wellbore density unshifted = " << unshifted.fluid_density
                       << ", shifted = " << shifted.fluid_density);
    BOOST_REQUIRE_GT(std::abs(shifted.fluid_density - unshifted.fluid_density), 1.0e-3);

    checkWellboreFlashDerivatives();
}
