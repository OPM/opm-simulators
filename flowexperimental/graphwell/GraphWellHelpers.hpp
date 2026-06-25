/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GRAPHWELL_HELPERS_HEADER_INCLUDED
#define OPM_GRAPHWELL_HELPERS_HEADER_INCLUDED

#include <opm/material/densead/Math.hpp>

#include <cmath>

namespace Opm::graphwellhelpers {

/// Haaland friction factor (copied from mswellhelpers; inlined so the GraphWell
/// works with arbitrary Evaluation sizes without extra instantiations).
template <typename ValueType, typename Scalar>
ValueType haalandFormula(const ValueType& re, const Scalar diameter, const Scalar roughness)
{
    const Scalar rel_roughness = roughness / diameter;
    const ValueType value = -3.6 * log10(6.9 / re + std::pow(rel_roughness / 3.7, 10. / 9.));
    return 1.0 / (value * value);
}

/// Friction pressure-loss magnitude over a pipe segment (>= 0). Mirrors
/// Opm::mswellhelpers::frictionPressureLoss.
template <typename ValueType, typename Scalar>
ValueType frictionPressureLoss(const Scalar l, const Scalar diameter,
                               const Scalar area, const Scalar roughness,
                               const ValueType& density,
                               const ValueType& w, const ValueType& mu)
{
    const ValueType re = abs(diameter * w / (area * mu));
    constexpr Scalar re1 = 2000.;
    constexpr Scalar re2 = 4000.;
    if (re < re1)
        return 32. * mu * l * abs(w) / (area * diameter * diameter * density);

    ValueType f;
    if (re > re2) {
        f = haalandFormula(re, diameter, roughness);
    } else {
        constexpr Scalar f1 = 16. / re1;
        const ValueType f2 = haalandFormula(ValueType{re2}, diameter, roughness);
        f = (f2 - f1) / (re2 - re1) * (re - re1) + f1;
    }
    return 2. * f * l * w * w / (area * area * diameter * density);
}

/// Valve constriction pressure-loss magnitude (>= 0). Mirrors
/// Opm::mswellhelpers::valveContrictionPressureLoss.
template <typename ValueType, typename Scalar>
ValueType valveContrictionPressureLoss(const ValueType& mass_rate, const ValueType& density,
                                       const Scalar area_con, const Scalar cv)
{
    const Scalar area = (area_con > 1.e-10 ? area_con : 1.e-10);
    return mass_rate * mass_rate / (2. * density * cv * cv * area * area);
}

/// Velocity head mass_rate^2 / (area^2 * density). Mirrors
/// Opm::mswellhelpers::velocityHead.
template <typename ValueType, typename Scalar>
ValueType velocityHead(const Scalar area, const ValueType& mass_rate, const ValueType& density)
{
    return mass_rate * mass_rate / (area * area * density);
}

/// Water-in-oil emulsion viscosity (copied from mswellhelpers).
template <typename ValueType, typename Scalar>
ValueType WIOEmulsionViscosity(const ValueType& oil_viscosity,
                               const ValueType& water_liquid_fraction,
                               const Scalar max_visco_ratio)
{
    const ValueType temp_value = 1. / (1. - (0.8415 / 0.7480 * water_liquid_fraction));
    const ValueType viscosity_ratio = pow(temp_value, 2.5);
    return (viscosity_ratio <= max_visco_ratio) ? oil_viscosity * viscosity_ratio
                                                : oil_viscosity * max_visco_ratio;
}

/// Oil-in-water emulsion viscosity (copied from mswellhelpers).
template <typename ValueType, typename Scalar>
ValueType OIWEmulsionViscosity(const ValueType& water_viscosity,
                               const ValueType& water_liquid_fraction,
                               const Scalar max_visco_ratio)
{
    const ValueType temp_value = 1. / (1. - (0.6019 / 0.6410) * (1. - water_liquid_fraction));
    const ValueType viscosity_ratio = pow(temp_value, 2.5);
    return (viscosity_ratio <= max_visco_ratio) ? water_viscosity * viscosity_ratio
                                                : water_viscosity * max_visco_ratio;
}

/// Emulsion viscosity for a spiral-ICD (copied from mswellhelpers::emulsionViscosity,
/// but taking the three SICD parameters as scalars instead of the SICD object).
template <typename ValueType, typename Scalar>
ValueType emulsionViscosity(const ValueType& water_fraction, const ValueType& water_viscosity,
                            const ValueType& oil_fraction, const ValueType& oil_viscosity,
                            const Scalar width_transition, const Scalar critical_value,
                            const Scalar max_visco_ratio)
{
    const ValueType transition_start_value = critical_value - width_transition / 2.0;
    const ValueType transition_end_value = critical_value + width_transition / 2.0;

    const ValueType liquid_fraction = water_fraction + oil_fraction;
    if (liquid_fraction == 0.)
        return ValueType{0.};

    const ValueType water_liquid_fraction = water_fraction / liquid_fraction;

    if (water_liquid_fraction <= transition_start_value) {
        return WIOEmulsionViscosity(oil_viscosity, water_liquid_fraction, max_visco_ratio);
    } else if (water_liquid_fraction >= transition_end_value) {
        return OIWEmulsionViscosity(water_viscosity, water_liquid_fraction, max_visco_ratio);
    } else {
        const ValueType v_start = WIOEmulsionViscosity(oil_viscosity, transition_start_value, max_visco_ratio);
        const ValueType v_end = OIWEmulsionViscosity(water_viscosity, transition_end_value, max_visco_ratio);
        return (v_start * (transition_end_value - water_liquid_fraction)
                + v_end * (water_liquid_fraction - transition_start_value)) / width_transition;
    }
}

} // namespace Opm::graphwellhelpers

#endif // OPM_GRAPHWELL_HELPERS_HEADER_INCLUDED
