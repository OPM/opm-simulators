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
/// prototype works with arbitrary Evaluation sizes without extra instantiations).
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

/// Velocity head mass_rate^2 / (area^2 * density). Mirrors
/// Opm::mswellhelpers::velocityHead.
template <typename ValueType, typename Scalar>
ValueType velocityHead(const Scalar area, const ValueType& mass_rate, const ValueType& density)
{
    return mass_rate * mass_rate / (area * area * density);
}

} // namespace Opm::graphwellhelpers

#endif // OPM_GRAPHWELL_HELPERS_HEADER_INCLUDED
