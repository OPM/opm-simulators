/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#include <opm/autodiff/LinearisedBlackoilResidual.hpp>

int Opm::LinearisedBlackoilResidual::sizeNonLinear() const
{
    int size = 0;

    std::vector<ADB>::const_iterator massBalanceIt = material_balance_eq.begin();
    const std::vector<ADB>::const_iterator endMassBalanceIt = material_balance_eq.end();
    for (; massBalanceIt != endMassBalanceIt; ++massBalanceIt) {
        size += (*massBalanceIt).size();
    }

    size += well_flux_eq.size();
    size += well_eq.size();

    return size;
}
