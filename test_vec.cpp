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

#include "AutoDiffVec.hpp"
#include <iostream>

int main()
{
    typedef AutoDiff::ForwardVec<double> ADV;
    ADV::V v(3);
    v << 1.0, 2.2, 3.4;
    std::cout << v << std::endl;
    ADV a = ADV::constant(v);
    ADV x = ADV::variable(v);
    ADV::M jac(3,3);
    jac.insert(0, 0) = 1.0;
    ADV f = ADV::function(v, jac);

    std::cout << a << "\n\n" << x << "\n\n" << f << std::endl;
}
