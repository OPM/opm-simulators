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

#include "AutoDiffBlock.hpp"
#include <iostream>

int main()
{
    typedef AutoDiff::ForwardBlock<double> ADV;
    std::vector<int> blocksizes = { 3, 1, 2 };
    int num_blocks = blocksizes.size();
    ADV::V v1(3);
    v1 << 0.2, 1.2, 13.4;
    ADV::V v2(3);
    v2 << 1.0, 2.2, 3.4;
    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };
    ADV a = ADV::constant(FirstVar, v1, blocksizes);
    ADV x = ADV::variable(FirstVar, v2, blocksizes);
    std::vector<ADV::M> jacs(num_blocks);
    for (int i = 0; i < num_blocks; ++i) {
        jacs[i] = ADV::M(blocksizes[FirstVar], blocksizes[i]);
        jacs[i].insert(0,0) = -1.0;
    }
    ADV f = ADV::function(FirstVar, v2, jacs);
    std::cout << a << x << f;
    /*
    ADV xpx = x + x;
    std::cout << xpx;
    ADV xpxpa = x + x + a;
    std::cout << xpxpa;

    std::cout << xpxpa - xpx;

    ADV sqx = x * x;

    std::cout << sqx;

    ADV sqxdx = sqx / x;

    std::cout << sqxdx;

    // std::cout << a << "\n\n" << x << "\n\n" << f << std::endl;
    */
}
