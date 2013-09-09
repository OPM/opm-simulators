/*===========================================================================
//
// File: test_ad.cpp
//
// Created: 2013-04-29 11:12:34+0200
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <config.h>

#include <opm/autodiff/AutoDiff.hpp>

#include <iostream>

int
main()
try
{
    typedef AutoDiff::Forward<double> AdFW;

    AdFW a = AdFW::variable(0.0);
    AdFW b = AdFW::variable(1.0);

    std::cout << "a:     " << a << '\n';
    std::cout << "b:     " << b << '\n';
    std::cout << "a + b: " << a + b << '\n';

    a = b;
    std::cout << "a:     " << a << '\n';
    std::cout << "b:     " << b << '\n';
    std::cout << "a + b: " << a + b << '\n';

    a = AdFW::variable(0.0);
    std::cout << "a:     " << a << '\n';

    a += 1;
    std::cout << "a:     " << a << '\n';
    std::cout << "a + 1: " << (a + 1) << '\n';
    std::cout << "1 + a: " << (1.0f + a) << '\n';

    a = AdFW::variable(1);
    std::cout << "a:     " << a << '\n';
    std::cout << "a - 1: " << (a - 1) << '\n';
    std::cout << "a - b: " << (a - b) << '\n';
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

