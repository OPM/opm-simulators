/*
  Copyright 2021 Equinor ASA

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

#include <config.h>

#define BOOST_TEST_MODULE SolveTransposed3x3
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <dune/istl/bcrsmatrix.hh>

#include <opm/simulators/linalg/bda/CPR.hpp>

BOOST_AUTO_TEST_CASE(testsolvetransposed3x3)
{
    std::vector<double> A = {-0.0973322, 1.18758e-08, -0.000185231,
                             0.125114, 1.24919e-11, 0,
                             -11.3854, 1.32995e-06, 0.065447};

    std::vector<double> b = {0, 1, 0};
    std::vector<double> x(3);
    std::vector<double> x_expected = {63886256.67, 66154278, 180813.715};

    Opm::Accelerator::solve_transposed_3x3(A.data(), b.data(), x.data());

    for (unsigned i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(x[i], x_expected[i], 1e-6);
    }
}
