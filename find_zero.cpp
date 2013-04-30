/*===========================================================================
//
// File: find_zero.cpp
//
// Created: 2013-04-29 11:58:29+0200
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

#include "AutoDiff.hpp"

#include <iostream>
#include <cmath>

struct Func
{
    template <typename T>
    T operator()(T x) const
    {
#if 1
        T r = std::sqrt(std::cos(x * x) + x) - 1.2;
        return r;
#else
        return x;
        // const int n = 6;
        // double xv[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
        // double yv[6] = { -0.5, -0.3, -0.1, 0.1, 0.3, 0.5 };
        // int interv = -1;
        // for (int i = 0; i < n; ++i) {
        //     if (x < xv[i]) {
        //         interv = i - 1;
        //         break;
        //     }
        // }
        // T t = (x - xv[interv])/(xv[interv+1] - xv[interv]);
        // return (1.0 - t)*yv[interv] + t*yv[interv+1];
#endif
    }
};

// template <class ErrorPolicy = ThrowOnError>
class Newton
{
public:
    /// Implements a scalar Newton solve.
    template <class Functor>
    inline static double solve(const Functor& f,
                               const double initial_guess,
                               const int max_iter,
                               const double tolerance,
                               int& iterations_used)
    {
        double x = initial_guess;
        iterations_used = 0;
        typedef AutoDiff::Forward<double> AD;
        while (std::abs(f(x)) > tolerance && ++iterations_used < max_iter) {
            AD xfad = AD::variable(x);
            AD rfad = f(xfad);
            x = x - rfad.val()/rfad.der();
        }
        return x;
    }
};


int main()
{
    int iter = 0;
    const double atol = 1.0e-13;
    const double soln = Newton::solve(Func(), 0.1, 30, atol, iter);

    std::cout.precision(16);
    std::cout << "Solution is: " << soln
              << "   using " << iter << " iterations." << '\n';
    std::cout << " f(x) = " << Func()(soln) << '\n';
}
