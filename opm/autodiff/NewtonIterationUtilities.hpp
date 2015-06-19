/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_NEWTONITERATIONUTILITIES_HEADER_INCLUDED
#define OPM_NEWTONITERATIONUTILITIES_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <vector>

namespace Opm
{

    /// Eliminate a variable via Schur complement.
    /// \param[in]  eqs  set of equations with Jacobians
    /// \param[in]  n    index of equation/variable to eliminate.
    /// \return          new set of equations, one smaller than eqs.
    /// Note: this method requires the eliminated variable to have the same size
    /// as the equation in the corresponding position (that also will be eliminated).
    std::vector< AutoDiffBlock<double> >
    eliminateVariable(const std::vector< AutoDiffBlock<double> >& eqs,
                      const int n);

    /// Recover that value of a variable previously eliminated.
    /// \param[in]  equation          previously eliminated equation.
    /// \param[in]  partial_solution  solution to the remainder system after elimination.
    /// \param[in]  n                 index of equation/variable that was eliminated.
    /// \return                       solution to complete system.
    AutoDiffBlock<double>::V recoverVariable(const AutoDiffBlock<double>& equation,
                                             const AutoDiffBlock<double>::V& partial_solution,
                                             const int n);

    /// Form an elliptic system of equations.
    /// \param[in]       num_phases  the number of fluid phases
    /// \param[in]       eqs         the equations
    /// \param[out]      A           the resulting full system matrix
    /// \param[out]      b           the right hand side
    /// This function will deal with the first num_phases
    /// equations in eqs, and return a matrix A for the full
    /// system that has a elliptic upper left corner, if possible.
    void formEllipticSystem(const int num_phases,
                            const std::vector< AutoDiffBlock<double> >& eqs,
                            Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                            AutoDiffBlock<double>::V& b);


} // namespace Opm

#endif // OPM_NEWTONITERATIONUTILITIES_HEADER_INCLUDED
