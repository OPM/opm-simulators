/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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


#ifndef OPM_MSWELLHELPERS_HEADER_INCLUDED
#define OPM_MSWELLHELPERS_HEADER_INCLUDED

#include <dune/istl/solvers.hh>

namespace Opm {

namespace mswellhelpers
{

    // obtain y = D^-1 * x
    template<typename MatrixType, typename VectorType>
    VectorType
    invDX(const MatrixType& D, VectorType x)
    {
        // the function will change the value of x, so we should not use reference of x here.

        // TODO: store some of the following information to avoid to call it again and again for
        // efficiency improvement.
        // Bassically, only the solve / apply step is different.

        VectorType y(x.size());
        y = 0.;

        Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(D);

        // Sequential incomplete LU decomposition as the preconditioner
        Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);

        // Preconditioned BICGSTAB solver
        Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
                                                   preconditioner,
                                                   1.e-6, // desired residual reduction factor
                                                   50, // maximum number of iterations
                                                   0); // verbosity of the solver

        // Object storing some statistics about the solving process
        Dune::InverseOperatorResult statistics ;

        // Solve
        linsolver.apply(y, x, statistics );

        return y;
    }


} // namespace mswellhelpers

}

#endif
