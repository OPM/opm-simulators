/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_WELL_EQUATIONS_HPP
#define OPM_COMP_WELL_EQUATIONS_HPP

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

namespace Opm {

// we should look into whether to use dynamic matrix and vector here
template <typename Scalar, int numWellEq, int numEq>
class CompWellEquations
{
public:
    // sparsity pattern for the matrices
    //[A C^T    [x       =  [ res
    // B  D ]   x_well]      res_well]

    CompWellEquations();

    // the vector type for the res_well and x_well
    using VectorBlockWellType = Dune::FieldVector<Scalar, numWellEq>;
    using BVectorWell = Dune::BlockVector<VectorBlockWellType>;

    // for res
    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;

    // the matrix type for the diagonal matrix D
    using DiagMatrixBlockWellType = Dune::FieldMatrix<Scalar, numWellEq>;
    using DiagMatWell = Dune::BCRSMatrix<DiagMatrixBlockWellType>;

    // the matrix type for the non-diagonal matrix B and C^T
    using OffDiagMatrixBlockWellType = Dune::FieldMatrix<Scalar, numWellEq, numEq>;
    using OffDiagMatWell = Dune::BCRSMatrix<OffDiagMatrixBlockWellType>;

    void init(const int num_conn, const std::vector<std::size_t>& cells);

    void clear();

    DiagMatWell& D()
    {
        return duneD_;
    }

    OffDiagMatWell& B()
    {
        return duneB_;
    }

    OffDiagMatWell& C()
    {
        return duneC_;
    }

    BVectorWell& residual()
    {
        return resWell_;
    }

    const BVectorWell& residual() const
    {
        return resWell_;
    }

    void solve(BVectorWell& dx_well) const;

    void invert();

    void apply(BVector& r) const;

    void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

private:
    // two off-diagonal matrices
    OffDiagMatWell duneB_;
    OffDiagMatWell duneC_;

    // diagonal matrix for the well
    DiagMatWell invDuneD_;
    DiagMatWell duneD_;

    // residuals of the well equations
    BVectorWell resWell_;

    // several vector used in the matrix calculation
    mutable BVectorWell Bx_;
    mutable BVectorWell invDrw_;


    // Store the global index of the well connection cells
    std::vector<std::size_t> cells_;
};

} // end of namespace Opm

#include "CompWellEquations_impl.hpp"

#endif // OPM_COMP_WELL_EQUATIONS_HPP
