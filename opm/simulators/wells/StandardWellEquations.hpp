/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_EQUATIONS_HEADER_INCLUDED
#define OPM_STANDARDWELL_EQUATIONS_HEADER_INCLUDED

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

namespace Opm
{

class BlackOilDefaultIndexTraits;
template<class Scalar, class Indices> class BlackOilFluidSystem;
template<class FluidSystem, class Indices, class Scalar> class StandardWellAssemble;
class WellContributions;

//! \brief Matrices and vectors for the well.
template<class Indices, class Scalar>
class StandardWellEquations {
public:
    // number of the conservation equations
    static constexpr int numWellConservationEq = Indices::numPhases + Indices::numSolvents;
    // number of the well control equations
    static constexpr int numWellControlEq = 1;
    // number of the well equations that will always be used
    // based on the solution strategy, there might be other well equations be introduced
    static constexpr int numStaticWellEq = numWellConservationEq + numWellControlEq;
    // the index for Bhp in primary variables and also the index of well control equation
    // they both will be the last one in their respective system.
    // TODO: we should have indices for the well equations and well primary variables separately
    static constexpr int Bhp = numStaticWellEq - numWellControlEq;

    // sparsity pattern for the matrices
    //[A C^T    [x       =  [ res
    // B  D ]   x_well]      res_well]

    // the vector type for the res_well and x_well
    using VectorBlockWellType = Dune::DynamicVector<Scalar>;
    using BVectorWell = Dune::BlockVector<VectorBlockWellType>;

    // the matrix type for the diagonal matrix D
    using DiagMatrixBlockWellType = Dune::DynamicMatrix<Scalar>;
    using DiagMatWell = Dune::BCRSMatrix<DiagMatrixBlockWellType>;

    // the matrix type for the non-diagonal matrix B and C^T
    using OffDiagMatrixBlockWellType = Dune::DynamicMatrix<Scalar>;
    using OffDiagMatWell = Dune::BCRSMatrix<OffDiagMatrixBlockWellType>;

    // block vector type
    using BVector = Dune::BlockVector<Dune::FieldVector<Scalar,Indices::numEq>>;

    //! \brief Constructor initializes matrix assembly mode and the parallel helper.
    StandardWellEquations(const ParallelWellInfo& parallel_well_info);

    //! \brief Setup sparsity pattern for the matrices.
    //! \param num_cells Total number of cells
    //! \param numWellEq Number of well equations
    //! \param numPerfs Number of perforations
    //! \param cells Cell indices for perforations
    void init(const int num_cells,
              const int numWellEq,
              const int numPerfs,
              const std::vector<int> cells);

    //! \brief Set all coefficients to 0.
    void clear();

    //! \brief Apply linear operator to vector.
    void apply(const BVector& x, BVector& Ax) const;

    //! \brief Apply linear operator to vector.
    void apply(BVector& r) const;

    //! \brief Apply inverted D matrix to vector.
    void solve(BVectorWell& dx_well) const;

    //! \brief Invert D matrix.
    void invert();

    // xw = inv(D)*(rw - C*x)
    void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

    //! \brief Add the contribution (C, D^-1, B matrices) of this Well to the WellContributions object
    void addWellContribution(WellContributions& wellContribs) const;

    //! \brief Add the contribution (C, D^-1, B matrices) of this Well to the SparseMatrixAdapter
    template<class SparseMatrixAdapter>
    void addWellContributions(SparseMatrixAdapter& jacobian) const;

    //! \brief Get the number of blocks of the C and B matrices, used to allocate memory in a WellContributions object
    unsigned int getNumBlocks() const
    {
        return duneB_.nonzeroes();
    }

    //! \brief Returns a const reference to the residual.
    const BVectorWell& getResidual() const
    {
        return resWell_;
    }

    //! \brief Sum with off-process contribution.
    void sumDistributed(Parallel::Communication comm);

protected:
    friend class StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                      Indices, Scalar>;
    // residuals of the well equations
    BVectorWell resWell_;

    // two off-diagonal matrices
    OffDiagMatWell duneB_;
    OffDiagMatWell duneC_;

    // diagonal matrix for the well
    DiagMatWell duneD_;
    DiagMatWell invDuneD_;

private:
    // several vector used in the matrix calculation
    mutable BVectorWell Bx_;
    mutable BVectorWell invDrw_;

    // Wrapper for the parallel application of B for distributed wells
    wellhelpers::ParallelStandardWellB<Scalar> parallelB_;
};

}

#endif // OPM_STANDARDWELL_EQUATIONS_HEADER_INCLUDED
