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

#ifndef OPM_MULTISEGMENTWELL_EQUATIONS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_EQUATIONS_HEADER_INCLUDED

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <memory>

namespace Dune {
template<class M> class UMFPack;
}

namespace Opm
{

template<typename FluidSystem, typename Indices, int numWellEq, int numEq> class MultisegmentWellEquationAccess;
template<typename FluidSystem, typename Indices> class MultisegmentWellGeneric;
#if COMPILE_GPU_BRIDGE
template<class Scalar> class WellContributions;
#endif
template<typename FluidSystem, typename Indices> class WellInterfaceGeneric;
template<typename FluidSystem, typename Indices> class WellState;

// TODO: numWelEq and numEq should be able to related to the Indices
template<typename FluidSystem, typename Indices, int numWellEq, int numEq>
class MultisegmentWellEquations
{
public:
    using Scalar = typename FluidSystem::Scalar;
    // sparsity pattern for the matrices
    // [A C^T    [x       =  [ res
    //  B  D ]   x_well]      res_well]

    // the vector type for the res_well and x_well
    using VectorBlockWellType = Dune::FieldVector<Scalar,numWellEq>;
    using BVectorWell = Dune::BlockVector<VectorBlockWellType>;

    using VectorBlockType = Dune::FieldVector<Scalar,numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;

    // the matrix type for the diagonal matrix D
    using DiagMatrixBlockWellType = Dune::FieldMatrix<Scalar,numWellEq,numWellEq>;
    using DiagMatWell = Dune::BCRSMatrix<DiagMatrixBlockWellType>;

    // the matrix type for the non-diagonal matrix B and C^T
    using OffDiagMatrixBlockWellType = Dune::FieldMatrix<Scalar,numWellEq,numEq>;
    using OffDiagMatWell = Dune::BCRSMatrix<OffDiagMatrixBlockWellType>;

    MultisegmentWellEquations(const MultisegmentWellGeneric<FluidSystem, Indices>& well, const ParallelWellInfo<Scalar>& pw_info);

    //! \brief Setup sparsity pattern for the matrices.
    //! \param numPerfs Number of perforations
    //! \param cells Cell indices for perforations
    //! \param segment_inlets Cell indices for segment inlets
    //! \param segment_perforations Cell indices for segment perforations
    void init(const int numPerfs,
              const std::vector<int>& cells,
              const std::vector<std::vector<int>>& segment_inlets,
              const std::vector<std::vector<int>>& segment_perforations);

    //! \brief Set all coefficients to 0.
    void clear();

    //! \brief Apply linear operator to vector.
    void apply(const BVector& x, BVector& Ax) const;

    //! \brief Apply linear operator to vector.
    void apply(BVector& r) const;

    //! \brief Compute the LU-decomposition of D matrix.
    void createSolver();

    //! \brief Apply inverted D matrix to residual and return result.
    BVectorWell solve() const;

    //! \brief Apply inverted D matrix to rhs and return result.
    BVectorWell solve(const BVectorWell& rhs) const;

    //! \brief Recover well solution.
    //! \details xw = inv(D)*(rw - C*x)
    void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

#if COMPILE_GPU_BRIDGE
    //! \brief Add the matrices of this well to the WellContributions object.
    void extract(WellContributions<Scalar>& wellContribs) const;
#endif

    //! \brief Add the matrices of this well to the sparse matrix adapter.
    template<class SparseMatrixAdapter>
    void extract(SparseMatrixAdapter& jacobian) const;

    //! \brief Extract CPR pressure matrix.
    template<class PressureMatrix>
    void extractCPRPressureMatrix(PressureMatrix& jacobian,
                                  const BVector& weights,
                                  const int pressureVarIndex,
                                  const bool /*use_well_weights*/,
                                  const WellInterfaceGeneric<FluidSystem, Indices>& well,
                                  const int seg_pressure_var_ind,
                                  const WellState<FluidSystem, Indices>& well_state) const;

    //! \brief Sum with off-process contribution.
    void sumDistributed(Parallel::Communication comm);

    //! \brief Returns a const reference to the residual.
    const BVectorWell& residual() const
    {
        return resWell_;
    }

  private:
    friend class MultisegmentWellEquationAccess<FluidSystem, Indices, numWellEq, numEq>;
    // two off-diagonal matrices
    OffDiagMatWell duneB_;
    OffDiagMatWell duneC_;
    // "diagonal" matrix for the well. It has offdiagonal entries for inlets and outlets.
    DiagMatWell duneD_;

    /// \brief solver for diagonal matrix
    ///
    /// This is a shared_ptr as MultisegmentWell is copied in computeWellPotentials...
    mutable std::shared_ptr<Dune::UMFPack<DiagMatWell>> duneDSolver_;

    // residuals of the well equations
    BVectorWell resWell_;

    const  MultisegmentWellGeneric<FluidSystem, Indices>& well_; //!< Reference to well

    // Store the global index of well perforated cells
    std::vector<int> cells_;

    const ParallelWellInfo<Scalar>& pw_info_;

    // Wrapper for the parallel application of B for distributed wells
    mswellhelpers::ParallellMSWellB<OffDiagMatWell> parallelB_;
};

}

#include "MultisegmentWellEquations.cpp"

#endif // OPM_MULTISEGMENTWELLWELL_EQUATIONS_HEADER_INCLUDED
