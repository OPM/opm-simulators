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

#include <config.h>
#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/wells/StandardWellEquations.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#if COMPILE_GPU_BRIDGE
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>
#endif

#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace Opm
{

template<typename Scalar, typename IndexTraits, int numEq>
StandardWellEquations<Scalar, IndexTraits, numEq>::
StandardWellEquations(const ParallelWellInfo<Scalar>& parallel_well_info)
    : parallelB_(duneB_, parallel_well_info)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
init(const int numWellEq,
     const int numPerfs,
     const std::vector<int>& cells)
{
    // setup sparsity pattern for the matrices
    //[A C^T    [x    =  [ res
    // B D] x_well]      res_well]
    // set the size of the matrices
    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, numPerfs, numPerfs);
    duneC_.setSize(1, numPerfs, numPerfs);

    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }
      // the block size is run-time determined now
    duneD_[0][0].resize(numWellEq, numWellEq);

    for (auto row = duneB_.createbegin(),
              end = duneB_.createend(); row != end; ++row) {
        for (int perf = 0 ; perf < numPerfs; ++perf) {
            row.insert(perf);
        }
    }

    for (int perf = 0 ; perf < numPerfs; ++perf) {
        // the block size is run-time determined now
        duneB_[0][perf].resize(numWellEq, numEq);
    }

         // make the C^T matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        for (int perf = 0; perf < numPerfs; ++perf) {
            row.insert(perf);
        }
    }

    for (int perf = 0; perf < numPerfs; ++perf) {
        duneC_[0][perf].resize(numWellEq, numEq);
    }

    resWell_.resize(1);
    // the block size of resWell_ is also run-time determined now
    resWell_[0].resize(numWellEq);

    // resize temporary class variables
    Bx_.resize(duneB_.N());
    for (unsigned i = 0; i < duneB_.N(); ++i) {
        Bx_[i].resize(numWellEq);
    }

    invDrw_.resize(duneD_.N());
    for (unsigned i = 0; i < duneD_.N(); ++i) {
        invDrw_[i].resize(numWellEq);
    }

    // Store the global index of well perforated cells
    cells_ = cells;
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::apply(const BVector& x, BVector& Ax) const
{
    assert(Bx_.size() == duneB_.N());
    assert(invDrw_.size() == invDuneD_.N());

    // Bx_ = duneB_ * x
    parallelB_.mv(x, Bx_);

    // invDBx = invDuneD_ * Bx_
    // TODO: with this, we modified the content of the invDrw_.
    // Is it necessary to do this to save some memory?
    auto& invDBx = invDrw_;
    invDuneD_.mv(Bx_, invDBx);

    // Ax = Ax - duneC_^T * invDBx
    duneC_.mmtv(invDBx, Ax);
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::apply(BVector& r) const
{
    assert(invDrw_.size() == invDuneD_.N());

    // invDrw_ = invDuneD_ * resWell_
    invDuneD_.mv(resWell_, invDrw_);
    // r = r - duneC_^T * invDrw_
    duneC_.mmtv(invDrw_, r);
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::invert()
{
    try {
        invDuneD_ = duneD_; // Not strictly need if not cpr with well contributions is used
        detail::invertMatrix(invDuneD_[0][0]);
    } catch (NumericalProblem&) {
        // for singular matrices, use identity as the inverse
        invDuneD_[0][0] = 0.0;
        for (std::size_t i = 0; i < invDuneD_[0][0].rows(); ++i) {
            invDuneD_[0][0][i][i] = 1.0;
        }
    }
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::solve(BVectorWell& dx_well) const
{
    invDuneD_.mv(resWell_, dx_well);
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::solve(const BVectorWell& rhs_well, BVectorWell& x_well) const
{
    invDuneD_.mv(rhs_well, x_well);
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    BVectorWell resWell = resWell_;
    // resWell = resWell - B * x
    parallelB_.mmv(x, resWell);
    // xw = D^-1 * resWell
    invDuneD_.mv(resWell, xw);
}

#if COMPILE_GPU_BRIDGE
template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
extract(const int numStaticWellEq,
        WellContributions<Scalar>& wellContribs) const
{
    std::vector<int> colIndices;
    std::vector<Scalar> nnzValues;
    colIndices.reserve(duneB_.nonzeroes());
    nnzValues.reserve(duneB_.nonzeroes() * numStaticWellEq * numEq);

    // duneC
    for (auto colC = duneC_[0].begin(),
              endC = duneC_[0].end(); colC != endC; ++colC )
    {
        colIndices.emplace_back(cells_[colC.index()]);
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < numEq; ++j) {
                nnzValues.emplace_back((*colC)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions<Scalar>::MatrixType::C,
                           colIndices.data(), nnzValues.data(), duneC_.nonzeroes());

    // invDuneD
    colIndices.clear();
    nnzValues.clear();
    colIndices.emplace_back(0);
    for (int i = 0; i < numStaticWellEq; ++i)
    {
        for (int j = 0; j < numStaticWellEq; ++j) {
            nnzValues.emplace_back(invDuneD_[0][0][i][j]);
        }
    }
    wellContribs.addMatrix(WellContributions<Scalar>::MatrixType::D,
                           colIndices.data(), nnzValues.data(), 1);

    // duneB
    colIndices.clear();
    nnzValues.clear();
    for (auto colB = duneB_[0].begin(),
              endB = duneB_[0].end(); colB != endB; ++colB )
    {
        colIndices.emplace_back(cells_[colB.index()]);
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < numEq; ++j) {
                nnzValues.emplace_back((*colB)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions<Scalar>::MatrixType::B,
                           colIndices.data(), nnzValues.data(), duneB_.nonzeroes());
}
#endif

template<typename Scalar, typename IndexTraits, int numEq>
template<class SparseMatrixAdapter>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
extract(SparseMatrixAdapter& jacobian) const
{
    // We need to change matrx A as follows
    // A -= C^T D^-1 B
    // D is diagonal
    // B and C have 1 row, nc colums and nonzero
    // at (0,j) only if this well has a perforation at cell j.
    typename SparseMatrixAdapter::MatrixBlock tmpMat;
    Dune::DynamicMatrix<Scalar> tmp;
    for (auto colC = duneC_[0].begin(),
              endC = duneC_[0].end(); colC != endC; ++colC)
    {
        // map the well perforated cell index to global cell index
        const auto row_index = this->cells_[colC.index()];

        for (auto colB = duneB_[0].begin(),
                  endB = duneB_[0].end(); colB != endB; ++colB)
        {
            // map the well perforated cell index to global cell index
            const auto col_index = this->cells_[colB.index()];
            detail::multMatrix(invDuneD_[0][0], (*colB), tmp);
            detail::negativeMultMatrixTransposed((*colC), tmp, tmpMat);
            jacobian.addToBlock(row_index, col_index, tmpMat);
        }
    }
}

template<typename Scalar, typename IndexTraits, int numEq>
unsigned int StandardWellEquations<Scalar, IndexTraits, numEq>::
getNumBlocks() const
{
    return duneB_.nonzeroes();
}

template<typename Scalar, typename IndexTraits, int numEq>
template<class PressureMatrix>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
extractCPRPressureMatrix(PressureMatrix& jacobian,
                         const BVector& weights,
                         const int pressureVarIndex,
                         const bool use_well_weights,
                         const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                         const int bhp_var_index,
                         const WellState<Scalar, IndexTraits>& well_state) const
{
    // This adds pressure quation for cpr
    // For use_well_weights=true
    //    weights lamda = inv(D)'v  v = 0 v(bhpInd) = 1
    //    the well equations are summed i lambda' B(:,pressureVarINd) -> B  lambda'*D(:,bhpInd) -> D
    // For use_well_weights = false
    //    weights lambda = \sum_i w /n where ths sum is over weights of all perforation cells
    //    in the case of pressure controlled trivial equations are used and bhp  C=B=0
    //    then the flow part of the well equations are summed lambda'*B(1:n,pressureVarInd) -> B lambda'*D(1:n,bhpInd) -> D
    // For bouth
    //    C -> w'C(:,bhpInd) where w is weights of the perforation cell

    // Add the well contributions in cpr
    // use_well_weights is a quasiimpes formulation which is not implemented in multisegment
    int nperf = 0;
    auto cell_weights = weights[0];// not need for not(use_well_weights)
    cell_weights = 0.0;
    const int number_cells = weights.size();
    const int welldof_ind = number_cells + well.indexOfWell();
    // do not assume anything about pressure controlled with use_well_weights (work fine with the assumtion also)
    if (!well.isPressureControlled(well_state) || use_well_weights) {
        // make coupling for reservoir to well
        for (auto colC = duneC_[0].begin(),
                  endC = duneC_[0].end(); colC != endC; ++colC) {
            // map the well perforated cell index to global cell index
            const auto row_index = cells_[colC.index()];
            const auto& bw = weights[row_index];
            Scalar matel = 0;
            assert((*colC).M() == bw.size());
            for (std::size_t i = 0; i < bw.size(); ++i) {
                matel += (*colC)[bhp_var_index][i] * bw[i];
            }

            jacobian[row_index][welldof_ind] = matel;
            cell_weights += bw;
            nperf += 1;
        }
    }
    if (nperf != 0)
        cell_weights /= nperf;
    else {
        // duneC_.size==0, which can happen e.g. if a well has no active perforations on this rank.
        // Add positive weight to diagonal to regularize Jacobian (other row entries are 0).
        // Row's variable has no observable effect, since there are no perforations.
        cell_weights = 1.;
    }

    BVectorWell  bweights(1);
    std::size_t blockSz = duneD_[0][0].N();
    bweights[0].resize(blockSz);
    bweights[0] = 0.0;
    Scalar diagElem = 0;
    if (use_well_weights ) {
        // calculate weighs and set diagonal element
        //NB! use this options without treating pressure controlled separated
        //NB! calculate quasiimpes well weights NB do not work well with trueimpes reservoir weights
        Scalar abs_max = 0;
        BVectorWell rhs(1);
        rhs[0].resize(blockSz);
        rhs[0][bhp_var_index] = 1.0;
        DiagMatrixBlockWellType inv_diag_block = invDuneD_[0][0];
        DiagMatrixBlockWellType inv_diag_block_transpose =
            Opm::wellhelpers::transposeDenseDynMatrix(inv_diag_block);
        for (std::size_t i = 0; i < blockSz; ++i) {
            bweights[0][i] = 0;
            for (std::size_t j = 0; j < blockSz; ++j) {
                bweights[0][i] += inv_diag_block_transpose[i][j] * rhs[0][j];
            }
            abs_max = std::max(abs_max, std::fabs(bweights[0][i]));
        }
        assert(abs_max > 0.0);
        for (std::size_t i = 0; i < blockSz; ++i) {
            bweights[0][i] /= abs_max;
        }
        diagElem = 1.0 / abs_max;
    } else {
        // set diagonal element
        if (well.isPressureControlled(well_state)) {
            bweights[0][blockSz-1] = 1.0;
            diagElem = 1.0; // better scaling could have used the calculation below if weights were calculated
        } else {
            // The first (blockSz - 1) block weights will scale the
            // conservation well equations, and is therefore set equal
            // to the cell weights. The last one will scale the control
            // equation, we set that to zero.
            for (std::size_t i = 0; i < blockSz - 1; ++i) {
                bweights[0][i] = cell_weights[i];
            }
            bweights[0][blockSz-1] = 0.0;
            diagElem = 0.0;
            const auto& locmat = duneD_[0][0];
            // For some models such as MICP, cell_weights.size() is larger than
            // (blockSz - 1) since the model has more conserved quantities than
            // the well model treats. We assume that the first (blockSz - 1)
            // conserved quantities correspond to those treated in the well model.
            for (std::size_t i = 0; i < blockSz - 1; ++i) {
                diagElem += locmat[i][bhp_var_index] * cell_weights[i];
            }

        }
    }
    //
    jacobian[welldof_ind][welldof_ind] = diagElem;
    // set the matrix elements for well reservoir coupling
    if (!well.isPressureControlled(well_state) || use_well_weights) {
        for (auto colB = duneB_[0].begin(),
                  endB = duneB_[0].end(); colB != endB; ++colB) {
            // map the well perforated cell index to global cell index
            const auto col_index = cells_[colB.index()];
            const auto& bw = bweights[0];
            Scalar matel = 0;
            for (std::size_t i = 0; i < bw.size(); ++i) {
                 matel += (*colB)[i][pressureVarIndex] * bw[i];
            }
            jacobian[welldof_ind][col_index] = matel;
        }
    }
}

template<typename Scalar, typename IndexTraits, int numEq>
void StandardWellEquations<Scalar, IndexTraits, numEq>::
sumDistributed(Parallel::Communication comm)
{
  // accumulate resWell_ and duneD_ in parallel to get effects of all perforations (might be distributed)
    wellhelpers::sumDistributedWellEntries(duneD_[0][0], resWell_[0], comm);
}

#define INSTANTIATE(T,N)                                                              \
    template class StandardWellEquations<T,BlackOilDefaultFluidSystemIndices,N>;                                        \
    template void StandardWellEquations<T,BlackOilDefaultFluidSystemIndices,N>::                                        \
        extract(Linear::IstlSparseMatrixAdapter<MatrixBlock<T,N,N>>&) const;          \
    template void StandardWellEquations<T,BlackOilDefaultFluidSystemIndices,N>::                                        \
        extractCPRPressureMatrix(Dune::BCRSMatrix<MatrixBlock<T,1,1>>&,               \
                                 const typename StandardWellEquations<T,BlackOilDefaultFluidSystemIndices,N>::BVector&, \
                                 const int,                                           \
                                 const bool,                                          \
                                 const WellInterfaceGeneric<T,BlackOilDefaultFluidSystemIndices>&,                      \
                                 const int,                                           \
                                 const WellState<T,BlackOilDefaultFluidSystemIndices>&) const;

#define INSTANTIATE_TYPE(T) \
    INSTANTIATE(T,1)        \
    INSTANTIATE(T,2)        \
    INSTANTIATE(T,3)        \
    INSTANTIATE(T,4)        \
    INSTANTIATE(T,5)        \
    INSTANTIATE(T,6)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
