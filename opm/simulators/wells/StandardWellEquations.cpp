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

#if COMPILE_BDA_BRIDGE
#include <opm/simulators/linalg/bda/WellContributions.hpp>
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

template<class Scalar, int numEq>
StandardWellEquations<Scalar,numEq>::
StandardWellEquations(const ParallelWellInfo& parallel_well_info)
    : parallelB_(duneB_, parallel_well_info)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::
init(const int num_cells,
     const int numWellEq,
     const int numPerfs,
     const std::vector<int>& cells)
{
    // setup sparsity pattern for the matrices
    //[A C^T    [x    =  [ res
    // B D] x_well]      res_well]
    // set the size of the matrices
    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, num_cells, numPerfs);
    duneC_.setSize(1, num_cells, numPerfs);

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
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0 ; perf < numPerfs; ++perf) {
        const int cell_idx = cells[perf];
        // the block size is run-time determined now
        duneB_[0][cell_idx].resize(numWellEq, numEq);
    }

         // make the C^T matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        for (int perf = 0; perf < numPerfs; ++perf) {
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0; perf < numPerfs; ++perf) {
        const int cell_idx = cells[perf];
        duneC_[0][cell_idx].resize(numWellEq, numEq);
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
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::apply(const BVector& x, BVector& Ax) const
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

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::apply(BVector& r) const
{
    assert(invDrw_.size() == invDuneD_.N());

    // invDrw_ = invDuneD_ * resWell_
    invDuneD_.mv(resWell_, invDrw_);
    // r = r - duneC_^T * invDrw_
    duneC_.mmtv(invDrw_, r);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::invert()
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

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::solve(BVectorWell& dx_well) const
{
    invDuneD_.mv(resWell_, dx_well);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::solve(const BVectorWell& rhs_well, BVectorWell& x_well) const
{
    invDuneD_.mv(rhs_well, x_well);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    BVectorWell resWell = resWell_;
    // resWell = resWell - B * x
    parallelB_.mmv(x, resWell);
    // xw = D^-1 * resWell
    invDuneD_.mv(resWell, xw);
}

#if COMPILE_BDA_BRIDGE
template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::
extract(const int numStaticWellEq,
        WellContributions& wellContribs) const
{
    std::vector<int> colIndices;
    std::vector<double> nnzValues;
    colIndices.reserve(duneB_.nonzeroes());
    nnzValues.reserve(duneB_.nonzeroes() * numStaticWellEq * numEq);

    // duneC
    for (auto colC = duneC_[0].begin(),
              endC = duneC_[0].end(); colC != endC; ++colC )
    {
        colIndices.emplace_back(colC.index());
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < numEq; ++j) {
                nnzValues.emplace_back((*colC)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions::MatrixType::C,
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
    wellContribs.addMatrix(WellContributions::MatrixType::D,
                           colIndices.data(), nnzValues.data(), 1);

    // duneB
    colIndices.clear();
    nnzValues.clear();
    for (auto colB = duneB_[0].begin(),
              endB = duneB_[0].end(); colB != endB; ++colB )
    {
        colIndices.emplace_back(colB.index());
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < numEq; ++j) {
                nnzValues.emplace_back((*colB)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions::MatrixType::B,
                           colIndices.data(), nnzValues.data(), duneB_.nonzeroes());
}
#endif

template<class Scalar, int numEq>
template<class SparseMatrixAdapter>
void StandardWellEquations<Scalar,numEq>::
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
        const auto row_index = colC.index();

        for (auto colB = duneB_[0].begin(),
                  endB = duneB_[0].end(); colB != endB; ++colB)
        {
            detail::multMatrix(invDuneD_[0][0], (*colB), tmp);
            detail::negativeMultMatrixTransposed((*colC), tmp, tmpMat);
            jacobian.addToBlock(row_index, colB.index(), tmpMat);
        }
    }
}

template<class Scalar, int numEq>
unsigned int StandardWellEquations<Scalar,numEq>::
getNumBlocks() const
{
    return duneB_.nonzeroes();
}

template<class Scalar, int numEq>
template<class PressureMatrix>
void StandardWellEquations<Scalar,numEq>::
extractCPRPressureMatrix(PressureMatrix& jacobian,
                         const BVector& weights,
                         const int pressureVarIndex,
                         const bool use_well_weights,
                         const WellInterfaceGeneric& well,
                         const int bhp_var_index,
                         const WellState& well_state) const
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
    assert(duneC_.M() == weights.size());
    const int welldof_ind = duneC_.M() + well.indexOfWell();
    // do not assume anything about pressure controlled with use_well_weights (work fine with the assumtion also)
    if (!well.isPressureControlled(well_state) || use_well_weights) {
        // make coupling for reservoir to well
        for (auto colC = duneC_[0].begin(),
                  endC = duneC_[0].end(); colC != endC; ++colC) {
            const auto row_ind = colC.index();
            const auto& bw = weights[row_ind];
            double matel = 0;
            assert((*colC).M() == bw.size());
            for (std::size_t i = 0; i < bw.size(); ++i) {
                matel += (*colC)[bhp_var_index][i] * bw[i];
            }

            jacobian[row_ind][welldof_ind] = matel;
            cell_weights += bw;
            nperf += 1;
        }
    }
    cell_weights /= nperf;

    BVectorWell  bweights(1);
    std::size_t blockSz = duneD_[0][0].N();
    bweights[0].resize(blockSz);
    bweights[0] = 0.0;
    double diagElem = 0;
    if (use_well_weights ) {
        // calculate weighs and set diagonal element
        //NB! use this options without treating pressure controlled separated
        //NB! calculate quasiimpes well weights NB do not work well with trueimpes reservoir weights
        double abs_max = 0;
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
            for (std::size_t i = 0; i < cell_weights.size(); ++i) {
                bweights[0][i] = cell_weights[i];
            }
            bweights[0][blockSz-1] = 0.0;
            diagElem = 0.0;
            const auto& locmat = duneD_[0][0];
            for (std::size_t i = 0; i < cell_weights.size(); ++i) {
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
            const auto col_index = colB.index();
            const auto& bw = bweights[0];
            double matel = 0;
            for (std::size_t i = 0; i < bw.size(); ++i) {
                 matel += (*colB)[i][pressureVarIndex] * bw[i];
            }
            jacobian[welldof_ind][col_index] = matel;
        }
    }
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::
sumDistributed(Parallel::Communication comm)
{
  // accumulate resWell_ and duneD_ in parallel to get effects of all perforations (might be distributed)
    wellhelpers::sumDistributedWellEntries(duneD_[0][0], resWell_[0], comm);
}

#define INSTANCE(N) \
template class StandardWellEquations<double,N>; \
template void StandardWellEquations<double,N>:: \
    extract(Linear::IstlSparseMatrixAdapter<MatrixBlock<double,N,N>>&) const; \
template void StandardWellEquations<double,N>:: \
    extractCPRPressureMatrix(Dune::BCRSMatrix<MatrixBlock<double,1,1>>&, \
                             const typename StandardWellEquations<double,N>::BVector&, \
                             const int, \
                             const bool, \
                             const WellInterfaceGeneric&, \
                             const int, \
                             const WellState&) const;

INSTANCE(1)
INSTANCE(2)
INSTANCE(3)
INSTANCE(4)
INSTANCE(5)
INSTANCE(6)

}
