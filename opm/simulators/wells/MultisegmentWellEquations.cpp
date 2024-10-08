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
#include <opm/simulators/wells/MultisegmentWellEquations.hpp>

#include <dune/istl/umfpack.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#if COMPILE_BDA_BRIDGE
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#endif

#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>

#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <cstddef>
#include <stdexcept>
//define COMMENTS

namespace Opm {

template<class Scalar, int numWellEq, int numEq>
MultisegmentWellEquations<Scalar,numWellEq,numEq>::
MultisegmentWellEquations(const MultisegmentWellGeneric<Scalar>& well, const ParallelWellInfo<Scalar>& pw_info)
    : well_(well)
    , pw_info_(pw_info)
{
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
init(const int num_cells_this_process,
     const int num_cells_all_processes_of_this_well,
     const int num_perfs_this_process,
     const int num_perfs_whole_mswell,
     const std::vector<int>& cells,
     const std::vector<int>& cells_whole_mswell,
     const std::vector<std::vector<int>>& segment_inlets,
     const std::vector<std::vector<int>>& perforations)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneBGlobal_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    duneCGlobal_.setBuildMode(OffDiagMatWell::row_wise);
    duneD_.setBuildMode(DiagMatWell::row_wise);

    // set the size and patterns for all the matrices and vectors
    // [A C^T    [x    =  [ res
    //  B D] x_well]      res_well]

    // calculating the NNZ for duneD_
    // NNZ = number_of_segments + 2 * (number_of_inlets / number_of_outlets)
    {
        int nnz_d = well_.numberOfSegments();
        for (const std::vector<int>& inlets : segment_inlets) {
            nnz_d += 2 * inlets.size();
        }
        duneD_.setSize(well_.numberOfSegments(), well_.numberOfSegments(), nnz_d);
    }

    // The matrices B and C have *space for* segments x num_cells_all_processes_of_this_well entries and num_perfs_whole_mswell nonzero entries
    // We will later fill these spaces on the respective processes, then collect everything on the lowest rank that the well perforates,
    // solve the Schur Complement corresponding to the MSWellequations on that rank and distribute the solution back to the other
    // ranks of the well.
    duneB_.setSize(well_.numberOfSegments(), num_cells_this_process, num_perfs_this_process);
    // Todo: Make the global matrix smaller than: num_cells_all_processes_of_this_well, because that includes the overlap part as well!
    duneBGlobal_.setSize(well_.numberOfSegments(), num_cells_all_processes_of_this_well, num_perfs_whole_mswell); 
    duneC_.setSize(well_.numberOfSegments(), num_cells_this_process, num_perfs_this_process);
    duneCGlobal_.setSize(well_.numberOfSegments(), num_cells_all_processes_of_this_well, num_perfs_whole_mswell); 

    // we need to add the off diagonal ones
    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // the number of the row corrspnds to the segment now
        const int seg = row.index();
        // adding the item related to outlet relation
        const Segment& segment = well_.segmentSet()[seg];
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) { // if there is a outlet_segment
            const int outlet_segment_index = well_.segmentNumberToIndex(outlet_segment_number);
            row.insert(outlet_segment_index);
        }

        // Add nonzeros for diagonal
        row.insert(seg);

        // insert the item related to its inlets
        for (const int& inlet : segment_inlets[seg]) {
            row.insert(inlet);
        }
    }

    // make the C matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : perforations[row.index()]) {
            const int local_perf_index = pw_info_.getGlobalPerfContainerFactory().globalToLocal(perf);
            if (local_perf_index < 0) // then the perforation is not on this process
                continue;
            const int cell_idx = cells[local_perf_index];
            row.insert(cell_idx);
        }
    }

    // make the global C matrix
    for (auto row = duneCGlobal_.createbegin(),
              end = duneCGlobal_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : perforations[row.index()]) {
            const int cell_idx = cells_whole_mswell[perf];
            row.insert(cell_idx);
        }
    }

    // make the B^T matrix
    for (auto row = duneB_.createbegin(),
              end = duneB_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : perforations[row.index()]) {
            const int local_perf_index = pw_info_.getGlobalPerfContainerFactory().globalToLocal(perf);
            if (local_perf_index < 0) // then the perforation is not on this process
                continue;
            const int cell_idx = cells[local_perf_index];
            row.insert(cell_idx);
        }
    }
    // make the global B matrix
    for (auto row = duneBGlobal_.createbegin(),
              end = duneBGlobal_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : perforations[row.index()]) {
            const int cell_idx = cells_whole_mswell[perf];
            row.insert(cell_idx);
        }
    }

    resWell_.resize(well_.numberOfSegments());
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::clear()
{
    duneB_ = 0.0;
    duneBGlobal_ = 0.0;
    duneC_ = 0.0;
    duneCGlobal_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
    duneDSolver_.reset();
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
apply(const BVector& x, BVector& Ax) const
{

#ifdef COMMENTS
    std::cout << "In apply, x:" << std::endl;
    for (int i = 0; i < x.size(); i++)
        std::cout << i << ": " << x[i] << std::endl;
    std::cout << "Beginning: Ax:" << std::endl; 
    for (int i = 0; i < Ax.size(); i++)
        std::cout << i << ": " << Ax[i] << std::endl;
#endif

    BVectorWell Bx(duneB_.N());

    duneB_.mv(x, Bx);

    // We need to communicate here to get the contributions from all segments
    this->pw_info_.communication().sum(Bx.data(), Bx.size());

    // invDBx = duneD^-1 * Bx_
    // TODO: do this only on one process and distribute the solution
    const BVectorWell invDBx = mswellhelpers::applyUMFPack(*duneDSolver_, Bx);

    // Ax = Ax - duneC_^T * invDBx
    duneC_.mmtv(invDBx,Ax);

#ifdef COMMENTS
    std::cout << "At the end: Ax:" << std::endl; 
    for (int i = 0; i < Ax.size(); i++)
        std::cout << i << ": " << Ax[i] << std::endl;
#endif
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
apply(BVector& r) const
{
    // invDrw_ = duneD^-1 * resWell_
    // TODO: do this only on one process and distribute the solution
    const BVectorWell invDrw = mswellhelpers::applyUMFPack(*duneDSolver_, resWell_);
    // r = r - duneC_^T * invDrw
    duneC_.mmtv(invDrw, r);
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::createSolver()
{
#if HAVE_UMFPACK
    if (duneDSolver_) {
        return;
    }

    if constexpr (std::is_same_v<Scalar,float>) {
        OPM_THROW(std::runtime_error, "MultisegmentWell support requires UMFPACK, "
                                      "and UMFPACK does not support float");
    } else {
        duneDSolver_ = std::make_shared<Dune::UMFPack<DiagMatWell>>(duneD_, 0);
    }
#else
    OPM_THROW(std::runtime_error, "MultisegmentWell support requires UMFPACK. "
              "Reconfigure opm-simulators with SuiteSparse/UMFPACK support and recompile.");
#endif
}

template<class Scalar, int numWellEq, int numEq>
typename MultisegmentWellEquations<Scalar,numWellEq,numEq>::BVectorWell
MultisegmentWellEquations<Scalar,numWellEq,numEq>::solve() const
{
    // TODO: do this only on one process and distribute the solution
    return mswellhelpers::applyUMFPack(*duneDSolver_, resWell_);
}

template<class Scalar, int numWellEq, int numEq>
typename MultisegmentWellEquations<Scalar,numWellEq,numEq>::BVectorWell
MultisegmentWellEquations<Scalar,numWellEq,numEq>::solve(const BVectorWell& rhs) const
{
    // TODO: do this only on one process and distribute the solution
    return mswellhelpers::applyUMFPack(*duneDSolver_, rhs);
}

template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    BVectorWell resWell = resWell_;
    BVectorWell Bx(duneB_.N());
    duneB_.mv(x, Bx);
    // We need to communicate here to get the contributions from all segments
    this->pw_info_.communication().sum(Bx.data(), Bx.size());

    // resWell = resWell - B * x
    resWell -= Bx;

    // xw = D^-1 * resWell
    // TODO: do this only on one process and distribute the solution
    xw = mswellhelpers::applyUMFPack(*duneDSolver_, resWell);
}

#if COMPILE_BDA_BRIDGE
template<class Scalar, int numWellEq, int numEq>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
extract(WellContributions<Scalar>& wellContribs) const
{
    unsigned int Mb = duneB_.N();       // number of blockrows in duneB_, duneC_ and duneD_
    unsigned int BnumBlocks = duneB_.nonzeroes();
    unsigned int DnumBlocks = duneD_.nonzeroes();

    // duneC
    std::vector<unsigned int> Ccols;
    std::vector<double> Cvals;
    Ccols.reserve(BnumBlocks);
    Cvals.reserve(BnumBlocks * numEq * numWellEq);
    for (auto rowC = duneC_.begin(); rowC != duneC_.end(); ++rowC) {
        for (auto colC = rowC->begin(), endC = rowC->end(); colC != endC; ++colC) {
            Ccols.emplace_back(colC.index());
            for (int i = 0; i < numWellEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    Cvals.emplace_back((*colC)[i][j]);
                }
            }
        }
    }

    // duneD
    if constexpr (std::is_same_v<Scalar,float>) {
        OPM_THROW(std::runtime_error, "Cannot use UMFPack with floats");
    } else {
        Dune::UMFPack<DiagMatWell> umfpackMatrix(duneD_, 0);
        double* Dvals = umfpackMatrix.getInternalMatrix().getValues();
        auto* Dcols = umfpackMatrix.getInternalMatrix().getColStart();
        auto* Drows = umfpackMatrix.getInternalMatrix().getRowIndex();

        // duneB
        std::vector<unsigned int> Bcols;
        std::vector<unsigned int> Brows;
        std::vector<double> Bvals;
        Bcols.reserve(BnumBlocks);
        Brows.reserve(Mb+1);
        Bvals.reserve(BnumBlocks * numEq * numWellEq);
        Brows.emplace_back(0);
        unsigned int sumBlocks = 0;
        for (auto rowB = duneB_.begin(); rowB != duneB_.end(); ++rowB) {
            int sizeRow = 0;
            for (auto colB = rowB->begin(), endB = rowB->end(); colB != endB; ++colB) {
                Bcols.emplace_back(colB.index());
                for (int i = 0; i < numWellEq; ++i) {
                    for (int j = 0; j < numEq; ++j) {
                        Bvals.emplace_back((*colB)[i][j]);
                    }
                }
                sizeRow++;
            }
            sumBlocks += sizeRow;
            Brows.emplace_back(sumBlocks);
        }

        wellContribs.addMultisegmentWellContribution(numEq,
                                                     numWellEq,
                                                     Mb,
                                                     Bvals,
                                                     Bcols,
                                                     Brows,
                                                     DnumBlocks,
                                                     Dvals,
                                                     Dcols,
                                                     Drows,
                                                     Cvals);
    }
}
#endif

template<class Scalar, int numWellEq, int numEq>
template<class SparseMatrixAdapter>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
extract(SparseMatrixAdapter& jacobian) const
{
    const auto invDuneD = mswellhelpers::invertWithUMFPack<BVectorWell>(duneD_.M(),
                                                                        numWellEq,
                                                                        *duneDSolver_);

    // We need to change matrix A as follows
    // A -= C^T D^-1 B
    // D is a (nseg x nseg) block matrix with (4 x 4) blocks.
    // B and C are (nseg x ncells) block matrices with (4 x 4 blocks).
    // They have nonzeros at (i, j) only if this well has a
    // perforation at cell j connected to segment i.  The code
    // assumes that no cell is connected to more than one segment,
    // i.e. the columns of B/C have no more than one nonzero.
    for (std::size_t rowC = 0; rowC < duneC_.N(); ++rowC) {
        for (auto colC = duneC_[rowC].begin(),
                  endC = duneC_[rowC].end(); colC != endC; ++colC) {
            const auto row_index = colC.index();
            for (std::size_t rowB = 0; rowB < duneB_.N(); ++rowB) {
                for (auto colB = duneB_[rowB].begin(),
                          endB = duneB_[rowB].end(); colB != endB; ++colB) {
                    const auto col_index = colB.index();
                    OffDiagMatrixBlockWellType tmp1;
                    detail::multMatrixImpl(invDuneD[rowC][rowB], (*colB), tmp1, std::true_type());
                    typename SparseMatrixAdapter::MatrixBlock tmp2;
                    detail::multMatrixTransposedImpl((*colC), tmp1, tmp2, std::false_type());
                    jacobian.addToBlock(row_index, col_index, tmp2);
                }
            }
        }
    }
}

template<class Scalar, int numWellEq, int numEq>
template<class PressureMatrix>
void MultisegmentWellEquations<Scalar,numWellEq,numEq>::
extractCPRPressureMatrix(PressureMatrix& jacobian,
                         const BVector& weights,
                         const int pressureVarIndex,
                         const bool /*use_well_weights*/,
                         const WellInterfaceGeneric<Scalar>& well,
                         const int seg_pressure_var_ind,
                         const WellState<Scalar>& well_state) const
{
    // Add the pressure contribution to the cpr system for the well

    // Add for coupling from well to reservoir
    const int welldof_ind = duneC_.M() + well.indexOfWell();
    if (!well.isPressureControlled(well_state)) {
        for (std::size_t rowC = 0; rowC < duneC_.N(); ++rowC) {
            for (auto colC = duneC_[rowC].begin(),
                      endC = duneC_[rowC].end(); colC != endC; ++colC) {
                const auto row_index = colC.index();
                const auto& bw = weights[row_index];
                double matel = 0.0;

                for (std::size_t i = 0; i< bw.size(); ++i) {
                    matel += bw[i]*(*colC)[seg_pressure_var_ind][i];
                }
                jacobian[row_index][welldof_ind] += matel;
            }
        }
    }

    // make cpr weights for well by pure avarage of reservoir weights of the perforations
    if (!well.isPressureControlled(well_state)) {
        auto well_weight = weights[0];
        well_weight = 0.0;
        int num_perfs = 0;
        for (std::size_t rowB = 0; rowB < duneB_.N(); ++rowB) {
            for (auto colB = duneB_[rowB].begin(),
                      endB = duneB_[rowB].end(); colB != endB; ++colB) {
                const auto col_index = colB.index();
                const auto& bw = weights[col_index];
                well_weight += bw;
                num_perfs += 1;
            }
        }

        well_weight /= num_perfs;
        assert(num_perfs > 0);

        // Add for coupling from reservoir to well and caclulate diag elelement corresping to incompressible standard well
        double diag_ell = 0.0;
        for (std::size_t rowB = 0; rowB < duneB_.N(); ++rowB) {
            const auto& bw = well_weight;
            for (auto colB = duneB_[rowB].begin(),
                      endB = duneB_[rowB].end(); colB != endB; ++colB) {
                const auto col_index = colB.index();
                double matel = 0.0;
                for (std::size_t i = 0; i< bw.size(); ++i) {
                    matel += bw[i] *(*colB)[i][pressureVarIndex];
                }
                jacobian[welldof_ind][col_index] += matel;
                diag_ell -= matel;
            }
        }

#define EXTRA_DEBUG_MSW 0
#if EXTRA_DEBUG_MSW
        if (diag_ell <= 0.0) {
            std::stringstream msg;
            msg << "Diagonal element for cprw on "
                      << this->name()
                      << " is " << diag_ell;
            OpmLog::warning(msg.str());
        }
#endif
#undef EXTRA_DEBUG_MSW
        jacobian[welldof_ind][welldof_ind] = diag_ell;
    } else {
        jacobian[welldof_ind][welldof_ind] = 1.0; // maybe we could have used diag_ell if calculated
    }
}

#define INSTANTIATE(T, numWellEq, numEq)                                                       \
    template class MultisegmentWellEquations<T,numWellEq,numEq>;                               \
    template void MultisegmentWellEquations<T,numWellEq,numEq>::                               \
        extract(Linear::IstlSparseMatrixAdapter<MatrixBlock<T,numEq,numEq>>&) const;           \
    template void MultisegmentWellEquations<T,numWellEq,numEq>::                               \
        extractCPRPressureMatrix(Dune::BCRSMatrix<MatrixBlock<T,1,1>>&,                        \
                                 const MultisegmentWellEquations<T,numWellEq,numEq>::BVector&, \
                                 const int,                                                    \
                                 const bool,                                                   \
                                 const WellInterfaceGeneric<T>&,                               \
                                 const int,                                                    \
                                 const WellState<T>&) const;

#define INSTANTIATE_TYPE(T) \
    INSTANTIATE(T,2,1)      \
    INSTANTIATE(T,2,2)      \
    INSTANTIATE(T,2,6)      \
    INSTANTIATE(T,3,2)      \
    INSTANTIATE(T,3,3)      \
    INSTANTIATE(T,3,4)      \
    INSTANTIATE(T,4,3)      \
    INSTANTIATE(T,4,4)      \
    INSTANTIATE(T,4,5)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
