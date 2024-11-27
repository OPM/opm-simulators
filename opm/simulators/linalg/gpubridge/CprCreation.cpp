/*
  Copyright 2024 Equinor ASA

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
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <dune/common/shared_ptr.hh>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/simulators/linalg/gpubridge/GpuBridge.hpp>
#include <opm/simulators/linalg/gpubridge/BlockedMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/CprCreation.hpp>

#include <opm/simulators/linalg/gpubridge/Misc.hpp>

namespace Opm::Accelerator {

using Opm::OpmLog;
using Dune::Timer;

template <class Scalar, unsigned int block_size>
CprCreation<Scalar, block_size>::CprCreation()
{
    diagIndices.resize(1);
}

template <class Scalar, unsigned int block_size>
void CprCreation<Scalar, block_size>::
create_preconditioner_amg(BlockedMatrix<Scalar> *mat_)
{
    mat = mat_;
    cprNb = mat_->Nb;
    cprnnzb = mat_->nnzbs;
    cprN = cprNb * block_size;
    cprnnz = cprnnzb * block_size * block_size;

    coarse_vals.resize(cprnnzb);
    coarse_x.resize(cprNb);
    coarse_y.resize(cprNb);
    weights.resize(cprN);

    try{
        Scalar rhs[] = {0, 0, 0};
        rhs[pressure_idx] = 1;

        // find diagonal index for each row
        if (diagIndices[0].empty()) {
            diagIndices[0].resize(cprNb);
            for (int row = 0; row < cprNb; ++row) {
                int start = mat->rowPointers[row];
                int end = mat->rowPointers[row + 1];
                auto candidate = std::find(mat->colIndices + start, mat->colIndices + end, row);
                assert(candidate != mat->colIndices + end);
                diagIndices[0][row] = candidate - mat->colIndices;
            }
        }

        // calculate weights for each row
        for (int row = 0; row < cprNb; ++row) {
            // solve to find weights
            Scalar *row_weights = weights.data() + block_size * row; // weights for this row
            solve_transposed_3x3(mat->nnzValues + block_size * block_size * diagIndices[0][row], rhs, row_weights);

            // normalize weights for this row
            Scalar abs_max = get_absmax(row_weights, block_size);
            for(unsigned int i = 0; i < block_size; i++){
                row_weights[i] /= abs_max;
            }
        }

        // extract pressure
        // transform blocks to scalars to create scalar linear system
        for (int row = 0; row < cprNb; ++row) {
            int start = mat->rowPointers[row];
            int end = mat->rowPointers[row + 1];
            for (int idx = start; idx < end; ++idx) {
                Scalar *block = mat->nnzValues + idx * block_size * block_size;
                Scalar *row_weights = weights.data() + block_size * row;
                Scalar value = 0.0;
                for (unsigned int i = 0; i < block_size; ++i) {
                    value += block[block_size * i + pressure_idx] * row_weights[i];
                }
                coarse_vals[idx] = value;
            }
        }

#if HAVE_MPI
        using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
        using Communication = Dune::Amg::SequentialInformation;
#endif
        using OverlapFlags = Dune::NegateSet<Communication::OwnerSet>;
        if (recalculate_aggregates) {
            dune_coarse = std::make_unique<DuneMat>(cprNb, cprNb, cprnnzb, DuneMat::row_wise);

//             typedef DuneMat::CreateIterator Iter;
            using Iter = typename DuneMat::CreateIterator;

            // setup sparsity pattern
            for(Iter row = dune_coarse->createbegin(); row != dune_coarse->createend(); ++row){
                int start = mat->rowPointers[row.index()];
                int end = mat->rowPointers[row.index() + 1];
                for (int idx = start; idx < end; ++idx) {
                    int col = mat->colIndices[idx];
                    row.insert(col);
                }
            }

            // set values
            for (int row = 0; row < cprNb; ++row) {
                int start = mat->rowPointers[row];
                int end = mat->rowPointers[row + 1];
                for (int idx = start; idx < end; ++idx) {
                    int col = mat->colIndices[idx];
                    (*dune_coarse)[row][col] = coarse_vals[idx];
                }
            }

            dune_op = std::make_shared<MatrixOperator>(*dune_coarse);
            Dune::Amg::SequentialInformation seqinfo;
            dune_amg = std::make_unique<DuneAmg>(dune_op, Dune::stackobject_to_shared_ptr(seqinfo));

            Opm::PropertyTree property_tree;
            property_tree.put("alpha", 0.333333333333);

            // The matrix has a symmetric sparsity pattern, but the values are not symmetric
            // Yet a SymmetricDependency is used in AMGCPR
            // An UnSymmetricCriterion is also available
            // using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::UnSymmetricCriterion<DuneMat, Dune::Amg::FirstDiagonal> >;
            using CriterionBase = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<DuneMat, Dune::Amg::FirstDiagonal>>;
            using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;
            const Criterion c = Opm::AMGHelper<MatrixOperator,Dune::Amg::SequentialInformation,DuneMat,DuneVec>::criterion(property_tree);
            num_pre_smooth_steps = c.getNoPreSmoothSteps();
            num_post_smooth_steps = c.getNoPostSmoothSteps();

            dune_amg->template build<OverlapFlags>(c);

            analyzeHierarchy();
            analyzeAggregateMaps();

            recalculate_aggregates = false;
        } else {
            // update values of coarsest level in AMG
            // this works because that level is actually a reference to the DuneMat held by dune_coarse
            for (int row = 0; row < cprNb; ++row) {
                int start = mat->rowPointers[row];
                int end = mat->rowPointers[row + 1];
                for (int idx = start; idx < end; ++idx) {
                    int col = mat->colIndices[idx];
                    (*dune_coarse)[row][col] = coarse_vals[idx];
                }
            }

            // update the rest of the AMG hierarchy
            dune_amg->recalculateGalerkin(OverlapFlags());
            analyzeHierarchy();
        }

    } catch (const std::exception& ex) {
        std::cerr << "Caught exception: " << ex.what() << std::endl;
        throw ex;
    }
}

template <class Scalar, unsigned int block_size>
void CprCreation<Scalar, block_size>::
analyzeHierarchy()
{
    const typename DuneAmg::ParallelMatrixHierarchy& matrixHierarchy = dune_amg->matrices();

    // store coarsest AMG level in umfpack format, also performs LU decomposition
    if constexpr (std::is_same_v<Scalar,float>) {
        OPM_THROW(std::runtime_error, "Cannot use CPR with float Scalar due to UMFPACK");
    } else {
        umfpack.setMatrix((*matrixHierarchy.coarsest()).getmat());
    }

    num_levels = dune_amg->levels();
    level_sizes.resize(num_levels);
    diagIndices.resize(num_levels);

    Amatrices.reserve(num_levels);
    Rmatrices.reserve(num_levels - 1);  // coarsest level does not need one
    invDiags.reserve(num_levels);

    Amatrices.clear();
    invDiags.clear();

    // matrixIter.dereference() returns MatrixAdapter
    // matrixIter.dereference().getmat() returns BCRSMat
    typename DuneAmg::ParallelMatrixHierarchy::ConstIterator matrixIter = matrixHierarchy.finest();
    for(int level = 0; level < num_levels; ++matrixIter, ++level) {
        const auto& A = matrixIter.dereference().getmat();
        level_sizes[level] = A.N();
        diagIndices[level].reserve(A.N());

        // extract matrix A
        Amatrices.emplace_back(A.N(), A.nonzeroes());
        // contiguous copy is not possible
        // std::copy(&(A[0][0][0][0]), &(A[0][0][0][0]) + A.nonzeroes(), Amatrices.back().nnzValues.data());
        // also update diagonal indices if needed, level 0 is already filled in create_preconditioner()
        int nnz_idx = 0;
        const bool fillDiagIndices = diagIndices[level].empty();
        for (typename DuneMat::const_iterator r = A.begin(); r != A.end(); ++r) {
            for (auto c = r->begin(); c != r->end(); ++c) {
                Amatrices.back().nnzValues[nnz_idx] = A[r.index()][c.index()];
                if (fillDiagIndices && r.index() == c.index()) {
                    diagIndices[level].emplace_back(nnz_idx);
                }
                nnz_idx++;
            }
        }

        Opm::GpuBridge<DuneMat, DuneVec, 1>::copySparsityPatternFromISTL(A, Amatrices.back().rowPointers, Amatrices.back().colIndices);

        // compute inverse diagonal values for current level
        invDiags.emplace_back(A.N());
        for (unsigned int row = 0; row < A.N(); ++row) {
            invDiags.back()[row] = 1 / Amatrices.back().nnzValues[diagIndices[level][row]];
        }
    }
}


template <class Scalar, unsigned int block_size>
void CprCreation<Scalar, block_size>::
analyzeAggregateMaps() 
{
    PcolIndices.resize(num_levels - 1);
    Rmatrices.clear();

    const typename DuneAmg::AggregatesMapList& aggregatesMaps = dune_amg->aggregatesMaps();

    typename DuneAmg::AggregatesMapList::const_iterator mapIter = aggregatesMaps.begin();
    for(int level = 0; level < num_levels - 1; ++mapIter, ++level) {
        typename DuneAmg::AggregatesMap *map = *mapIter;

        Rmatrices.emplace_back(level_sizes[level+1], level_sizes[level], level_sizes[level]);
        std::fill(Rmatrices.back().nnzValues.begin(), Rmatrices.back().nnzValues.end(), 1.0);

        // get indices for each row of P and R
        std::vector<std::vector<unsigned> > indicesR(level_sizes[level+1]);
        PcolIndices[level].resize(level_sizes[level]);

        using AggregateIterator = typename DuneAmg::AggregatesMap::const_iterator;
        for(AggregateIterator ai = map->begin(); ai != map->end(); ++ai){
            if (*ai != DuneAmg::AggregatesMap::ISOLATED) {
                const long int diff = ai - map->begin();
                PcolIndices[level][diff] = *ai;
                indicesR[*ai].emplace_back(diff);
            }
        }

        int col_idx = 0;
        // set sparsity pattern of R
        Rmatrices.back().rowPointers[0] = 0;
        for (unsigned int i = 0; i < indicesR.size(); ++i) {
            Rmatrices.back().rowPointers[i + 1] = Rmatrices.back().rowPointers[i] + indicesR[i].size();
            for (auto it = indicesR[i].begin(); it != indicesR[i].end(); ++it) {
                Rmatrices.back().colIndices[col_idx++] = *it;
            }
        }
    }
}

#define INSTANTIATE_TYPE(T)          \
    template class CprCreation<T,1>; \
    template class CprCreation<T,2>; \
    template class CprCreation<T,3>; \
    template class CprCreation<T,4>; \
    template class CprCreation<T,5>; \
    template class CprCreation<T,6>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm


