/*
  Copyright 2021 Equinor ASA

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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <dune/istl/umfpack.hh>
#include <dune/common/shared_ptr.hh>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/CPR.hpp>


namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
CPR<block_size>::CPR(int verbosity_, ILUReorder opencl_ilu_reorder_) :
    verbosity(verbosity_), opencl_ilu_reorder(opencl_ilu_reorder_)
{}


template <unsigned int block_size>
void CPR<block_size>::init(int Nb_, int nnzb_, std::shared_ptr<cl::Context>& context_, std::shared_ptr<cl::CommandQueue>& queue_)
{
    this->Nb = Nb_;
    this->nnzb = nnzb_;
    this->N = Nb_ * block_size;
    this->nnz = nnzb_ * block_size * block_size;

    context = context_;
    queue = queue_;

    coarse_vals.resize(nnzb);
    coarse_x.resize(Nb);
    coarse_y.resize(Nb);
    weights.resize(N);
    diagIndices.resize(1);
}


double get_absmax(const double *data, const int B) {
    double abs_max = 0.0;
    for (int i = 0; i < B; ++i) {
        if (std::fabs(data[i]) > abs_max) {
            abs_max = std::fabs(data[i]);
        }
    }
    return abs_max;
}


// solve A^T * x = b
void solve_transposed_3x3(const double *A, const double *b, double *x) {
    const int B = 3;
    // from dune-common/densematrix.hh, but transposed, so replace [r*B+c] with [r+c*B]
    double t4  = A[0+0*B] * A[1+1*B];
    double t6  = A[0+0*B] * A[1+2*B];
    double t8  = A[0+1*B] * A[1+0*B];
    double t10 = A[0+2*B] * A[1+0*B];
    double t12 = A[0+1*B] * A[2+0*B];
    double t14 = A[0+2*B] * A[2+0*B];

    double d = (t4*A[2+2*B]-t6*A[2+1*B]-t8*A[2+2*B]+
          t10*A[2+1*B]+t12*A[1+2*B]-t14*A[1+1*B]); //determinant

    x[0] = (b[0]*A[1+1*B]*A[2+2*B] - b[0]*A[2+1*B]*A[1+2*B]
          - b[1] *A[0+1*B]*A[2+2*B] + b[1]*A[2+1*B]*A[0+2*B]
          + b[2] *A[0+1*B]*A[1+2*B] - b[2]*A[1+1*B]*A[0+2*B]) / d;

    x[1] = (A[0+0*B]*b[1]*A[2+2*B] - A[0+0*B]*b[2]*A[1+2*B]
          - A[1+0*B] *b[0]*A[2+2*B] + A[1+0*B]*b[2]*A[0+2*B]
          + A[2+0*B] *b[0]*A[1+2*B] - A[2+0*B]*b[1]*A[0+2*B]) / d;

    x[2] = (A[0+0*B]*A[1+1*B]*b[2] - A[0+0*B]*A[2+1*B]*b[1]
          - A[1+0*B] *A[0+1*B]*b[2] + A[1+0*B]*A[2+1*B]*b[0]
          + A[2+0*B] *A[0+1*B]*b[1] - A[2+0*B]*A[1+1*B]*b[0]) / d;
}


template <unsigned int block_size>
void CPR<block_size>::create_preconditioner(BlockedMatrix *mat_) {
    this->mat = mat_;

    try{
        double rhs[] = {0, 0, 0};
        rhs[pressure_idx] = 1;

        // find diagonal index for each row
        if (diagIndices[0].empty()) {
            diagIndices[0].resize(Nb);
            for (int row = 0; row < Nb; ++row) {
                int start = mat->rowPointers[row];
                int end = mat->rowPointers[row + 1];
                auto candidate = std::find(mat->colIndices + start, mat->colIndices + end, row);
                assert(candidate != mat->colIndices + end);
                diagIndices[0][row] = candidate - mat->colIndices;
            }
        }

        // calculate weights for each row
        for (int row = 0; row < Nb; ++row) {
            // solve to find weights
            double *row_weights = weights.data() + block_size * row; // weights for this row
            solve_transposed_3x3(mat->nnzValues + block_size * block_size * diagIndices[0][row], rhs, row_weights);

            // normalize weights for this row
            double abs_max = get_absmax(row_weights, block_size);
            for(unsigned int i = 0; i < block_size; i++){
                row_weights[i] /= abs_max;
            }
        }

        // transform blocks to scalars
        for (int row = 0; row < Nb; ++row) {
            int start = mat->rowPointers[row];
            int end = mat->rowPointers[row + 1];
            for (int idx = start; idx < end; ++idx) {
                double *block = mat->nnzValues + idx * block_size * block_size;
                double *row_weights = weights.data() + block_size * row;
                double value = 0.0;
                for (unsigned int i = 0; i < block_size; ++i) {
                    value += block[block_size * i + pressure_idx] * row_weights[i];
                }
                coarse_vals[idx] = value;
            }
        }

        using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
        using OverlapFlags = Dune::NegateSet<Communication::OwnerSet>;
        if (recalculate_aggregates) {
            dune_coarse = std::make_unique<DuneMat>(Nb, Nb, nnzb, DuneMat::row_wise);

            typedef DuneMat::CreateIterator Iter;

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
            for (int row = 0; row < Nb; ++row) {
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
            const Criterion c = Opm::PreconditionerFactory<MatrixOperator, Dune::Amg::SequentialInformation>::amgCriterion(property_tree);

            dune_amg->build<OverlapFlags>(c);

            analyzeHierarchy();
            analyzeAggregateMaps();

            recalculate_aggregates = false;
        } else {
            // update values of coarsest level in AMG
            // this works because that level is actually a reference to the DuneMat held by dune_coarse
            for (int row = 0; row < Nb; ++row) {
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

        // initialize OpenclMatrices and Buffers if needed
        std::call_once(opencl_buffers_allocated, [&](){
            d_Amatrices.reserve(num_levels);
            d_Pmatrices.reserve(num_levels - 1);
            d_Rmatrices.reserve(num_levels - 1);
            d_invDiags.reserve(num_levels-1);
            for (Matrix& m : Amatrices) {
                d_Amatrices.emplace_back(context.get(), m.N, m.M, m.nnzs);
            }
            for (Matrix& m : Pmatrices) {
                d_Pmatrices.emplace_back(context.get(), m.N, m.M, m.nnzs);
                d_invDiags.emplace_back(*context, CL_MEM_READ_WRITE, sizeof(double) * m.N); // create a cl::Buffer
                d_t.emplace_back(*context, CL_MEM_READ_WRITE, sizeof(double) * m.N);
            }
            for (Matrix& m : Rmatrices) {
                d_Rmatrices.emplace_back(context.get(), m.N, m.M, m.nnzs);
                d_f.emplace_back(*context, CL_MEM_READ_WRITE, sizeof(double) * m.N);
                d_u.emplace_back(*context, CL_MEM_READ_WRITE, sizeof(double) * m.N);
            }
            d_weights = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
            d_rs = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
            d_mat = std::make_unique<OpenclMatrix<block_size> >(context.get(), Nb, Nb, nnzb);
            d_coarse_y = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * Nb);
            d_coarse_x = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * Nb);
        });

        // upload matrices and vectors to GPU
        d_mat->upload(queue.get(), mat);

        err = CL_SUCCESS;
        events.resize(Pmatrices.size() + 1);
        for (unsigned int i = 0; i < Pmatrices.size(); ++i) {
            d_Amatrices[i].upload(queue.get(), &Amatrices[i]);

            err |= queue->enqueueWriteBuffer(d_invDiags[i], CL_FALSE, 0, sizeof(double) * Amatrices[i].N, invDiags[i].data(), nullptr, &events[i]);
        }
        err |= queue->enqueueWriteBuffer(*d_weights, CL_FALSE, 0, sizeof(double) * N, weights.data(), nullptr, &events[Pmatrices.size()]);
        cl::WaitForEvents(events);
        events.clear();
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "CPR OpenCL enqueueWriteBuffer error");
        }
        for (unsigned int i = 0; i < Pmatrices.size(); ++i) {
            d_Pmatrices[i].upload(queue.get(), &Pmatrices[i]);
        }
        for (unsigned int i = 0; i < Rmatrices.size(); ++i) {
            d_Rmatrices[i].upload(queue.get(), &Rmatrices[i]);
        }

    } catch (const std::exception& ex) {
        std::cerr << "Caught exception: " << ex.what() << std::endl;
        throw ex;
    }
}


template <unsigned int block_size>
void CPR<block_size>::analyzeHierarchy() {
    const DuneAmg::ParallelMatrixHierarchy& matrixHierarchy = dune_amg->matrices();

    num_levels = dune_amg->levels();
    level_sizes.resize(num_levels);
    diagIndices.resize(num_levels);

    Amatrices.reserve(num_levels);
    Pmatrices.reserve(num_levels - 1);  // coarsest level does not need one
    Rmatrices.reserve(num_levels - 1);
    invDiags.reserve(num_levels);

    Amatrices.clear();
    invDiags.clear();

    // matrixIter.dereference() returns MatrixAdapter
    // matrixIter.dereference().getmat() returns BCRSMat
    DuneAmg::ParallelMatrixHierarchy::ConstIterator matrixIter = matrixHierarchy.finest();
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

        Opm::BdaBridge<DuneMat, DuneVec, 1>::getSparsityPattern(A, Amatrices.back().rowPointers, Amatrices.back().colIndices);

        // compute inverse diagonal values for current level
        invDiags.emplace_back(A.N());
        for (unsigned int row = 0; row < A.N(); ++row) {
            invDiags.back()[row] = 1 / Amatrices.back().nnzValues[diagIndices[level][row]];
        }
    }
}


template <unsigned int block_size>
void CPR<block_size>::analyzeAggregateMaps() {

    Pmatrices.clear();
    Rmatrices.clear();

    const DuneAmg::AggregatesMapList& aggregatesMaps = dune_amg->aggregatesMaps();

    DuneAmg::AggregatesMapList::const_iterator mapIter = aggregatesMaps.begin();
    for(int level = 0; level < num_levels - 1; ++mapIter, ++level) {
        DuneAmg::AggregatesMap *map = *mapIter;

        // get indices for each row of P and R
        std::vector<std::set<long int> > indicesP(level_sizes[level]);
        std::vector<std::set<long int> > indicesR(level_sizes[level+1]);
        using AggregateIterator = DuneAmg::AggregatesMap::const_iterator;
        for(AggregateIterator ai = map->begin(); ai != map->end(); ++ai){
            if (*ai != DuneAmg::AggregatesMap::ISOLATED) {
                const long int diff = ai - map->begin();
                indicesP[diff].insert(*ai);
                indicesR[*ai].insert(diff);
            }
        }

        Pmatrices.emplace_back(level_sizes[level], level_sizes[level+1], level_sizes[level]);
        Rmatrices.emplace_back(level_sizes[level+1], level_sizes[level], level_sizes[level]);
        std::fill(Pmatrices.back().nnzValues.begin(), Pmatrices.back().nnzValues.end(), 1.0);
        std::fill(Rmatrices.back().nnzValues.begin(), Rmatrices.back().nnzValues.end(), 1.0);

        // set sparsity pattern of P
        int col_idx = 0;
        Pmatrices.back().rowPointers[0] = 0;
        for (unsigned int i = 0; i < indicesP.size(); ++i) {
            Pmatrices.back().rowPointers[i + 1] = Pmatrices.back().rowPointers[i] + indicesP[i].size();
            for (auto it = indicesP[i].begin(); it != indicesP[i].end(); ++it) {
                Pmatrices.back().colIndices[col_idx] = *it;
                col_idx++;
            }
        }
        // set sparsity pattern of R
        col_idx = 0;
        Rmatrices.back().rowPointers[0] = 0;
        for (unsigned int i = 0; i < indicesR.size(); ++i) {
            Rmatrices.back().rowPointers[i + 1] = Rmatrices.back().rowPointers[i] + indicesR[i].size();
            for (auto it = indicesR[i].begin(); it != indicesR[i].end(); ++it) {
                Rmatrices.back().colIndices[col_idx] = *it;
                col_idx++;
            }
        }
    }

}


void solve_coarse_umfpack(const Matrix *A, std::vector<double> &y, std::vector<double> &x) {
    const int N = A->N;
    const int nnzs = A->nnzs;

    using DuneMat = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1> >;
    DuneMat B(N, N, nnzs, DuneMat::row_wise);

    typedef DuneMat::CreateIterator Iter;

    // setup sparsity pattern
    for(Iter row = B.createbegin(); row != B.createend(); ++row){
        int start = A->rowPointers[row.index()];
        int end = A->rowPointers[row.index() + 1];
        for (int idx = start; idx < end; ++idx) {
            int col = A->colIndices[idx];
            row.insert(col);
        }
    }

    // set values
    for (int row = 0; row < N; ++row) {
        int start = A->rowPointers[row];
        int end = A->rowPointers[row + 1];
        for (int idx = start; idx < end; ++idx) {
            int col = A->colIndices[idx];
            B[row][col] = A->nnzValues[idx];
        }
    }

    // create umfpack object
    Dune::UMFPack<DuneMat> umfpack(B, 0);

    umfpack.apply(x.data(), y.data());
}


template <unsigned int block_size>
void CPR<block_size>::amg_cycle_gpu(const int level, cl::Buffer &y, cl::Buffer &x) {
    OpenclMatrix<1> *A = &d_Amatrices[level];
    OpenclMatrix<1> *P = &d_Pmatrices[level];
    OpenclMatrix<1> *R = &d_Rmatrices[level];
    int Ncur = A->Nb;

    if (level == num_levels - 1) {
        // solve coarsest level
        std::vector<double> h_y(Ncur), h_x(Ncur, 0);

        events.resize(1);
        err = queue->enqueueReadBuffer(y, CL_FALSE, 0, sizeof(double) * Ncur, h_y.data(), nullptr, &events[0]);
        cl::WaitForEvents(events);
        events.clear();
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "CPR OpenCL enqueueReadBuffer error");
        }

        solve_coarse_umfpack(&Amatrices[level], h_y, h_x);

        events.resize(1);
        err = queue->enqueueWriteBuffer(x, CL_FALSE, 0, sizeof(double) * Ncur, h_x.data(), nullptr, &events[0]);
        cl::WaitForEvents(events);
        events.clear();
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "CPR OpenCL enqueueWriteBuffer error");
        }
        return;
    }
    int Nnext = d_Amatrices[level+1].Nb;

    cl::Buffer& t = d_t[level];
    cl::Buffer& f = d_f[level];
    cl::Buffer& u = d_u[level]; // u was 0-initialized earlier

    // presmooth
    double jacobi_damping = 0.72; // default value in amgcl: 0.72
    OpenclKernels::residual(A->nnzValues, A->colIndices, A->rowPointers, x, y, t, Ncur, 1);
    OpenclKernels::vmul(jacobi_damping, d_invDiags[level], t, x, Ncur);

    // move to coarser level
    OpenclKernels::residual(A->nnzValues, A->colIndices, A->rowPointers, x, y, t, Ncur, 1);
    OpenclKernels::spmv(R->nnzValues, R->colIndices, R->rowPointers, t, f, Nnext, 1, true);
    amg_cycle_gpu(level + 1, f, u);
    OpenclKernels::spmv(P->nnzValues, P->colIndices, P->rowPointers, u, x, Ncur, 1, false);

    // postsmooth
    OpenclKernels::residual(A->nnzValues, A->colIndices, A->rowPointers, x, y, t, Ncur, 1);
    OpenclKernels::vmul(jacobi_damping, d_invDiags[level], t, x, Ncur);
}


// x = prec(y)
template <unsigned int block_size>
void CPR<block_size>::apply(const cl::Buffer& y, cl::Buffer& x) {
    // 0-initialize u and x vectors
    events.resize(d_u.size() + 1);
    err = queue->enqueueFillBuffer(*d_coarse_x, 0, 0, sizeof(double) * Nb, nullptr, &events[0]);
    for (unsigned int i = 0; i < d_u.size(); ++i) {
        err |= queue->enqueueFillBuffer(d_u[i], 0, 0, sizeof(double) * Rmatrices[i].N, nullptr, &events[i + 1]);
    }
    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "CPR OpenCL enqueueWriteBuffer error");
    }

    OpenclKernels::residual(d_mat->nnzValues, d_mat->colIndices, d_mat->rowPointers, x, y, *d_rs, Nb, block_size);
    OpenclKernels::move_to_coarse(*d_rs, *d_weights, *d_coarse_y, Nb);

    amg_cycle_gpu(0, *d_coarse_y, *d_coarse_x);

    OpenclKernels::move_to_fine(*d_coarse_x, x, pressure_idx, Nb);
}



#define INSTANTIATE_BDA_FUNCTIONS(n)  \
template class CPR<n>;

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm


