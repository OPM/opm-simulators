/*
  Copyright 2022 Equinor ASA

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

#include <numeric>
#include <iterator>
#include <algorithm>

#include <dune/common/timer.hh>
#include <dune/istl/scaledidmatrix.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/BSPAI.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
BSPAI<block_size>::BSPAI(ILUReorder opencl_ilu_reorder_, int verbosity_) :
    Preconditioner<block_size>(verbosity_)
{
    assert(opencl_ilu_reorder_ == ILUReorder::NONE);
}

template <unsigned int block_size>
BSPAI<block_size>::~BSPAI()
{
    solver.free();
}

template <unsigned int block_size>
void BSPAI<block_size>::setOpencl(std::shared_ptr<cl::Context>& context_, std::shared_ptr<cl::CommandQueue>& queue_)
{
    context = context_;
    queue = queue_;
}

template <unsigned int block_size>
void BSPAI<block_size>::buildIJSets(int tcol)
{
    jset.clear();

    auto fcol = colIndices.begin() + rowPointers[tcol];
    auto lcol = colIndices.begin() + rowPointers[tcol + 1];
    jset.insert(fcol, lcol);

    for(int f = 0; f <= fill_in; f++){
        iset.clear();

        for(auto it = jset.begin(); it != jset.end(); ++it){
            auto frow = colIndices.begin() + rowPointers[*it];
            auto lrow = colIndices.begin() + rowPointers[*it + 1];
            iset.insert(frow, lrow);
        }

        if(f < fill_in){
            jset = iset;
        }
    }
}

template <unsigned int block_size>
void BSPAI<block_size>::gatherSubmatIndices()
{
    std::vector<int> tmp(submatIndices);
    std::transform(tmp.begin(), tmp.end(), submatIndices.begin(),
        [=](int i){return std::distance(jset.begin(), std::find(jset.begin(), jset.end(), i));});
}

template <unsigned int block_size>
bool BSPAI<block_size>::analyze_matrix(BlockedMatrix *mat)
{
    const unsigned int bs = block_size;
    this->Nb = mat->Nb;
    this->N = mat->Nb * bs;
    this->nnzb = mat->nnzbs;
    this->nnz = mat->nnzbs * bs * bs;

    submat.resize(Nb);
    eyeBlockIndices.resize(Nb);
    submatValsPositions.resize(Nb);
    rowPointers.resize(Nb + 1);
    spaiColPointers.resize(Nb + 1);
    colIndices.resize(nnzb);

    std::copy(mat->rowPointers, mat->rowPointers + Nb + 1, rowPointers.begin());
    std::copy(mat->colIndices, mat->colIndices + nnzb, colIndices.begin());

    Dune::Timer t_analyze_matrix;

    for(int tcol = 0; tcol < Nb; tcol++){
        buildIJSets(tcol);

        submatIndices.clear();
        submatPointers.assign(iset.size() + 1, 0);

        unsigned int i = 1;
        for(auto rit = iset.begin(); rit != iset.end(); ++rit){
            auto fcol = colIndices.begin() + rowPointers[*rit];
            auto lcol = colIndices.begin() + rowPointers[*rit + 1];

            for(auto cit = fcol; cit != lcol; ++cit){
                if(jset.count(*cit)){
                    submatIndices.push_back(*cit);
                    submatValsPositions[tcol].resize(submatValsPositions[tcol].size() + bs * bs);
                    std::iota(submatValsPositions[tcol].end() - bs * bs, submatValsPositions[tcol].end(), (rowPointers[*rit] + cit - fcol) * bs * bs);
                }
            }

            submatPointers[i] = submatIndices.size();
            i++;
        }

        submat[tcol].setSize(iset.size(), jset.size(), submatPointers.back());
        submat[tcol].setBuildMode(DuneMat::row_wise);

        gatherSubmatIndices();

        for(typename DuneMat::CreateIterator row = submat[tcol].createbegin(); row != submat[tcol].createend(); ++row){
            for(int p = submatPointers[row.index()]; p < submatPointers[row.index() + 1]; ++p){
                row.insert(submatIndices[p]);
            }
        }

        spaiColPointers[tcol + 1] = spaiColPointers[tcol] + jset.size();
        spaiRowIndices.insert(spaiRowIndices.end(), jset.begin(), jset.end());
        eyeBlockIndices[tcol] = std::distance(iset.begin(), std::find(iset.begin(), iset.end(), tcol));
    }

    spaiNnzValues.resize(spaiColPointers.back() * bs * bs);

    if(verbosity >= 4){
        std::ostringstream out;
        out << "BSPAI::analyze_matrix time: " << t_analyze_matrix.stop() << " s";
        OpmLog::info(out.str());
    }

    std::call_once(ocl_init, [&]() {
        d_spaiColPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * spaiColPointers.size());
        d_spaiRowIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * spaiRowIndices.size());
        d_spaiNnzValues = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * spaiNnzValues.size());

        events.resize(2);
        err = queue->enqueueWriteBuffer(d_spaiColPointers, CL_FALSE, 0, sizeof(int) * spaiColPointers.size(), spaiColPointers.data(), nullptr, &events[0]);
        err |= queue->enqueueWriteBuffer(d_spaiRowIndices, CL_FALSE, 0, sizeof(int) * spaiRowIndices.size(), spaiRowIndices.data(), nullptr, &events[1]);

        cl::WaitForEvents(events);
        events.clear();

        if (err != CL_SUCCESS) {
            OPM_THROW(std::logic_error, "BSPAI::analyze_matrix OpenCL enqueueWriteBuffer error");
        }
    });

    return true;
}

template <unsigned int block_size>
bool BSPAI<block_size>::create_preconditioner(BlockedMatrix *mat)
{
    const unsigned int bs = block_size;
    int count;

    if (bs != 3) {
        OPM_THROW(std::logic_error, "Currently, creation of BSPAI preconditioner on GPU only supports block_size = 3");
    }

    Dune::Timer t_preconditioner;

    solver.setBlocked(bs > 1);

    for(int tcol = 0; tcol < Nb; tcol++){
        count = 0;
        for(auto row = submat[tcol].begin(); row != submat[tcol].end(); ++row){
            for(auto col = (*row).begin(); col != (*row).end(); ++col){
                for(auto br = (*col).begin(); br != (*col).end(); ++br){
                    for(auto bc = (*br).begin(); bc != (*br).end(); ++bc){
                        (*bc) = mat->nnzValues[submatValsPositions[tcol][count]];
                        ++count;
                    }
                }
            }
        }

        sol.resize(submat[tcol].M());
        rhs.resize(submat[tcol].N());
        rhs = 0;
        rhs[eyeBlockIndices[tcol]] = Dune::ScaledIdentityMatrix<double, bs>(1);

        solver.setMatrix(submat[tcol]);
        solver.apply(sol, rhs, res);

        for(unsigned int i = 0; i < submat[tcol].M(); i++){
            for(unsigned int j = 0; j < bs; j++){
                for(unsigned int k = 0; k < bs; k++){
                    spaiNnzValues[(spaiColPointers[tcol] + i) * bs * bs + j * bs + k] = sol[i][j][k];
                }
            }
        }
    }

    events.resize(1);
    err = queue->enqueueWriteBuffer(d_spaiNnzValues, CL_FALSE, 0, sizeof(double) * spaiNnzValues.size(), spaiNnzValues.data(), nullptr, &events[0]);

    cl::WaitForEvents(events);
    events.clear();

    if (err != CL_SUCCESS) {
        OPM_THROW(std::logic_error, "BSPAI::create_preconditioner OpenCL enqueueWriteBuffer error");
    }

    return true;
}

template <unsigned int block_size>
void BSPAI<block_size>::apply(const cl::Buffer& d_x, cl::Buffer& d_y){
    const unsigned int bs = block_size;

    events.resize(1);
    err = queue->enqueueFillBuffer(d_y, 0, 0, sizeof(double) * N, nullptr, &events[0]);

    cl::WaitForEvents(events);
    events.clear();

    if (err != CL_SUCCESS) {
        OPM_THROW(std::logic_error, "BSPAI::apply OpenCL enqueueWriteBuffer error");
    }

    OpenclKernels::csc_spmv_blocked(d_spaiNnzValues, d_spaiColPointers, d_spaiRowIndices, d_x, d_y, Nb, bs);
}

#define INSTANTIATE_BDA_FUNCTIONS(n)  \
template class BSPAI<n>;

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

}
}
