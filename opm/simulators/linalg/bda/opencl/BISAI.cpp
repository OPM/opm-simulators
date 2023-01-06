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
#include <algorithm>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>
#include <opm/simulators/linalg/bda/opencl/BISAI.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>
#include <opm/simulators/linalg/bda/opencl/ChowPatelIlu.hpp> // disable BISAI if ChowPatel is selected

#include <sstream>

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
BISAI<block_size>::BISAI(bool opencl_ilu_parallel_, int verbosity_) :
    Preconditioner<block_size>(verbosity_)
{
#if CHOW_PATEL
    OPM_THROW(std::logic_error, "Error --linear-solver=isai cannot be used if ChowPatelIlu is used, probably defined by CMake\n");
#endif
    bilu0 = std::make_unique<BILU0<block_size> >(opencl_ilu_parallel_, verbosity_);
}

template <unsigned int block_size>
void BISAI<block_size>::setOpencl(std::shared_ptr<cl::Context>& context_, std::shared_ptr<cl::CommandQueue>& queue_)
{
    context = context_;
    queue = queue_;

    bilu0->setOpencl(context, queue);
}

std::vector<int> buildCsrToCscOffsetMap(std::vector<int> colPointers, std::vector<int> rowIndices){
    std::vector<int> aux(colPointers); // colPointers must be copied to this vector
    std::vector<int> csrToCscOffsetMap(rowIndices.size()); // map must have the same size as the indices vector

    for(unsigned int row = 0; row < colPointers.size() - 1; row++){
        for(int jj = colPointers[row]; jj < colPointers[row+1]; jj++){
            int col = rowIndices[jj];
            int dest = aux[col];
            csrToCscOffsetMap[dest] = jj;
            aux[col]++;
        }
    }

    return csrToCscOffsetMap;
}

template <unsigned int block_size>
bool BISAI<block_size>::analyze_matrix(BlockedMatrix *mat)
{
    return analyze_matrix(mat, nullptr);
}

template <unsigned int block_size>
bool BISAI<block_size>::analyze_matrix(BlockedMatrix *mat, BlockedMatrix *jacMat)
{
    const unsigned int bs = block_size;
    auto *m = mat;

    if (jacMat) {
        m = jacMat;
    }

    this->N = m->Nb * bs;
    this->Nb = m->Nb;
    this->nnz = m->nnzbs * bs * bs;
    this->nnzb = m->nnzbs;

    if (jacMat) {
        return bilu0->analyze_matrix(mat, jacMat);
    } else {
        return bilu0->analyze_matrix(mat);
    }
}

template <unsigned int block_size>
void BISAI<block_size>::buildLowerSubsystemsStructures(){
    lower.subsystemPointers.assign(Nb + 1, 0);

    Dune::Timer t_buildLowerSubsystemsStructures;

    for(int tcol = 0; tcol < Nb; tcol++){
        int frow = diagIndex[tcol] + 1;
        int lrow = colPointers[tcol + 1];
        int nx = lrow - frow;
        int nv = 0;

        for(int sweep = 0; sweep < nx - 1; sweep++){
            for(int xid = sweep + 1; xid < nx; xid++){
                for(int ptr = diagIndex[rowIndices[frow + sweep]] + 1; ptr < colPointers[rowIndices[frow + sweep + 1]]; ptr++){
                    if(rowIndices[ptr] == rowIndices[frow + xid]){
                        lower.nzIndices.push_back(csrToCscOffsetMap[ptr]);
                        lower.knownRhsIndices.push_back(csrToCscOffsetMap[frow + sweep]);
                        lower.unknownRhsIndices.push_back(csrToCscOffsetMap[frow + xid]);
                        nv++;
                    }
                }
            }
        }

        lower.subsystemPointers[tcol + 1] = lower.subsystemPointers[tcol] + nv;
    }

    if(verbosity >= 4){
        std::ostringstream out;
        out << "BISAI buildLowerSubsystemsStructures time: " << t_buildLowerSubsystemsStructures.stop() << " s";
        OpmLog::info(out.str());
    }
}

template <unsigned int block_size>
void BISAI<block_size>::buildUpperSubsystemsStructures(){
    upper.subsystemPointers.assign(Nb + 1, 0);

    Dune::Timer t_buildUpperSubsystemsStructures;

    for(int tcol = 0; tcol < Nb; tcol++){
        int frow = colPointers[tcol];
        int lrow = diagIndex[tcol];
        int nx = lrow - frow + 1;
        int nv = 0;

        for(int sweep = 0; sweep < nx - 1; sweep++){
            for(int xid = 0; xid < nx; xid++){
                for(int ptr = colPointers[rowIndices[lrow - sweep]]; ptr < diagIndex[rowIndices[lrow - sweep]]; ptr++){
                    if(rowIndices[ptr] == rowIndices[lrow - xid]){
                        upper.nzIndices.push_back(csrToCscOffsetMap[ptr]);
                        upper.knownRhsIndices.push_back(csrToCscOffsetMap[lrow - sweep]);
                        upper.unknownRhsIndices.push_back(csrToCscOffsetMap[lrow - xid]);
                        nv++;
                    }
                }
            }
        }

        upper.subsystemPointers[tcol + 1] = upper.subsystemPointers[tcol] + nv;
    }

    if(verbosity >= 4){
        std::ostringstream out;
        out << "BISAI buildUpperSubsystemsStructures time: " << t_buildUpperSubsystemsStructures.stop() << " s";
        OpmLog::info(out.str());
    }
}

template <unsigned int block_size>
bool BISAI<block_size>::create_preconditioner(BlockedMatrix *mat, BlockedMatrix *jacMat)
{
    const unsigned int bs = block_size;

    if (bs != 3) {
        OPM_THROW(std::logic_error, "Creation of ISAI preconditioner on GPU only supports block_size = 3");
    }

    Dune::Timer t_preconditioner;

    if (jacMat) {
        bilu0->create_preconditioner(mat, jacMat);
    } else {
        bilu0->create_preconditioner(mat);
    }

    std::call_once(initialize, [&]() {
        std::tie(colPointers, rowIndices, diagIndex) = bilu0->get_preconditioner_structure();

        csrToCscOffsetMap = buildCsrToCscOffsetMap(colPointers, rowIndices);
        buildLowerSubsystemsStructures();
        buildUpperSubsystemsStructures();

        d_colPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * colPointers.size());
        d_rowIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * rowIndices.size());
        d_csrToCscOffsetMap = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * csrToCscOffsetMap.size());
        d_diagIndex = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * diagIndex.size());
        d_invLvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * nnzb * bs * bs);
        d_invUvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * nnzb * bs * bs);
        d_invL_x = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * Nb * bs);
        d_lower.subsystemPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * lower.subsystemPointers.size());
        d_upper.subsystemPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * upper.subsystemPointers.size());

        if(!lower.nzIndices.empty()){ // knownRhsIndices and unknownRhsIndices will also be empty if nzIndices is empty
            d_lower.nzIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * lower.nzIndices.size());
            d_lower.knownRhsIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * lower.knownRhsIndices.size());
            d_lower.unknownRhsIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * lower.unknownRhsIndices.size());
        }

        if(!upper.nzIndices.empty()){ // knownRhsIndices and unknownRhsIndices will also be empty if nzIndices is empty
            d_upper.nzIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * upper.nzIndices.size());
            d_upper.knownRhsIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * upper.knownRhsIndices.size());
            d_upper.unknownRhsIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * upper.unknownRhsIndices.size());
        }

        events.resize(6);
        err = queue->enqueueWriteBuffer(d_colPointers, CL_FALSE, 0, colPointers.size() * sizeof(int), colPointers.data(), nullptr, &events[0]);
        err |= queue->enqueueWriteBuffer(d_rowIndices, CL_FALSE, 0, rowIndices.size() * sizeof(int), rowIndices.data(), nullptr, &events[1]);
        err |= queue->enqueueWriteBuffer(d_csrToCscOffsetMap, CL_FALSE, 0, csrToCscOffsetMap.size() * sizeof(int), csrToCscOffsetMap.data(), nullptr, &events[2]);
        err |= queue->enqueueWriteBuffer(d_diagIndex, CL_FALSE, 0, diagIndex.size() * sizeof(int), diagIndex.data(), nullptr, &events[3]);
        err |= queue->enqueueWriteBuffer(d_lower.subsystemPointers, CL_FALSE, 0, sizeof(int) * lower.subsystemPointers.size(), lower.subsystemPointers.data(), nullptr, &events[4]);
        err |= queue->enqueueWriteBuffer(d_upper.subsystemPointers, CL_FALSE, 0, sizeof(int) * upper.subsystemPointers.size(), upper.subsystemPointers.data(), nullptr, &events[5]);

        if(!lower.nzIndices.empty()){
            events.resize(events.size() + 3);
            err |= queue->enqueueWriteBuffer(d_lower.nzIndices, CL_FALSE, 0, sizeof(int) * lower.nzIndices.size(), lower.nzIndices.data(), nullptr, &events[events.size() - 3]);
            err |= queue->enqueueWriteBuffer(d_lower.knownRhsIndices, CL_FALSE, 0, sizeof(int) * lower.knownRhsIndices.size(), lower.knownRhsIndices.data(), nullptr, &events[events.size() - 2]);
            err |= queue->enqueueWriteBuffer(d_lower.unknownRhsIndices, CL_FALSE, 0, sizeof(int) * lower.unknownRhsIndices.size(), lower.unknownRhsIndices.data(), nullptr, &events[events.size() - 1]);
        }

        if(!upper.nzIndices.empty()){
            events.resize(events.size() + 3);
            err |= queue->enqueueWriteBuffer(d_upper.nzIndices, CL_FALSE, 0, sizeof(int) * upper.nzIndices.size(), upper.nzIndices.data(), nullptr, &events[events.size() - 3]);
            err |= queue->enqueueWriteBuffer(d_upper.knownRhsIndices, CL_FALSE, 0, sizeof(int) * upper.knownRhsIndices.size(), upper.knownRhsIndices.data(), nullptr, &events[events.size() - 2]);
            err |= queue->enqueueWriteBuffer(d_upper.unknownRhsIndices, CL_FALSE, 0, sizeof(int) * upper.unknownRhsIndices.size(), upper.unknownRhsIndices.data(), nullptr, &events[events.size() - 1]);
        }

        cl::WaitForEvents(events);
        events.clear();

        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "BISAI OpenCL enqueueWriteBuffer error");
        }
    });

    std::tie(d_LUvals, d_invDiagVals) = bilu0->get_preconditioner_data();

    events.resize(2);
    err = queue->enqueueFillBuffer(d_invLvals, 0, 0, sizeof(double) * nnzb * bs * bs, nullptr, &events[0]);
    err |= queue->enqueueFillBuffer(d_invUvals, 0, 0, sizeof(double) * nnzb * bs * bs, nullptr, &events[1]);
    cl::WaitForEvents(events);
    events.clear();

    OpenclKernels::isaiL(d_diagIndex, d_colPointers, d_csrToCscOffsetMap, d_lower.subsystemPointers, d_lower.nzIndices, d_lower.unknownRhsIndices, d_lower.knownRhsIndices, d_LUvals, d_invLvals, Nb);
    OpenclKernels::isaiU(d_diagIndex, d_colPointers, d_rowIndices, d_csrToCscOffsetMap, d_upper.subsystemPointers, d_upper.nzIndices, d_upper.unknownRhsIndices, d_upper.knownRhsIndices, d_LUvals,
            d_invDiagVals, d_invUvals, Nb);

    if(verbosity >= 4){
        std::ostringstream out;
        out << "BISAI createPreconditioner time: " << t_preconditioner.stop() << " s";
        OpmLog::info(out.str());
    }

    return true;
}

template <unsigned int block_size>
bool BISAI<block_size>::create_preconditioner(BlockedMatrix *mat)
{
    return create_preconditioner(mat, nullptr);
}

template <unsigned int block_size>
void BISAI<block_size>::apply(const cl::Buffer& x, cl::Buffer& y){
    const unsigned int bs = block_size;

    OpenclKernels::spmv(d_invLvals, d_rowIndices, d_colPointers, x, d_invL_x, Nb, bs, true, true); // application of isaiL is a simple spmv with addition
                                                                                                   // (to compensate for the unitary diagonal that is not
                                                                                                   // included in isaiL, for simplicity)
    OpenclKernels::spmv(d_invUvals, d_rowIndices, d_colPointers, d_invL_x, y, Nb, bs); // application of isaiU is a simple spmv
}

#define INSTANTIATE_BDA_FUNCTIONS(n)  \
template class BISAI<n>;

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

}
}
