/*
  Copyright 2019 Equinor ASA

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

#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>
#include <opm/simulators/linalg/bda/opencl/ChowPatelIlu.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>

#include <sstream>

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
BILU0<block_size>::BILU0(bool opencl_ilu_parallel_, int verbosity_) :
    Preconditioner<block_size>(verbosity_), opencl_ilu_parallel(opencl_ilu_parallel_)
{
#if CHOW_PATEL
    chowPatelIlu.setVerbosity(verbosity);
#endif
}


template <unsigned int block_size>
bool BILU0<block_size>::analyze_matrix(BlockedMatrix *mat)
{
    return analyze_matrix(mat, nullptr);
}


template <unsigned int block_size>
bool BILU0<block_size>::analyze_matrix(BlockedMatrix *mat, BlockedMatrix *jacMat)
{
    const unsigned int bs = block_size;

    this->N = mat->Nb * block_size;
    this->Nb = mat->Nb;
    this->nnz = mat->nnzbs * block_size * block_size;
    this->nnzb = mat->nnzbs;

    std::vector<int> CSCRowIndices;
    std::vector<int> CSCColPointers;

    auto *matToDecompose = jacMat ? jacMat : mat; // decompose jacMat if valid, otherwise decompose mat

    if (opencl_ilu_parallel) {
        toOrder.resize(Nb);
        fromOrder.resize(Nb);
        CSCRowIndices.resize(matToDecompose->nnzbs);
        CSCColPointers.resize(Nb + 1);

        LUmat = std::make_unique<BlockedMatrix>(*matToDecompose);

        Timer t_convert;
        csrPatternToCsc(matToDecompose->colIndices, matToDecompose->rowPointers, CSCRowIndices.data(), CSCColPointers.data(), Nb);
        if(verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 convert CSR to CSC: " << t_convert.stop() << " s";
            OpmLog::info(out.str());
        }
    } else {
        LUmat = std::make_unique<BlockedMatrix>(*matToDecompose);
    }

    Timer t_analysis;
    std::ostringstream out;
    if (opencl_ilu_parallel) {
        out << "opencl_ilu_parallel: true (level_scheduling)\n";
        findLevelScheduling(matToDecompose->colIndices, matToDecompose->rowPointers, CSCRowIndices.data(), CSCColPointers.data(), Nb, &numColors, toOrder.data(), fromOrder.data(), rowsPerColor);
    } else {
        out << "opencl_ilu_parallel: false\n";
        // numColors = 1;
        // rowsPerColor.emplace_back(Nb);
        numColors = Nb;
        for(int i = 0; i < Nb; ++i){
            rowsPerColor.emplace_back(1);
        }
    }

    if (verbosity >= 1) {
        out << "BILU0 analysis took: " << t_analysis.stop() << " s, " << numColors << " colors\n";
    }

#if CHOW_PATEL
    out << "BILU0 CHOW_PATEL: " << CHOW_PATEL << ", CHOW_PATEL_GPU: " << CHOW_PATEL_GPU;
#endif
    OpmLog::info(out.str());

    diagIndex.resize(mat->Nb);
    invDiagVals.resize(mat->Nb * bs * bs);

#if CHOW_PATEL
    Lmat = std::make_unique<BlockedMatrix>(mat->Nb, (mat->nnzbs - mat->Nb) / 2, block_size);
    Umat = std::make_unique<BlockedMatrix>(mat->Nb, (mat->nnzbs - mat->Nb) / 2, block_size);
#endif

    s.invDiagVals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * mat->Nb);
    s.rowsPerColor = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (numColors + 1));
    s.diagIndex = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * LUmat->Nb);
    s.rowIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(unsigned) * LUmat->Nb);
#if CHOW_PATEL
    s.Lvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * Lmat->nnzbs);
    s.Lcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Lmat->nnzbs);
    s.Lrows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Lmat->Nb + 1));
    s.Uvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * Lmat->nnzbs);
    s.Ucols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Lmat->nnzbs);
    s.Urows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Lmat->Nb + 1));
#else
    s.LUvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * LUmat->nnzbs);
    s.LUcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * LUmat->nnzbs);
    s.LUrows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (LUmat->Nb + 1));
#endif

    events.resize(3);
    err = queue->enqueueWriteBuffer(s.invDiagVals, CL_FALSE, 0, mat->Nb * sizeof(double) * bs * bs, invDiagVals.data(), nullptr, &events[0]);

    rowsPerColorPrefix.resize(numColors + 1); // resize initializes value 0.0
    for (int i = 0; i < numColors; ++i) {
        rowsPerColorPrefix[i + 1] = rowsPerColorPrefix[i] + rowsPerColor[i];
    }

    err |= queue->enqueueWriteBuffer(s.rowsPerColor, CL_FALSE, 0, (numColors + 1) * sizeof(int), rowsPerColorPrefix.data(), nullptr, &events[1]);

    if (opencl_ilu_parallel) {
        err |= queue->enqueueWriteBuffer(s.rowIndices, CL_FALSE, 0, Nb * sizeof(unsigned), fromOrder.data(), nullptr, &events[2]);
    } else {
        // fromOrder is not initialized, so use something else to fill s.rowIndices
        // s.rowIndices[i] == i must hold, since every rowidx is mapped to itself (i.e. no actual mapping)
        // rowsPerColorPrefix is misused here, it contains an increasing sequence (0, 1, 2, ...)
        err |= queue->enqueueWriteBuffer(s.rowIndices, CL_FALSE, 0, Nb * sizeof(unsigned), rowsPerColorPrefix.data(), nullptr, &events[2]);
    }

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "BILU0 OpenCL enqueueWriteBuffer error");
    }

    return true;
}



template <unsigned int block_size>
bool BILU0<block_size>::create_preconditioner(BlockedMatrix *mat)
{
    return create_preconditioner(mat, nullptr);
}


template <unsigned int block_size>
bool BILU0<block_size>::create_preconditioner(BlockedMatrix *mat, BlockedMatrix *jacMat)
{
    const unsigned int bs = block_size;

    auto *matToDecompose = jacMat ? jacMat : mat;

    // TODO: remove this copy by replacing inplace ilu decomp by out-of-place ilu decomp
    Timer t_copy;
    memcpy(LUmat->nnzValues, matToDecompose->nnzValues, sizeof(double) * bs * bs * matToDecompose->nnzbs);

    if (verbosity >= 3){
        std::ostringstream out;
        out << "BILU0 memcpy: " << t_copy.stop() << " s";
        OpmLog::info(out.str());
    }

#if CHOW_PATEL
    chowPatelIlu.decomposition(queue.get(), context.get(),
                               LUmat.get(), Lmat.get(), Umat.get(),
                               invDiagVals.data(), diagIndex,
                               s.diagIndex, s.invDiagVals,
                               s.Lvals, s.Lcols, s.Lrows,
                               s.Uvals, s.Ucols, s.Urows);
#else
    Timer t_copyToGpu;

    events.resize(1);
    queue->enqueueWriteBuffer(s.LUvals, CL_FALSE, 0, LUmat->nnzbs * bs * bs * sizeof(double), LUmat->nnzValues, nullptr, &events[0]);

    std::call_once(pattern_uploaded, [&](){
        // find the positions of each diagonal block
        for (int row = 0; row < Nb; ++row) {
            int rowStart = LUmat->rowPointers[row];
            int rowEnd = LUmat->rowPointers[row+1];

            auto candidate = std::find(LUmat->colIndices + rowStart, LUmat->colIndices + rowEnd, row);
            assert(candidate != LUmat->colIndices + rowEnd);
            diagIndex[row] = candidate - LUmat->colIndices;
        }
        events.resize(4);
        queue->enqueueWriteBuffer(s.diagIndex, CL_FALSE, 0, Nb * sizeof(int), diagIndex.data(), nullptr, &events[1]);
        queue->enqueueWriteBuffer(s.LUcols, CL_FALSE, 0, LUmat->nnzbs * sizeof(int), LUmat->colIndices, nullptr, &events[2]);
        queue->enqueueWriteBuffer(s.LUrows, CL_FALSE, 0, (LUmat->Nb + 1) * sizeof(int), LUmat->rowPointers, nullptr, &events[3]);
    });

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "BILU0 OpenCL enqueueWriteBuffer error");
    }

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "BILU0 copy to GPU: " << t_copyToGpu.stop() << " s";
        OpmLog::info(out.str());
    }

    Timer t_decomposition;
    std::ostringstream out;
    for (int color = 0; color < numColors; ++color) {
        const unsigned int firstRow = rowsPerColorPrefix[color];
        const unsigned int lastRow = rowsPerColorPrefix[color + 1];
        if (verbosity >= 5) {
            out << "color " << color << ": " << firstRow << " - " << lastRow << " = " << lastRow - firstRow << "\n";
        }
        OpenclKernels::ILU_decomp(firstRow, lastRow, s.rowIndices,
                                  s.LUvals, s.LUcols, s.LUrows, s.diagIndex,
                                  s.invDiagVals, rowsPerColor[color], block_size);
    }

    if (verbosity >= 3) {
        queue->finish();
        out << "BILU0 decomposition: " << t_decomposition.stop() << " s";
        OpmLog::info(out.str());
    }
#endif // CHOW_PATEL

    return true;
} // end create_preconditioner()


// kernels are blocking on an NVIDIA GPU, so waiting for events is not needed
// however, if individual kernel calls are timed, waiting for events is needed
// behavior on other GPUs is untested
template <unsigned int block_size>
void BILU0<block_size>::apply(const cl::Buffer& y, cl::Buffer& x)
{
    const double relaxation = 0.9;
    cl::Event event;
    Timer t_apply;

    for (int color = 0; color < numColors; ++color) {
#if CHOW_PATEL
        OpenclKernels::ILU_apply1(s.rowIndices, s.Lvals, s.Lcols, s.Lrows,
                                  s.diagIndex, y, x, s.rowsPerColor,
                                  color, rowsPerColor[color], block_size);
#else
        OpenclKernels::ILU_apply1(s.rowIndices, s.LUvals, s.LUcols, s.LUrows,
                                  s.diagIndex, y, x, s.rowsPerColor,
                                  color, rowsPerColor[color], block_size);
#endif
    }

    for (int color = numColors - 1; color >= 0; --color) {
#if CHOW_PATEL
        OpenclKernels::ILU_apply2(s.rowIndices, s.Uvals, s.Ucols, s.Urows,
                                  s.diagIndex, s.invDiagVals, x, s.rowsPerColor,
                                  color, rowsPerColor[color], block_size);
#else
        OpenclKernels::ILU_apply2(s.rowIndices, s.LUvals, s.LUcols, s.LUrows,
                                  s.diagIndex, s.invDiagVals, x, s.rowsPerColor,
                                  color, rowsPerColor[color], block_size);
#endif
    }

    // apply relaxation
    OpenclKernels::scale(x, relaxation, N);

    if (verbosity >= 4) {
        std::ostringstream out;
        out << "BILU0 apply: " << t_apply.stop() << " s";
        OpmLog::info(out.str());
    }
}



#define INSTANTIATE_BDA_FUNCTIONS(n) \
template class BILU0<n>;


INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
