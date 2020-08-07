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
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/MatrixBlock.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/BILU0.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


namespace bda
{

    using Opm::OpmLog;
    using Dune::Timer;

    template <unsigned int block_size>
    BILU0<block_size>::BILU0(bool level_scheduling_, bool graph_coloring_, int verbosity_) :
        verbosity(verbosity_), level_scheduling(level_scheduling_), graph_coloring(graph_coloring_)
    {
        if (level_scheduling == graph_coloring) {
            OPM_THROW(std::logic_error, "Error, either level_scheduling or graph_coloring must be true, not both\n");
        }
    }

    template <unsigned int block_size>
    BILU0<block_size>::~BILU0()
    {
        delete[] invDiagVals;
        delete[] diagIndex;
        delete[] toOrder;
        delete[] fromOrder;
    }

    template <unsigned int block_size>
    bool BILU0<block_size>::init(BlockedMatrix<block_size> *mat)
    {
        const unsigned int bs = block_size;

        this->N = mat->Nb * block_size;
        this->Nb = mat->Nb;
        this->nnz = mat->nnzbs * block_size * block_size;
        this->nnzbs = mat->nnzbs;

        toOrder = new int[Nb];
        fromOrder = new int[Nb];

        int *CSCRowIndices = new int[nnzbs];
        int *CSCColPointers = new int[Nb + 1];

        Timer t_convert;
        csrPatternToCsc(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb);
        if(verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 convert CSR to CSC: " << t_convert.stop() << " s";
            OpmLog::info(out.str());
        }

        Timer t_analysis;
        rmat = std::make_shared<BlockedMatrix<block_size> >(mat->Nb, mat->nnzbs);
        LUmat = std::make_unique<BlockedMatrix<block_size> >(*rmat);
        if (level_scheduling) {
            findLevelScheduling(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        } else if (graph_coloring) {
            findGraphColoring<block_size>(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, mat->Nb, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        }
        if(verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 analysis took: " << t_analysis.stop() << " s, " << numColors << " colors";
            OpmLog::info(out.str());
        }

        delete[] CSCRowIndices;
        delete[] CSCColPointers;

        diagIndex = new int[mat->Nb];
        invDiagVals = new double[mat->Nb * bs * bs];

        Lmat = std::make_unique<BlockedMatrix<block_size> >(mat->Nb, (mat->nnzbs - mat->Nb) / 2);
        Umat = std::make_unique<BlockedMatrix<block_size> >(mat->Nb, (mat->nnzbs - mat->Nb) / 2);

        LUmat->nnzValues = new double[mat->nnzbs * bs * bs];

        s.Lvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * Lmat->nnzbs);
        s.Uvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * Umat->nnzbs);
        s.Lcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Lmat->nnzbs);
        s.Ucols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Umat->nnzbs);
        s.Lrows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Lmat->Nb + 1));
        s.Urows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Umat->Nb + 1));
        s.invDiagVals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * mat->Nb);
        s.rowsPerColor = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (numColors + 1));

        queue->enqueueWriteBuffer(s.Lvals, CL_TRUE, 0, Lmat->nnzbs * sizeof(double) * bs * bs, Lmat->nnzValues);
        queue->enqueueWriteBuffer(s.Uvals, CL_TRUE, 0, Umat->nnzbs * sizeof(double) * bs * bs, Umat->nnzValues);
        queue->enqueueWriteBuffer(s.Lcols, CL_TRUE, 0, Lmat->nnzbs * sizeof(int), Lmat->colIndices);
        queue->enqueueWriteBuffer(s.Ucols, CL_TRUE, 0, Umat->nnzbs * sizeof(int), Umat->colIndices);
        queue->enqueueWriteBuffer(s.Lrows, CL_TRUE, 0, (Lmat->Nb + 1) * sizeof(int), Lmat->rowPointers);
        queue->enqueueWriteBuffer(s.Urows, CL_TRUE, 0, (Umat->Nb + 1) * sizeof(int), Umat->rowPointers);
        queue->enqueueWriteBuffer(s.invDiagVals, CL_TRUE, 0, mat->Nb * sizeof(double) * bs * bs, invDiagVals);

        int *rowsPerColorPrefix = new int[numColors + 1];
        rowsPerColorPrefix[0] = 0;
        for (int i = 0; i < numColors; ++i) {
            rowsPerColorPrefix[i+1] = rowsPerColorPrefix[i] + rowsPerColor[i];
        }
        queue->enqueueWriteBuffer(s.rowsPerColor, CL_TRUE, 0, (numColors + 1) * sizeof(int), rowsPerColorPrefix);
        delete[] rowsPerColorPrefix;

        return true;
    } // end init()


    template <unsigned int block_size>
    bool BILU0<block_size>::create_preconditioner(BlockedMatrix<block_size> *mat)
    {
        const unsigned int bs = block_size;
        Timer t_reorder;
        reorderBlockedMatrixByPattern<block_size>(mat, toOrder, fromOrder, rmat.get());

        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 reorder matrix: " << t_reorder.stop() << " s";
            OpmLog::info(out.str());
        }

        // TODO: remove this copy by replacing inplace ilu decomp by out-of-place ilu decomp
        Timer t_copy;
        memcpy(LUmat->nnzValues, rmat->nnzValues, sizeof(double) * bs * bs * rmat->nnzbs);

        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 memcpy: " << t_copy.stop() << " s";
            OpmLog::info(out.str());
        }

        int i, j, ij, ik, jk;
        int iRowStart, iRowEnd, jRowEnd;
        double pivot[bs * bs];

        int LSize = 0;
        Opm::Detail::Inverter<bs> inverter;   // reuse inverter to invert blocks

        Timer t_decomposition;

        // go through all rows
        for (i = 0; i < LUmat->Nb; i++) {
            iRowStart = LUmat->rowPointers[i];
            iRowEnd = LUmat->rowPointers[i + 1];

            // go through all elements of the row
            for (ij = iRowStart; ij < iRowEnd; ij++) {
                j = LUmat->colIndices[ij];
                // if the element is the diagonal, store the index and go to next row
                if (j == i) {
                    diagIndex[i] = ij;
                    break;
                }
                // if an element beyond the diagonal is reach, no diagonal was found
                // throw an error now. TODO: perform reordering earlier to prevent this
                if (j > i) {
                    std::ostringstream out;
                    out << "BILU0 Error could not find diagonal value in row: " << i;
                    OpmLog::error(out.str());
                    return false;
                }

                LSize++;
                // calculate the pivot of this row
                blockMult<bs>(LUmat->nnzValues + ij * bs * bs, invDiagVals + j * bs * bs, &pivot[0]);

                memcpy(LUmat->nnzValues + ij * bs * bs, &pivot[0], sizeof(double) * block_size * block_size);

                jRowEnd = LUmat->rowPointers[j + 1];
                jk = diagIndex[j] + 1;
                ik = ij + 1;
                // substract that row scaled by the pivot from this row.
                while (ik < iRowEnd && jk < jRowEnd) {
                    if (LUmat->colIndices[ik] == LUmat->colIndices[jk]) {
                        blockMultSub<bs>(LUmat->nnzValues + ik * bs * bs, pivot, LUmat->nnzValues + jk * bs * bs);
                        ik++;
                        jk++;
                    } else {
                        if (LUmat->colIndices[ik] < LUmat->colIndices[jk])
                        { ik++; }
                        else
                        { jk++; }
                    }
                }
            }
            // store the inverse in the diagonal!
            inverter(LUmat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs);

            memcpy(LUmat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs, sizeof(double) * bs * bs);
        }

        Lmat->rowPointers[0] = 0;
        Umat->rowPointers[0] = 0;

        // Split the LU matrix into two by comparing column indices to diagonal indices
        for (i = 0; i < LUmat->Nb; i++) {
            int offsetL = Lmat->rowPointers[i];
            int rowSize = diagIndex[i] - LUmat->rowPointers[i];
            int offsetLU = LUmat->rowPointers[i];
            memcpy(Lmat->nnzValues + offsetL * bs * bs, LUmat->nnzValues + offsetLU * bs * bs, sizeof(double) * bs * bs * rowSize);
            memcpy(Lmat->colIndices + offsetL, LUmat->colIndices + offsetLU, sizeof(int) * rowSize);
            offsetL += rowSize;
            Lmat->rowPointers[i + 1] = offsetL;
        }
        // Reverse the order or the (blocked) rows for the U matrix,
        // because the rows are accessed in reverse order when applying the ILU0
        int URowIndex = 0;
        for (i = LUmat->Nb - 1; i >= 0; i--) {
            int offsetU = Umat->rowPointers[URowIndex];
            int rowSize = LUmat->rowPointers[i + 1] - diagIndex[i] - 1;
            int offsetLU = diagIndex[i] + 1;
            memcpy(Umat->nnzValues + offsetU * bs * bs, LUmat->nnzValues + offsetLU * bs * bs, sizeof(double) * bs * bs * rowSize);
            memcpy(Umat->colIndices + offsetU, LUmat->colIndices + offsetLU, sizeof(int) * rowSize);
            offsetU += rowSize;
            Umat->rowPointers[URowIndex + 1] = offsetU;
            URowIndex++;
        }
        if (verbosity >= 3) {
            std::ostringstream out;
            out << "BILU0 decomposition: " << t_decomposition.stop() << " s";
            OpmLog::info(out.str());
        }

        Timer t_copyToGpu;
        if (pattern_uploaded == false) {
            queue->enqueueWriteBuffer(s.Lcols, CL_TRUE, 0, Lmat->nnzbs * sizeof(int), Lmat->colIndices);
            queue->enqueueWriteBuffer(s.Ucols, CL_TRUE, 0, Umat->nnzbs * sizeof(int), Umat->colIndices);
            queue->enqueueWriteBuffer(s.Lrows, CL_TRUE, 0, (Lmat->Nb + 1) * sizeof(int), Lmat->rowPointers);
            queue->enqueueWriteBuffer(s.Urows, CL_TRUE, 0, (Umat->Nb + 1) * sizeof(int), Umat->rowPointers);
            pattern_uploaded = true;
        }
        queue->enqueueWriteBuffer(s.Lvals, CL_TRUE, 0, Lmat->nnzbs * sizeof(double) * bs * bs, Lmat->nnzValues);
        queue->enqueueWriteBuffer(s.Uvals, CL_TRUE, 0, Umat->nnzbs * sizeof(double) * bs * bs, Umat->nnzValues);
        queue->enqueueWriteBuffer(s.invDiagVals, CL_TRUE, 0, Nb * sizeof(double) * bs * bs, invDiagVals);
        if (verbosity >= 3) {
            std::ostringstream out;
            out << "BILU0 copy to GPU: " << t_copyToGpu.stop() << " s";
            OpmLog::info(out.str());
        }

        return true;
    } // end create_preconditioner()

    // kernels are blocking on an NVIDIA GPU, so waiting for events is not needed
    // however, if individual kernel calls are timed, waiting for events is needed
    // behavior on other GPUs is untested
    template <unsigned int block_size>
    void BILU0<block_size>::apply(cl::Buffer& x, cl::Buffer& y)
    {
        cl::Event event;
        Timer t_apply;

        for(int color = 0; color < numColors; ++color){
            event = (*ILU_apply1)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.Lvals, s.Lcols, s.Lrows, (unsigned int)Nb, x, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
            // event.wait();
        }
        for(int color = numColors-1; color >= 0; --color){
            event = (*ILU_apply2)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.Uvals, s.Ucols, s.Urows, (unsigned int)Nb, s.invDiagVals, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
            // event.wait();
        }

        if (verbosity >= 3) {
            event.wait();
            std::ostringstream out;
            out << "BILU0 apply: " << t_apply.stop() << " s";
            OpmLog::info(out.str());
        }
    }


    template <unsigned int block_size>
    void BILU0<block_size>::setOpenCLContext(cl::Context *context_){
        this->context = context_;
    }
    template <unsigned int block_size>
    void BILU0<block_size>::setOpenCLQueue(cl::CommandQueue *queue_){
        this->queue = queue_;
    }
    template <unsigned int block_size>
    void BILU0<block_size>::setKernelParameters(const unsigned int work_group_size_, const unsigned int total_work_items_, const unsigned int lmem_per_work_group_){
        this->work_group_size = work_group_size_;
        this->total_work_items = total_work_items_;
        this->lmem_per_work_group = lmem_per_work_group_;
    }
    template <unsigned int block_size>
    void BILU0<block_size>::setKernels(
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply1_,
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply2_
    ){
        this->ILU_apply1 = ILU_apply1_;
        this->ILU_apply2 = ILU_apply2_;
    }


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                     \
template BILU0<n>::BILU0(bool, bool, int);                                               \
template BILU0<n>::~BILU0();                                                             \
template bool BILU0<n>::init(BlockedMatrix<n>*);                                         \
template bool BILU0<n>::create_preconditioner(BlockedMatrix<n>*);                        \
template void BILU0<n>::apply(cl::Buffer& x, cl::Buffer& y);                             \
template void BILU0<n>::setOpenCLContext(cl::Context*);                                  \
template void BILU0<n>::setOpenCLQueue(cl::CommandQueue*);                               \
template void BILU0<n>::setKernelParameters(unsigned int, unsigned int, unsigned int);   \
template void BILU0<n>::setKernels(                                                      \
    cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *, \
    cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *  \
    );

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // end namespace bda


