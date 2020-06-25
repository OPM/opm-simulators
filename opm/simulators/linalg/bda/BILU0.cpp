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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/BILU0.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


namespace bda
{

    using Opm::OpmLog;

    // define 'second' as 'BdaSolver<>::second', this allows usage of the second() function for timing
    // typedefs cannot handle templates
    const auto second = BdaSolver<>::second;

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
        delete[] rowsPerColor;
        freeBlockedMatrix(&LUMat);
        freeBlockedMatrix(&LMat);
        freeBlockedMatrix(&UMat);
        delete[] toOrder;
        delete[] fromOrder;
        freeBlockedMatrix(&rMat);
    }

    template <unsigned int block_size>
    bool BILU0<block_size>::init(BlockedMatrix *mat)
    {
        double t1 = 0.0, t2 = 0.0;
        BlockedMatrix *CSCmat = nullptr;

        this->N = mat->Nb * block_size;
        this->Nb = mat->Nb;
        this->nnz = mat->nnzbs * block_size * block_size;
        this->nnzbs = mat->nnzbs;

        toOrder = new int[N];
        fromOrder = new int[N];
        if (level_scheduling) {
            CSCmat = new BlockedMatrix();
            CSCmat->Nb = Nb;
            CSCmat->nnzbs = nnzbs;
            CSCmat->nnzValues = new Block[nnzbs];
            CSCmat->colIndices = new int[nnzbs];
            CSCmat->rowPointers = new int[Nb + 1];
            if(verbosity >= 3){
                t1 = second();
            }
            bcsr_to_bcsc(mat->nnzValues, mat->colIndices, mat->rowPointers, CSCmat->nnzValues, CSCmat->colIndices, CSCmat->rowPointers, mat->Nb);
            if(verbosity >= 3){
                t2 = second();
                std::ostringstream out;
                out << "BILU0 convert CSR to CSC: " << t2 - t1 << " s";
                OpmLog::info(out.str());
            }
        }

        if(verbosity >= 3){
            t1 = second();
        }
        rMat = allocateBlockedMatrix(mat->Nb, mat->nnzbs);
        LUMat = soft_copyBlockedMatrix(rMat);
        if (level_scheduling) {
            rowsPerColor = findLevelScheduling(mat->colIndices, mat->rowPointers, CSCmat->colIndices, CSCmat->rowPointers, mat->Nb, &numColors, toOrder, fromOrder);
        } else if (graph_coloring) {
            rowsPerColor = findGraphColoring(mat->colIndices, mat->rowPointers, mat->Nb, mat->Nb, mat->Nb, &numColors, toOrder, fromOrder);
        }
        if (rowsPerColor == nullptr) {
            return false;
        }
        if(verbosity >= 3){
            t2 = second();
            std::ostringstream out;
            out << "BILU0 analysis took: " << t2 - t1 << " s, " << numColors << " colors";
            OpmLog::info(out.str());
        }

        diagIndex = new int[mat->Nb];
        invDiagVals = new Block[mat->Nb];

        LMat = allocateBlockedMatrix(mat->Nb, (mat->nnzbs - mat->Nb) / 2);
        UMat = allocateBlockedMatrix(mat->Nb, (mat->nnzbs - mat->Nb) / 2);

        LUMat->nnzValues = new Block[mat->nnzbs];

        if (level_scheduling) {
            delete[] CSCmat->nnzValues;
            delete[] CSCmat->colIndices;
            delete[] CSCmat->rowPointers;
            delete CSCmat;
        }


        s.Lvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Block) * LMat->nnzbs);
        s.Uvals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Block) * UMat->nnzbs);
        s.Lcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * LMat->nnzbs);
        s.Ucols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * UMat->nnzbs);
        s.Lrows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (LMat->Nb + 1));
        s.Urows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (UMat->Nb + 1));
        s.invDiagVals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(Block) * mat->Nb);
        s.rowsPerColor = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (numColors + 1));

        queue->enqueueWriteBuffer(s.Lvals, CL_TRUE, 0, LMat->nnzbs * sizeof(Block), LMat->nnzValues);
        queue->enqueueWriteBuffer(s.Uvals, CL_TRUE, 0, UMat->nnzbs * sizeof(Block), UMat->nnzValues);
        queue->enqueueWriteBuffer(s.Lcols, CL_TRUE, 0, LMat->nnzbs * sizeof(int), LMat->colIndices);
        queue->enqueueWriteBuffer(s.Ucols, CL_TRUE, 0, UMat->nnzbs * sizeof(int), UMat->colIndices);
        queue->enqueueWriteBuffer(s.Lrows, CL_TRUE, 0, (LMat->Nb + 1) * sizeof(int), LMat->rowPointers);
        queue->enqueueWriteBuffer(s.Urows, CL_TRUE, 0, (UMat->Nb + 1) * sizeof(int), UMat->rowPointers);
        queue->enqueueWriteBuffer(s.invDiagVals, CL_TRUE, 0, mat->Nb * sizeof(Block), invDiagVals);

        int *rowsPerColorPrefix = new int[numColors + 1];
        int prefix = 0;
        for (int i = 0; i < numColors + 1; ++i) {
            rowsPerColorPrefix[i] = prefix;
            prefix += rowsPerColor[i];
        }
        queue->enqueueWriteBuffer(s.rowsPerColor, CL_TRUE, 0, (numColors + 1) * sizeof(int), rowsPerColorPrefix);
        delete[] rowsPerColorPrefix;

        return true;
    } // end init()


    template <unsigned int block_size>
    bool BILU0<block_size>::create_preconditioner(BlockedMatrix *mat)
    {
        double t1 = 0.0, t2 = 0.0;
        if (verbosity >= 3){
            t1 = second();
        }
        blocked_reorder_matrix_by_pattern(mat, toOrder, fromOrder, rMat);
        if (verbosity >= 3){
            t2 = second();
            std::ostringstream out;
            out << "BILU0 reorder matrix: " << t2 - t1 << " s";
            OpmLog::info(out.str());
        }

        // TODO: remove this copy by replacing inplace ilu decomp by out-of-place ilu decomp
        if (verbosity >= 3){
            t1 = second();
        }
        memcpy(LUMat->nnzValues, rMat->nnzValues, sizeof(Block) * rMat->nnzbs);
        if (verbosity >= 3){
            t2 = second();
            std::ostringstream out;
            out << "BILU0 memcpy: " << t2 - t1 << " s";
            OpmLog::info(out.str());
        }

        int i, j, ij, ik, jk;
        int iRowStart, iRowEnd, jRowEnd;
        Block pivot;

        int LSize = 0;
        const int blockSquare = block_size * block_size;

        if (verbosity >= 3){
            t1 = second();
        }
        // go through all rows
        for (i = 0; i < LUMat->Nb; i++) {
            iRowStart = LUMat->rowPointers[i];
            iRowEnd = LUMat->rowPointers[i + 1];

            // go through all elements of the row
            for (ij = iRowStart; ij < iRowEnd; ij++) {
                j = LUMat->colIndices[ij];
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
                blockMult(LUMat->nnzValues[ij], invDiagVals[j], pivot);

                memcpy(LUMat->nnzValues[ij], &pivot, sizeof(double)*block_size * block_size);

                jRowEnd = LUMat->rowPointers[j + 1];
                jk = diagIndex[j] + 1;
                ik = ij + 1;
                // substract that row scaled by the pivot from this row.
                while (ik < iRowEnd && jk < jRowEnd) {
                    if (LUMat->colIndices[ik] == LUMat->colIndices[jk]) {
                        blockMultSub(LUMat->nnzValues[ik], pivot, LUMat->nnzValues[jk]);
                        ik++;
                        jk++;
                    } else {
                        if (LUMat->colIndices[ik] < LUMat->colIndices[jk])
                        { ik++; }
                        else
                        { jk++; }
                    }
                }
            }
            // store the inverse in the diagonal!
            blockInvert3x3(LUMat->nnzValues[ij], invDiagVals[i]);

            memcpy(LUMat->nnzValues[ij], invDiagVals[i], sizeof(double)*block_size * block_size);
        }

        LMat->rowPointers[0] = 0;
        UMat->rowPointers[0] = 0;

        // Split the LU matrix into two by comparing column indices to diagonal indices
        for (i = 0; i < LUMat->Nb; i++) {
            int offsetL = LMat->rowPointers[i];
            int rowSize = diagIndex[i] - LUMat->rowPointers[i];
            int offsetLU = LUMat->rowPointers[i];
            memcpy(LMat->nnzValues[offsetL], LUMat->nnzValues[offsetLU], sizeof(double) * blockSquare * rowSize);
            memcpy(LMat->colIndices + offsetL, LUMat->colIndices + offsetLU, sizeof(int) * rowSize);
            offsetL += rowSize;
            LMat->rowPointers[i + 1] = offsetL;
        }
        // Reverse the order or the (blocked) rows for the U matrix,
        // because the rows are accessed in reverse order when applying the ILU0
        int URowIndex = 0;
        for (i = LUMat->Nb - 1; i >= 0; i--) {
            int offsetU = UMat->rowPointers[URowIndex];
            int rowSize = LUMat->rowPointers[i + 1] - diagIndex[i] - 1;
            int offsetLU = diagIndex[i] + 1;
            memcpy(UMat->nnzValues[offsetU], LUMat->nnzValues[offsetLU], sizeof(double) * blockSquare * rowSize);
            memcpy(UMat->colIndices + offsetU, LUMat->colIndices + offsetLU, sizeof(int) * rowSize);
            offsetU += rowSize;
            UMat->rowPointers[URowIndex + 1] = offsetU;
            URowIndex++;
        }
        if (verbosity >= 3) {
            t2 = second();
            std::ostringstream out;
            out << "BILU0 decomposition: " << t2 - t1 << " s";
            OpmLog::info(out.str());
        }

        if (verbosity >= 3) {
            t1 = second();
        }
        if (pattern_uploaded == false) {
            queue->enqueueWriteBuffer(s.Lcols, CL_TRUE, 0, LMat->nnzbs * sizeof(int), LMat->colIndices);
            queue->enqueueWriteBuffer(s.Ucols, CL_TRUE, 0, UMat->nnzbs * sizeof(int), UMat->colIndices);
            queue->enqueueWriteBuffer(s.Lrows, CL_TRUE, 0, (LMat->Nb + 1) * sizeof(int), LMat->rowPointers);
            queue->enqueueWriteBuffer(s.Urows, CL_TRUE, 0, (UMat->Nb + 1) * sizeof(int), UMat->rowPointers);
            pattern_uploaded = true;
        }
        queue->enqueueWriteBuffer(s.Lvals, CL_TRUE, 0, LMat->nnzbs * sizeof(Block), LMat->nnzValues);
        queue->enqueueWriteBuffer(s.Uvals, CL_TRUE, 0, UMat->nnzbs * sizeof(Block), UMat->nnzValues);
        queue->enqueueWriteBuffer(s.invDiagVals, CL_TRUE, 0, Nb * sizeof(Block), invDiagVals);
        if (verbosity >= 3) {
            t2 = second();
            std::ostringstream out;
            out << "BILU0 copy to GPU: " << t2 - t1 << " s";
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
        double t1 = 0.0, t2 = 0.0;
        if (verbosity >= 3) {
            t1 = second();
        }
        cl::Event event;

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
            t2 = second();
            std::ostringstream out;
            out << "BILU0 apply: " << t2 - t1 << " s";
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
template bool BILU0<n>::init(BlockedMatrix*);                                            \
template bool BILU0<n>::create_preconditioner(BlockedMatrix*);                           \
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


