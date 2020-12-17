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
#include <opm/simulators/linalg/bda/fgpilu.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


namespace bda
{

// if CHOW_PATEL is 0, exact ILU decomposition is performed on CPU
// if CHOW_PATEL is 1, iterative ILU decomposition (FGPILU) is done, as described in:
//    FINE-GRAINED PARALLEL INCOMPLETE LU FACTORIZATION, E. Chow and A. Patel, SIAM 2015, https://doi.org/10.1137/140968896
// if CHOW_PATEL_GPU is 0, the decomposition is done on CPU
// if CHOW_PATEL_GPU is 1, the decomposition is done by bda::FGPILU::decomposition() on GPU
#define CHOW_PATEL     1
#define CHOW_PATEL_GPU 1

    using Opm::OpmLog;
    using Dune::Timer;

    template <unsigned int block_size>
    BILU0<block_size>::BILU0(ILUReorder opencl_ilu_reorder_, int verbosity_) :
        verbosity(verbosity_), opencl_ilu_reorder(opencl_ilu_reorder_)
    {}

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
        std::ostringstream out;
        if (opencl_ilu_reorder == ILUReorder::LEVEL_SCHEDULING) {
            out << "BILU0 reordering strategy: " << "level_scheduling\n";
            findLevelScheduling(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        } else if (opencl_ilu_reorder == ILUReorder::GRAPH_COLORING) {
            out << "BILU0 reordering strategy: " << "graph_coloring\n";
            findGraphColoring<block_size>(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, mat->Nb, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        } else if (opencl_ilu_reorder == ILUReorder::NONE) {
            out << "BILU0 reordering strategy: none\n";
            numColors = 1;
            rowsPerColor.emplace_back(Nb);
            // numColors = Nb;
            // for(int i = 0; i < Nb; ++i){
            //     rowsPerColor.emplace_back(1);
            // }
            for(int i = 0; i < Nb; ++i){
                toOrder[i] = i;
                fromOrder[i] = i;
            }
        } else {
            OPM_THROW(std::logic_error, "Error ilu reordering strategy not set correctly\n");
        }
        if(verbosity >= 3){
            out << "BILU0 analysis took: " << t_analysis.stop() << " s, " << numColors << " colors\n";
        }
        if(verbosity >= 2){
            out << "BILU0 CHOW_PATEL: " << CHOW_PATEL << ", CHOW_PATEL_GPU: " << CHOW_PATEL_GPU;
        }
        OpmLog::info(out.str());

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

    // implements Fine-Grained Parallel ILU algorithm (FGPILU), Chow and Patel 2015
    template <unsigned int block_size>
    void BILU0<block_size>::chow_patel_decomposition()
    {
        const unsigned int bs = block_size;
        int num_sweeps = 6;

        // split matrix into L and U
        // also convert U into BSC format (Ut)
        // Ut stores diagonal for now
        int num_blocks_L = 0;

        // Ut is actually BSC format
        std::unique_ptr<BlockedMatrix<bs> > Ut = std::make_unique<BlockedMatrix<bs> >(Umat->Nb, Umat->nnzbs + Umat->Nb);

        Lmat->rowPointers[0] = 0;
        for (int i = 0; i < Nb+1; i++) {
            Ut->rowPointers[i] = 0;
        }

        // for every row
        for (int i = 0; i < Nb; i++) {
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block in this row
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                if (i <= j) {
                    Ut->rowPointers[j+1]++;   // actually colPointers
                } else {
                    Lmat->colIndices[num_blocks_L] = j;
                    memcpy(Lmat->nnzValues + num_blocks_L * bs * bs, LUmat->nnzValues + ij * bs * bs, sizeof(double) * bs * bs);
                    num_blocks_L++;
                }
            }
            Lmat->rowPointers[i+1] = num_blocks_L;
        }

        // prefix sum
        int sum = 0;
        for (int i = 1; i < Nb+1; i++) {
            sum += Ut->rowPointers[i];
            Ut->rowPointers[i] = sum;
        }

        // for every row
        for (int i = 0; i < Nb; i++) {
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block in this row
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                if (i <= j){
                    int idx = Ut->rowPointers[j]++;
                    Ut->colIndices[idx] = i;     // actually rowIndices
                    memcpy(Ut->nnzValues + idx * bs * bs, LUmat->nnzValues + ij * bs * bs, sizeof(double) * bs * bs);
                }
            }
        }

        // rotate
        // the Ut->rowPointers were increased in the last loop
        for (int i = Nb; i > 0; --i) {
            Ut->rowPointers[i] = Ut->rowPointers[i-1];
        }
        Ut->rowPointers[0] = 0;

        Opm::Detail::Inverter<bs> inverter;

        // Utmp is needed for CPU and GPU decomposition, because U is transposed, and reversed at the end
        // Ltmp is only needed for CPU decomposition, GPU creates GPU buffer for Ltmp
        double *Utmp = new double[Ut->nnzbs * block_size * block_size];

        // actual ILU decomposition
#if CHOW_PATEL_GPU
        fgpilu.decomposition(queue, context,
                    Ut->rowPointers, Ut->colIndices, Ut->nnzValues, Ut->nnzbs,
                    Lmat->rowPointers, Lmat->colIndices, Lmat->nnzValues, Lmat->nnzbs,
                    LUmat->rowPointers, LUmat->colIndices, LUmat->nnzValues, LUmat->nnzbs,
                    Nb, num_sweeps, verbosity);
#else
        double *Ltmp = new double[Lmat->nnzbs * block_size * block_size];
        for (int sweep = 0; sweep < num_sweeps; ++sweep) {

            // for every row
            for (int row = 0; row < Nb; row++) {
                // update U
                int jColStart = Ut->rowPointers[row];
                int jColEnd = Ut->rowPointers[row + 1];
                // for every block in this row
                for (int ij = jColStart; ij < jColEnd; ij++) {
                    int col = Ut->colIndices[ij];
                    // refine Uij element (or diagonal)
                    int i1 = LUmat->rowPointers[col];
                    int i2 = LUmat->rowPointers[col+1];
                    int kk = 0;
                    for(kk = i1; kk < i2; ++kk) {
                        ptrdiff_t c = LUmat->colIndices[kk];
                        if (c >= row) {
                            break;
                        }
                    }
                    double aij[bs*bs];
                    memcpy(&aij[0], LUmat->nnzValues + kk * bs * bs, sizeof(double) * bs * bs);
                    int jk = Lmat->rowPointers[col];
                    int ik = (jk < Lmat->rowPointers[col+1]) ? Lmat->colIndices[jk] : Nb;

                    for (int k = jColStart; k < ij; ++k) {
                        int ki = Ut->colIndices[k];
                        while (ik < ki) {
                            ++jk;
                            ik = Lmat->colIndices[jk];
                        }
                        if (ik == ki) {
                            blockMultSub<bs>(&aij[0], Lmat->nnzValues + jk * bs * bs, Ut->nnzValues + k * bs * bs);
                        }
                    }

                    memcpy(Utmp + ij * bs * bs, &aij[0], sizeof(double) * bs * bs);
                }

                // update L
                int iRowStart = Lmat->rowPointers[row];
                int iRowEnd = Lmat->rowPointers[row + 1];

                for (int ij = iRowStart; ij < iRowEnd; ij++) {
                    int j = Lmat->colIndices[ij];
                    // refine Lij element
                    int i1 = LUmat->rowPointers[row];
                    int i2 = LUmat->rowPointers[row+1];
                    int kk = 0;
                    for(kk = i1; kk < i2; ++kk) {
                        ptrdiff_t c = LUmat->colIndices[kk];
                        if (c >= j) {
                            break;
                        }
                    }
                    double aij[bs*bs];
                    memcpy(&aij[0], LUmat->nnzValues + kk * bs * bs, sizeof(double) * bs * bs);
                    int jk = Ut->rowPointers[j];
                    int ik = Ut->colIndices[jk];
                    for (int k = iRowStart; k < ij; ++k) {
                        int ki = Lmat->colIndices[k];
                        while(ik < ki) {
                            ++jk;
                            ik = Ut->colIndices[jk];
                        }

                        if(ik == ki) {
                            blockMultSub<bs>(&aij[0], Lmat->nnzValues + k * bs * bs , Ut->nnzValues + jk * bs * bs);
                        }
                    }
                    // calculate aij / ujj
                    double ujj[bs*bs];
                    inverter(Ut->nnzValues + (Ut->rowPointers[j+1] - 1) * bs * bs, &ujj[0]);
                    // lij = aij / ujj
                    blockMult<bs>(&aij[0], &ujj[0], Ltmp + ij * bs * bs);

                }
            }
            double *t = Lmat->nnzValues;
            Lmat->nnzValues = Ltmp;
            Ltmp = t;
            t = Ut->nnzValues;
            Ut->nnzValues = Utmp;
            Utmp = t;
        } // end sweep
        delete[] Ltmp;
#endif

        // convert Ut to BSR
        // diagonal stored separately
        std::vector<int> ptr(Nb+1, 0);
        std::vector<int> col(Ut->rowPointers[Nb]);

        // count blocks per row for U (BSR)
        // store diagonal in invDiagVals
        for(int i = 0; i < Nb; ++i) {
            for(int k = Ut->rowPointers[i]; k < Ut->rowPointers[i+1]; ++k) {
                int j = Ut->colIndices[k];
                if (j == i) {
                    inverter(Ut->nnzValues + k * bs * bs, invDiagVals + i * bs * bs);
                } else {
                    ++ptr[j+1];
                }
            }
        }

        // prefix sum
        int sumU = 0;
        for (int i = 1; i < Nb+1; i++) {
            sumU += ptr[i];
            ptr[i] = sumU;
        }

        // actually copy nonzero values for U
        for(int i = 0; i < Nb; ++i) {
            for(int k = Ut->rowPointers[i]; k < Ut->rowPointers[i+1]; ++k) {
                int j = Ut->colIndices[k];
                if (j != i) {
                    int head = ptr[j]++;
                    col[head]  = i;
                    memcpy(Utmp + head * bs * bs, Ut->nnzValues + k * bs * bs, sizeof(double) * bs * bs);
                }
            }
        }

        // the ptr[] were increased in the last loop
        std::rotate(ptr.begin(), ptr.end() - 1, ptr.end());
        ptr.front() = 0;

        // reversing the rows of U, because that is the order they are used in
        int URowIndex = 0;
        int offsetU = 0;   // number of nnz blocks that are already copied to Umat
        Umat->rowPointers[0] = 0;
        for (int i = LUmat->Nb - 1; i >= 0; i--) {
            int rowSize = ptr[i + 1] - ptr[i];   // number of blocks in this row
            memcpy(Umat->nnzValues + offsetU * bs * bs, Utmp + ptr[i] * bs * bs, sizeof(double) * bs * bs * rowSize);
            memcpy(Umat->colIndices + offsetU, col.data() + ptr[i], sizeof(int) * rowSize);
            offsetU += rowSize;
            Umat->rowPointers[URowIndex + 1] = offsetU;
            URowIndex++;
        }

        delete[] Utmp;

    }

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

        Timer t_decomposition;
#if CHOW_PATEL
        chow_patel_decomposition();
#else
        int i, j, ij, ik, jk;
        int iRowStart, iRowEnd, jRowEnd;
        double pivot[bs * bs];

        int LSize = 0;
        Opm::Detail::Inverter<bs> inverter;   // reuse inverter to invert blocks

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
#endif
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
template BILU0<n>::BILU0(ILUReorder, int);                                               \
template BILU0<n>::~BILU0();                                                             \
template bool BILU0<n>::init(BlockedMatrix<n>*);                                         \
template void BILU0<n>::chow_patel_decomposition();                                      \
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


