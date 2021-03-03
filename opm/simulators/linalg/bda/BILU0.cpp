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
#include <opm/simulators/linalg/bda/ChowPatelIlu.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


namespace bda
{

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
        if (opencl_ilu_reorder != ILUReorder::NONE) {
            delete[] toOrder;
            delete[] fromOrder;
        }
    }

    template <unsigned int block_size>
    bool BILU0<block_size>::init(BlockedMatrix<block_size> *mat)
    {
        const unsigned int bs = block_size;

        this->N = mat->Nb * block_size;
        this->Nb = mat->Nb;
        this->nnz = mat->nnzbs * block_size * block_size;
        this->nnzbs = mat->nnzbs;

        int *CSCRowIndices = nullptr;
        int *CSCColPointers = nullptr;

        if (opencl_ilu_reorder == ILUReorder::NONE) {
            LUmat = std::make_unique<BlockedMatrix<block_size> >(*mat);
        } else {
            toOrder = new int[Nb];
            fromOrder = new int[Nb];
            CSCRowIndices = new int[nnzbs];
            CSCColPointers = new int[Nb + 1];
            rmat = std::make_shared<BlockedMatrix<block_size> >(mat->Nb, mat->nnzbs);
            LUmat = std::make_unique<BlockedMatrix<block_size> >(*rmat);

            Timer t_convert;
            csrPatternToCsc(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb);
            if(verbosity >= 3){
                std::ostringstream out;
                out << "BILU0 convert CSR to CSC: " << t_convert.stop() << " s";
                OpmLog::info(out.str());
            }
        }

        Timer t_analysis;
        std::ostringstream out;
        if (opencl_ilu_reorder == ILUReorder::LEVEL_SCHEDULING) {
            out << "BILU0 reordering strategy: " << "level_scheduling\n";
            findLevelScheduling(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        } else if (opencl_ilu_reorder == ILUReorder::GRAPH_COLORING) {
            out << "BILU0 reordering strategy: " << "graph_coloring\n";
            findGraphColoring<block_size>(mat->colIndices, mat->rowPointers, CSCRowIndices, CSCColPointers, mat->Nb, mat->Nb, mat->Nb, &numColors, toOrder, fromOrder, rowsPerColor);
        } else if (opencl_ilu_reorder == ILUReorder::NONE) {
            out << "BILU0 reordering strategy: none\n";
            // numColors = 1;
            // rowsPerColor.emplace_back(Nb);
            numColors = Nb;
            for(int i = 0; i < Nb; ++i){
                rowsPerColor.emplace_back(1);
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

        if (opencl_ilu_reorder != ILUReorder::NONE) {
            delete[] CSCRowIndices;
            delete[] CSCColPointers;
        }

        diagIndex = new int[mat->Nb];
        invDiagVals = new double[mat->Nb * bs * bs];

#if CHOW_PATEL
        Lmat = std::make_unique<BlockedMatrix<block_size> >(mat->Nb, (mat->nnzbs - mat->Nb) / 2);
        Umat = std::make_unique<BlockedMatrix<block_size> >(mat->Nb, (mat->nnzbs - mat->Nb) / 2);
#endif

        LUmat->nnzValues = new double[mat->nnzbs * bs * bs];

        s.invDiagVals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * bs * bs * mat->Nb);
        s.rowsPerColor = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (numColors + 1));
        s.diagIndex = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * LUmat->Nb);
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

        events.resize(2);
        err = queue->enqueueWriteBuffer(s.invDiagVals, CL_FALSE, 0, mat->Nb * sizeof(double) * bs * bs, invDiagVals, nullptr, &events[0]);

        rowsPerColorPrefix.resize(numColors + 1); // resize initializes value 0.0
        for (int i = 0; i < numColors; ++i) {
            rowsPerColorPrefix[i+1] = rowsPerColorPrefix[i] + rowsPerColor[i];
        }
        err |= queue->enqueueWriteBuffer(s.rowsPerColor, CL_FALSE, 0, (numColors + 1) * sizeof(int), rowsPerColorPrefix.data(), nullptr, &events[1]);

        cl::WaitForEvents(events);
        events.clear();
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "BILU0 OpenCL enqueueWriteBuffer error");
        }

        return true;
    } // end init()


    // implements Fine-Grained Parallel ILU algorithm (FGPILU), Chow and Patel 2015
    template <unsigned int block_size>
    void BILU0<block_size>::chow_patel_decomposition()
    {
#if CHOW_PATEL
        const unsigned int bs = block_size;
        int num_sweeps = 6;

        // split matrix into L and U
        // also convert U into BSC format (Ut)
        // Ut stores diagonal for now
        // original matrix LUmat is assumed to be symmetric

#ifndef NDEBUG
        // verify that matrix is symmetric
        for (int i = 0; i < Nb; ++i){
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block (i, j) in this row, check if (j, i) also exists
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                int jRowStart = LUmat->rowPointers[j];
                int jRowEnd = LUmat->rowPointers[j + 1];
                bool blockFound = false;
                // check all blocks on row j
                // binary search is possible
                for (int ji = jRowStart; ji < jRowEnd; ji++) {
                    int row = LUmat->colIndices[ji];
                    if (i == row) {
                        blockFound = true;
                        break;
                    }
                }
                if (false == blockFound) {
                    OPM_THROW(std::logic_error, "Error sparsity pattern must be symmetric when using chow_patel_decomposition()");
                }
            }
        }
#endif

        Timer t_total, t_preprocessing;

        // Ut is actually BSC format
        std::unique_ptr<BlockedMatrix<bs> > Ut = std::make_unique<BlockedMatrix<bs> >(Nb, (nnzbs + Nb) / 2);

        Lmat->rowPointers[0] = 0;
        for (int i = 0; i < Nb+1; i++) {
            Ut->rowPointers[i] = 0;
        }

        Opm::Detail::Inverter<bs> inverter;

        // store inverted diagonal
        for (int i = 0; i < Nb; i++) {
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block in this row
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                if (i == j) {
                    inverter(LUmat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs);
                }
            }
        }

        // initialize initial guess for L: L_A * D
        // L_A is strictly lower triangular part of A
        // D is inv(diag(A))
        int num_blocks_L = 0;
        for (int i = 0; i < Nb; i++) {
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block in this row
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                if (i <= j) {
                    Ut->rowPointers[j+1]++;   // actually colPointers, now simply indicates how many blocks this col holds
                } else {
                    Lmat->colIndices[num_blocks_L] = j;
                    // multiply block of L with corresponding diag block
                    blockMult<bs>(LUmat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs, Lmat->nnzValues + num_blocks_L * bs * bs);
                    num_blocks_L++;
                }
            }
            // TODO: copy all blocks for L at once, instead of copying each block individually
            Lmat->rowPointers[i+1] = num_blocks_L;
        }

        // prefix sum to sum rowsizes into colpointers
        std::partial_sum(Ut->rowPointers, Ut->rowPointers+Nb+1, Ut->rowPointers);

        // initialize initial guess for U
        for (int i = 0; i < Nb; i++) {
            int iRowStart = LUmat->rowPointers[i];
            int iRowEnd = LUmat->rowPointers[i + 1];
            // for every block in this row
            for (int ij = iRowStart; ij < iRowEnd; ij++) {
                int j = LUmat->colIndices[ij];
                if (i <= j){
                    int idx = Ut->rowPointers[j]++; // rowPointers[i] is (mis)used as the write offset of the current row i
                    Ut->colIndices[idx] = i;     // actually rowIndices
                    memcpy(Ut->nnzValues + idx * bs * bs, LUmat->nnzValues + ij * bs * bs, sizeof(double) * bs * bs);
                }
            }
        }

        // rotate
        // the Ut->rowPointers were increased in the last loop
        // now Ut->rowPointers[i+1] is at the same position as Ut->rowPointers[i] should have for a crs matrix. reset to correct expected value
        for (int i = Nb; i > 0; --i) {
            Ut->rowPointers[i] = Ut->rowPointers[i-1];
        }
        Ut->rowPointers[0] = 0;


        // Utmp is needed for CPU and GPU decomposition, because U is transposed, and reversed after decomposition
        // U will be reversed because it is used with backwards substitution, the last row is used first
        // Ltmp is only needed for CPU decomposition, GPU creates GPU buffer for Ltmp
        double *Utmp = new double[Ut->nnzbs * block_size * block_size];

        if (verbosity >= 3) {
            std::ostringstream out;
            out << "BILU0 ChowPatel preprocessing: " << t_preprocessing.stop() << " s";
            OpmLog::info(out.str());
        }

        // actual ILU decomposition
        Timer t_decomposition;
#if CHOW_PATEL_GPU
        chowPatelIlu.decomposition(queue, context,
                    Ut->rowPointers, Ut->colIndices, Ut->nnzValues, Ut->nnzbs,
                    Lmat->rowPointers, Lmat->colIndices, Lmat->nnzValues, Lmat->nnzbs,
                    LUmat->rowPointers, LUmat->colIndices, LUmat->nnzValues, LUmat->nnzbs,
                    Nb, num_sweeps, verbosity);
#else
        double *Ltmp = new double[Lmat->nnzbs * block_size * block_size];
        for (int sweep = 0; sweep < num_sweeps; ++sweep) {

            // algorithm
            // for every block in A (LUmat):
            //     if i > j:
            //         Lij = (Aij - sum k=1 to j-1 {Lik*Ukj}) / Ujj
            //     else:
            //         Uij = (Aij - sum k=1 to i-1 {Lik*Ukj})

            // for every row
            for (int row = 0; row < Nb; row++) {
                // update U
                // Uij = (Aij - sum k=1 to i-1 {Lik*Ukj})
                int jColStart = Ut->rowPointers[row];
                int jColEnd = Ut->rowPointers[row + 1];
                int colU = row; // rename for clarity, next row in Ut means next col in U
                // for every block in this row
                for (int ij = jColStart; ij < jColEnd; ij++) {
                    int rowU1 = Ut->colIndices[ij]; // actually rowIndices for U
                    // refine Uij element (or diagonal)
                    int i1 = LUmat->rowPointers[rowU1];
                    int i2 = LUmat->rowPointers[rowU1+1];

                    // search on row rowU1, find blockIndex in LUmat of block with same col (colU) as Uij
                    // LUmat->nnzValues[kk] is block Aij
                    auto candidate = std::find(LUmat->colIndices + i1, LUmat->colIndices + i2, colU);
                    assert(candidate != LUmat->colIndices + i2);
                    auto kk = candidate - LUmat->colIndices;

                    double aij[bs*bs];
                    // copy block to Aij so operations can be done on it without affecting LUmat
                    memcpy(&aij[0], LUmat->nnzValues + kk * bs * bs, sizeof(double) * bs * bs);

                    int jk = Lmat->rowPointers[rowU1]; // points to row rowU1 in L
                    // if row rowU1 is empty, skip row
                    if (jk < Lmat->rowPointers[rowU1+1]) {
                        int colL = Lmat->colIndices[jk];
                        // only check until block U(i,j) is reached
                        for (int k = jColStart; k < ij; ++k) {
                            int rowU2 = Ut->colIndices[k];
                            while (colL < rowU2) {
                                ++jk; // check next block on row rowU1 of L
                                colL = Lmat->colIndices[jk];
                            }
                            if (colL == rowU2) {
                                // Aij -= (Lik * Ukj)
                                blockMultSub<bs>(&aij[0], Lmat->nnzValues + jk * bs * bs, Ut->nnzValues + k * bs * bs);
                            }
                        }
                    }

                    // Uij_new = Aij - sum
                    memcpy(Utmp + ij * bs * bs, &aij[0], sizeof(double) * bs * bs);
                }

                // update L
                // Lij = (Aij - sum k=1 to j-1 {Lik*Ukj}) / Ujj
                int iRowStart = Lmat->rowPointers[row];
                int iRowEnd = Lmat->rowPointers[row + 1];

                for (int ij = iRowStart; ij < iRowEnd; ij++) {
                    int j = Lmat->colIndices[ij];
                    // refine Lij element
                    // search on row 'row', find blockIndex in LUmat of block with same col (j) as Lij
                    // LUmat->nnzValues[kk] is block Aij
                    int i1 = LUmat->rowPointers[row];
                    int i2 = LUmat->rowPointers[row+1];

                    auto candidate = std::find(LUmat->colIndices + i1, LUmat->colIndices + i2, j);
                    assert(candidate != LUmat->colIndices + i2);
                    auto kk = candidate - LUmat->colIndices;

                    double aij[bs*bs];
                    // copy block to Aij so operations can be done on it without affecting LUmat
                    memcpy(&aij[0], LUmat->nnzValues + kk * bs * bs, sizeof(double) * bs * bs);

                    int jk = Ut->rowPointers[j];  // actually colPointers, jk points to col j in U
                    int rowU = Ut->colIndices[jk];  // actually rowIndices, rowU is the row of block jk
                    // only check until block L(i,j) is reached
                    for (int k = iRowStart; k < ij; ++k) {
                        int colL = Lmat->colIndices[k];
                        while (rowU < colL) {
                            ++jk; // check next block on col j of U
                            rowU = Ut->colIndices[jk];
                        }

                        if (rowU == colL) {
                            // Aij -= (Lik * Ukj)
                            blockMultSub<bs>(&aij[0], Lmat->nnzValues + k * bs * bs , Ut->nnzValues + jk * bs * bs);
                        }
                    }

                    // calculate (Aij - sum) / Ujj
                    double ujj[bs*bs];
                    inverter(Ut->nnzValues + (Ut->rowPointers[j+1] - 1) * bs * bs, &ujj[0]);
                    // Lij_new = (Aij - sum) / Ujj
                    blockMult<bs>(&aij[0], &ujj[0], Ltmp + ij * bs * bs);
                }
            }
            // 1st sweep writes to Ltmp
            // 2nd sweep writes to Lmat->nnzValues
            std::swap(Lmat->nnzValues, Ltmp);
            std::swap(Ut->nnzValues, Utmp);
        } // end sweep

        // if number of sweeps is even, swap again so data is in Lmat->nnzValues
        if (num_sweeps % 2 == 0) {
            std::swap(Lmat->nnzValues, Ltmp);
            std::swap(Ut->nnzValues, Utmp);
        }
        delete[] Ltmp;
#endif

        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 ChowPatel decomposition: " << t_decomposition.stop() << " s";
            OpmLog::info(out.str());
        }

        Timer t_postprocessing;

        // convert Ut to BSR
        // diagonal stored separately
        std::vector<int> ptr(Nb+1, 0);
        std::vector<int> cols(Ut->rowPointers[Nb]);

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
        std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

        // actually copy nonzero values for U
        for(int i = 0; i < Nb; ++i) {
            for(int k = Ut->rowPointers[i]; k < Ut->rowPointers[i+1]; ++k) {
                int j = Ut->colIndices[k];
                if (j != i) {
                    int head = ptr[j]++;
                    cols[head]  = i;
                    memcpy(Utmp + head * bs * bs, Ut->nnzValues + k * bs * bs, sizeof(double) * bs * bs);
                }
            }
        }

        // the ptr[] were increased in the last loop
        std::rotate(ptr.begin(), ptr.end() - 1, ptr.end());
        ptr.front() = 0;


        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 ChowPatel postprocessing: " << t_postprocessing.stop() << " s";
            OpmLog::info(out.str());
        }

        Timer t_copyToGpu;

        events.resize(3);
        queue->enqueueWriteBuffer(s.Lvals, CL_FALSE, 0, Lmat->nnzbs * bs * bs * sizeof(double), Lmat->nnzValues, nullptr, &events[0]);
        queue->enqueueWriteBuffer(s.Uvals, CL_FALSE, 0, Umat->nnzbs * bs * bs * sizeof(double), Utmp, nullptr, &events[1]);
        queue->enqueueWriteBuffer(s.invDiagVals, CL_FALSE, 0, LUmat->Nb * bs * bs * sizeof(double), invDiagVals, nullptr, &events[2]);

        std::call_once(pattern_uploaded, [&](){
            // find the positions of each diagonal block
            // must be done after reordering
            for (int row = 0; row < Nb; ++row) {
                int rowStart = LUmat->rowPointers[row];
                int rowEnd = LUmat->rowPointers[row+1];

                auto candidate = std::find(LUmat->colIndices + rowStart, LUmat->colIndices + rowEnd, row);
                assert(candidate != LUmat->colIndices + rowEnd);
                diagIndex[row] = candidate - LUmat->colIndices;
            }
            events.resize(8);
            queue->enqueueWriteBuffer(s.diagIndex, CL_FALSE, 0, Nb * sizeof(int), diagIndex, nullptr, &events[3]);
            queue->enqueueWriteBuffer(s.Lcols, CL_FALSE, 0, Lmat->nnzbs * sizeof(int), Lmat->colIndices, nullptr, &events[4]);
            queue->enqueueWriteBuffer(s.Lrows, CL_FALSE, 0, (Lmat->Nb + 1) * sizeof(int), Lmat->rowPointers, nullptr, &events[5]);
            queue->enqueueWriteBuffer(s.Ucols, CL_FALSE, 0, Umat->nnzbs * sizeof(int), cols.data(), nullptr, &events[6]);
            queue->enqueueWriteBuffer(s.Urows, CL_FALSE, 0, (Umat->Nb + 1) * sizeof(int), ptr.data(), nullptr, &events[7]);
        });

        cl::WaitForEvents(events);
        events.clear();
        if (err != CL_SUCCESS) {
            // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
            OPM_THROW(std::logic_error, "BILU0 OpenCL enqueueWriteBuffer error");
        }

        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 ChowPatel copy to GPU: " << t_copyToGpu.stop() << " s\n";
            out << "BILU0 ChowPatel total: " << t_total.stop() << " s";
            OpmLog::info(out.str());
        }

        delete[] Utmp;
#endif // CHOW_PATEL
    }

    template <unsigned int block_size>
    bool BILU0<block_size>::create_preconditioner(BlockedMatrix<block_size> *mat)
    {
        const unsigned int bs = block_size;
        auto *m = mat;

        if (opencl_ilu_reorder != ILUReorder::NONE) {
            m = rmat.get();
            Timer t_reorder;
            reorderBlockedMatrixByPattern<block_size>(mat, toOrder, fromOrder, rmat.get());

            if (verbosity >= 3){
                std::ostringstream out;
                out << "BILU0 reorder matrix: " << t_reorder.stop() << " s";
                OpmLog::info(out.str());
            }
        }

        // TODO: remove this copy by replacing inplace ilu decomp by out-of-place ilu decomp
        // this copy can have mat or rmat ->nnzValues as origin, depending on the reorder strategy
        Timer t_copy;
        memcpy(LUmat->nnzValues, m->nnzValues, sizeof(double) * bs * bs * m->nnzbs);

        if (verbosity >= 3){
            std::ostringstream out;
            out << "BILU0 memcpy: " << t_copy.stop() << " s";
            OpmLog::info(out.str());
        }

#if CHOW_PATEL
        chow_patel_decomposition();
#else
        Timer t_copyToGpu;

        events.resize(1);
        queue->enqueueWriteBuffer(s.LUvals, CL_FALSE, 0, LUmat->nnzbs * bs * bs * sizeof(double), LUmat->nnzValues, nullptr, &events[0]);

        std::call_once(pattern_uploaded, [&](){
            // find the positions of each diagonal block
            // must be done after reordering
            for (int row = 0; row < Nb; ++row) {
                int rowStart = LUmat->rowPointers[row];
                int rowEnd = LUmat->rowPointers[row+1];

                auto candidate = std::find(LUmat->colIndices + rowStart, LUmat->colIndices + rowEnd, row);
                assert(candidate != LUmat->colIndices + rowEnd);
                diagIndex[row] = candidate - LUmat->colIndices;
            }
            events.resize(4);
            queue->enqueueWriteBuffer(s.diagIndex, CL_FALSE, 0, Nb * sizeof(int), diagIndex, nullptr, &events[1]);
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
        cl::Event event;
        for (int color = 0; color < numColors; ++color) {
            const unsigned int firstRow = rowsPerColorPrefix[color];
            const unsigned int lastRow = rowsPerColorPrefix[color+1];
            const unsigned int work_group_size2 = 128;
            const unsigned int num_work_groups2 = 1024;
            const unsigned int total_work_items2 = num_work_groups2 * work_group_size2;
            const unsigned int num_hwarps_per_group = work_group_size2 / 16;
            const unsigned int lmem_per_work_group2 = num_hwarps_per_group * bs * bs * sizeof(double);           // each block needs a pivot
            if (verbosity >= 4) {
                out << "color " << color << ": " << firstRow << " - " << lastRow << " = " << lastRow - firstRow << "\n";
            }
            event = (*ilu_decomp_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items2), cl::NDRange(work_group_size2)), firstRow, lastRow, s.LUvals, s.LUcols, s.LUrows, s.invDiagVals, s.diagIndex, LUmat->Nb, cl::Local(lmem_per_work_group2));
            event.wait();
        }

        if (verbosity >= 3) {
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
    void BILU0<block_size>::apply(cl::Buffer& x, cl::Buffer& y)
    {
        cl::Event event;
        Timer t_apply;

        for(int color = 0; color < numColors; ++color){
#if CHOW_PATEL
            event = (*ILU_apply1)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.Lvals, s.Lcols, s.Lrows, s.diagIndex, x, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
#else
            event = (*ILU_apply1)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.LUvals, s.LUcols, s.LUrows, s.diagIndex, x, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
#endif
            // event.wait();
        }

        for(int color = numColors-1; color >= 0; --color){
#if CHOW_PATEL
            event = (*ILU_apply2)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.Uvals, s.Ucols, s.Urows, s.diagIndex, s.invDiagVals, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
#else
            event = (*ILU_apply2)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), s.LUvals, s.LUcols, s.LUrows, s.diagIndex, s.invDiagVals, y, s.rowsPerColor, color, block_size, cl::Local(lmem_per_work_group));
#endif
            // event.wait();
        }

        if (verbosity >= 4) {
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
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply1_,
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply2_,
        cl::make_kernel<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg> *ilu_decomp_k_
    ){
        this->ILU_apply1 = ILU_apply1_;
        this->ILU_apply2 = ILU_apply2_;
        this->ilu_decomp_k = ilu_decomp_k_;
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
    cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *, \
    cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *, \
    cl::make_kernel<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg> *                 \
    );

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // end namespace bda


