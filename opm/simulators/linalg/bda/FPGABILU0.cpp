/*
  Copyright 2020 Equinor ASA

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

#include <opm/simulators/linalg/bda/FPGABILU0.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>
#include <opm/simulators/linalg/bda/FPGAUtils.hpp>

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
FPGABILU0<block_size>::FPGABILU0(ILUReorder opencl_ilu_reorder_, int verbosity_, int maxRowsPerColor_, int maxColsPerColor_, int maxNNZsPerRow_, int maxNumColors_) :
    verbosity(verbosity_), opencl_ilu_reorder(opencl_ilu_reorder_), maxRowsPerColor(maxRowsPerColor_), maxColsPerColor(maxColsPerColor_), maxNNZsPerRow(maxNNZsPerRow_), maxNumColors(maxNumColors_)
{
    if (opencl_ilu_reorder == ILUReorder::LEVEL_SCHEDULING) {
        level_scheduling = true;
    } else if (opencl_ilu_reorder == ILUReorder::GRAPH_COLORING) {
        graph_coloring = true;
    } else {
        OPM_THROW(std::logic_error, "Error ilu reordering strategy not set correctly\n");
    }
}


template <unsigned int block_size>
FPGABILU0<block_size>::~FPGABILU0()
{
    delete[] invDiagVals;
}


template <unsigned int block_size>
bool FPGABILU0<block_size>::init(BlockedMatrix *mat)
{
    const unsigned int bs = block_size;

    resultPointers.resize(numResultPointers, nullptr);
    resultSizes.resize(numResultSizes);

    // Set nnzSplit as hardcoded constant until support for more than one nnzVals read array is added.
    const unsigned int nnzSplit = 1;

    this->N = mat->Nb * block_size;
    this->Nb = mat->Nb;
    this->nnz = mat->nnzbs * block_size * block_size;
    this->nnzbs = mat->nnzbs;

    toOrder.resize(Nb);
    fromOrder.resize(Nb);

    std::vector<int> CSCRowIndices(nnzbs);
    std::vector<int> CSCColPointers(Nb + 1);

    if (level_scheduling) {
        Timer t_convert;
        csrPatternToCsc(mat->colIndices, mat->rowPointers, CSCRowIndices.data(), CSCColPointers.data(), mat->Nb);
        if (verbosity >= 3) {
            std::ostringstream out;
            out << "FPGABILU0 convert CSR to CSC: " << t_convert.stop() << " s";
            OpmLog::info(out.str());
        }
    }

    Timer t_analysis;
    rMat = std::make_shared<BlockedMatrix>(mat->Nb, mat->nnzbs, block_size);
    LUMat = std::make_unique<BlockedMatrix>(*rMat);
    std::ostringstream out;
    if (level_scheduling) {
        out << "FPGABILU0 reordering strategy: " << "level_scheduling\n";
        findLevelScheduling(mat->colIndices, mat->rowPointers, CSCRowIndices.data(), CSCColPointers.data(), mat->Nb, &numColors, toOrder.data(), fromOrder.data(), rowsPerColor);
    } else if (graph_coloring) {
        out << "FPGABILU0 reordering strategy: " << "graph_coloring\n";
        findGraphColoring<bs>(mat->colIndices, mat->rowPointers, CSCRowIndices.data(), CSCColPointers.data(), mat->Nb, maxRowsPerColor, maxColsPerColor, &numColors, toOrder.data(), fromOrder.data(), rowsPerColor);
    }

    if (numColors > maxNumColors) {
        std::ostringstream errorstring;
        errorstring << "ERROR: the matrix was reordered into too many colors. Created " << numColors << " colors, while hardware only supports up to " << maxNumColors << "\n";
        OPM_THROW(std::logic_error, errorstring.str());
    }

    if (verbosity >= 3) {
        out << "FPGABILU0 analysis took: " << t_analysis.stop() << " s, " << numColors << " colors";
    }
    OpmLog::info(out.str());

    int colorRoundedValSize = 0, LColorRoundedValSize = 0, UColorRoundedValSize = 0;
    int NROffsetSize = 0, LNROffsetSize = 0, UNROffsetSize = 0;
    int blockDiagSize = 0;
    // This reordering is needed here only to te result can be used to calculate worst-case scenario array sizes
    reorderBlockedMatrixByPattern(mat, toOrder.data(), fromOrder.data(), rMat.get());
    int doneRows = 0;
    for (int c = 0; c < numColors; c++) {
        for (int i = doneRows; i < doneRows + rowsPerColor[c]; i++) {
            for (int j = rMat->rowPointers[i]; j < rMat->rowPointers[i + 1]; j++) {
                int columnIndex = rMat->colIndices[j];
                if (columnIndex < i) {
                    LColorRoundedValSize += 9;
                    LNROffsetSize += 9;
                }
                if (columnIndex > i) {
                    UColorRoundedValSize += 9;
                    UNROffsetSize += 9;
                }
                colorRoundedValSize += 9;
                NROffsetSize += 9;
            }
            blockDiagSize += 12;
        }
        // End of color: round all sizes to nearest cacheline
        colorRoundedValSize = roundUpTo(colorRoundedValSize, 32);
        LColorRoundedValSize = roundUpTo(LColorRoundedValSize, 32);
        UColorRoundedValSize = roundUpTo(UColorRoundedValSize, 32);
        NROffsetSize = roundUpTo(NROffsetSize, 64);
        LNROffsetSize = roundUpTo(LNROffsetSize, 64);
        UNROffsetSize = roundUpTo(UNROffsetSize, 64);
        blockDiagSize = roundUpTo(blockDiagSize, 8);

        doneRows += rowsPerColor[c];
    }
    int colorSizesNum = 8 + roundUpTo(4 * numColors, 16);
    int worstCaseColumnAccessNum = numColors * maxColsPerColor;

    nnzValues.resize(nnzSplit, std::vector<double>(colorRoundedValSize));
    LnnzValues.resize(nnzSplit, std::vector<double>(LColorRoundedValSize));
    UnnzValues.resize(nnzSplit, std::vector<double>(UColorRoundedValSize));
    // initial number of nnz, used to allocate
    nnzValsSizes.resize(nnzSplit, colorRoundedValSize);
    LnnzValsSizes.resize(nnzSplit, LColorRoundedValSize);
    UnnzValsSizes.resize(nnzSplit, UColorRoundedValSize);
    colIndices.resize(colorRoundedValSize);
    LColIndices.resize(LColorRoundedValSize);
    UColIndices.resize(UColorRoundedValSize);
    NROffsets.resize(NROffsetSize);
    LNROffsets.resize(LNROffsetSize);
    UNROffsets.resize(UNROffsetSize);
    PIndicesAddr.resize(worstCaseColumnAccessNum);
    LPIndicesAddr.resize(worstCaseColumnAccessNum);
    UPIndicesAddr.resize(worstCaseColumnAccessNum);
    colorSizes.resize(colorSizesNum);
    LColorSizes.resize(colorSizesNum);
    UColorSizes.resize(colorSizesNum);
    blockDiag.resize(blockDiagSize);
    colIndicesInColor.resize(numColors, std::vector<int>(rMat->Nb * block_size, 0));
    LColIndicesInColor.resize(numColors, std::vector<int>(rMat->Nb * block_size, 0));
    UColIndicesInColor.resize(numColors, std::vector<int>(rMat->Nb * block_size, 0));

    int err = rMat->findPartitionColumns(numColors, rowsPerColor.data(),
                                         maxRowsPerColor, maxColsPerColor,
                                         colIndicesInColor, PIndicesAddr.data(), colorSizes.data(),
                                         LColIndicesInColor, LPIndicesAddr.data(), LColorSizes.data(),
                                         UColIndicesInColor, UPIndicesAddr.data(), UColorSizes.data());
    if (err != 0) {
        std::ostringstream errorstring;
        errorstring << "ERROR: findPartitionColumns failed, code " << err << "\n";
        OPM_THROW(std::logic_error, errorstring.str());
    }

    diagIndex.resize(mat->Nb, 0);
    invDiagVals = new double[mat->Nb * bs * bs];
    LMat = std::make_unique<BlockedMatrix>(mat->Nb, (mat->nnzbs - mat->Nb) / 2, block_size);
    UMat = std::make_unique<BlockedMatrix>(mat->Nb, (mat->nnzbs - mat->Nb) / 2, block_size);
    resultPointers[0] = (void *) colorSizes.data();
    resultPointers[1] = (void *) PIndicesAddr.data();
    resultPointers[2] = (void *) nnzValues.data();
    resultPointers[3] = (void *) colIndices.data();
    resultPointers[4] = (void *) NROffsets.data();
    resultPointers[5] = (void *) nnzValsSizes.data();
    resultPointers[6] = (void *) LColorSizes.data();
    resultPointers[7] = (void *) LPIndicesAddr.data();
    resultPointers[8] = (void *) LnnzValues.data();
    resultPointers[9] = (void *) LColIndices.data();
    resultPointers[10] = (void *) LNROffsets.data();
    resultPointers[11] = (void *) LnnzValsSizes.data();
    resultPointers[12] = (void *) UColorSizes.data();
    resultPointers[13] = (void *) UPIndicesAddr.data();
    resultPointers[14] = (void *) UnnzValues.data();
    resultPointers[15] = (void *) UColIndices.data();
    resultPointers[16] = (void *) UNROffsets.data();
    resultPointers[17] = (void *) UnnzValsSizes.data();
    resultPointers[18] = (void *) blockDiag.data();
    //resultPointers[19] and [20] are set by the caller
    resultSizes[0] = mat->Nb * block_size;
    resultSizes[1] = colorRoundedValSize; // zeropadded valSize;
    resultSizes[2] = numColors;
    resultSizes[3] = worstCaseColumnAccessNum; //totalCols
    resultSizes[4] = NROffsetSize; //NRFlagSize
    resultSizes[5] = blockDiagSize; //diagValsSize
    resultSizes[6] = mat->Nb * block_size;
    resultSizes[7] = LColorRoundedValSize; // zeropadded LValSize;
    resultSizes[8] = numColors;
    resultSizes[9] = worstCaseColumnAccessNum; //LTotalCols
    resultSizes[10] = LNROffsetSize; //LNRFlagSize
    resultSizes[11] = blockDiagSize; //LDiagValsSize
    resultSizes[12] = mat->Nb * block_size;
    resultSizes[13] = UColorRoundedValSize; // zeropadded UValSize;
    resultSizes[14] = numColors;
    resultSizes[15] = worstCaseColumnAccessNum; //UTotalCols
    resultSizes[16] = UNROffsetSize; //UNRFlagSize
    resultSizes[17] = blockDiagSize; //UDiagValsSize
    return true;
} // end init()


template <unsigned int block_size>
bool FPGABILU0<block_size>::create_preconditioner(BlockedMatrix *mat)
{
    const unsigned int bs = block_size;
    Timer t_reorder;
    reorderBlockedMatrixByPattern(mat, toOrder.data(), fromOrder.data(), rMat.get());

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "FPGABILU0 reorder matrix: " << t_reorder.stop() << " s";
        OpmLog::info(out.str());
    }

    // TODO: remove this copy by replacing inplace ilu decomp by out-of-place ilu decomp
    Timer t_memcpy;
    memcpy(LUMat->nnzValues, rMat->nnzValues, sizeof(double) * bs * bs * rMat->nnzbs);

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "FPGABILU0 memcpy: " << t_memcpy.stop() << " s";
        OpmLog::info(out.str());
    }

    int i, j, ij, ik, jk;
    int iRowStart, iRowEnd, jRowEnd;
    double pivot[bs * bs];
    int LSize = 0;
    Opm::Detail::Inverter<bs> inverter;   // reuse inverter to invert blocks

    Timer t_decomposition;

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
            blockMult<bs>(LUMat->nnzValues + ij * bs * bs, invDiagVals + j * bs * bs, &pivot[0]);

            memcpy(LUMat->nnzValues + ij * bs * bs, &pivot[0], sizeof(double) * bs * bs);

            jRowEnd = LUMat->rowPointers[j + 1];
            jk = diagIndex[j] + 1;
            ik = ij + 1;
            // subtract that row scaled by the pivot from this row.
            while (ik < iRowEnd && jk < jRowEnd) {
                if (LUMat->colIndices[ik] == LUMat->colIndices[jk]) {
                    blockMultSub<bs>(LUMat->nnzValues + ik * bs * bs, pivot, LUMat->nnzValues + jk * bs * bs);
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
        inverter(LUMat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs);
        memcpy(LUMat->nnzValues + ij * bs * bs, invDiagVals + i * bs * bs, sizeof(double) * bs * bs);
    }

    LMat->rowPointers[0] = 0;
    UMat->rowPointers[0] = 0;

    // Split the LU matrix into two by comparing column indices to diagonal indices
    for (i = 0; i < LUMat->Nb; i++) {
        LMat->rowPointers[i + 1] = LMat->rowPointers[i];
        for (j = LUMat->rowPointers[i]; j < LUMat->rowPointers[i + 1]; j++) {
            if (j < diagIndex[i]) {
                memcpy(LMat->nnzValues + (LMat->rowPointers[i + 1]) * bs * bs, LUMat->nnzValues + j * bs * bs, sizeof(double) * bs * bs);
                LMat->colIndices[LMat->rowPointers[i + 1]] = LUMat->colIndices[j];
                LMat->rowPointers[i + 1] = LMat->rowPointers[i + 1] + 1;
            }
        }
    }
    // Reverse the order or the (blocked) rows for the U matrix,
    // because the rows are accessed in reverse order when applying the ILU0
    int URowIndex = 0;
    for (i = LUMat->Nb - 1; i >= 0; i--) {
        UMat->rowPointers[URowIndex + 1] = UMat->rowPointers[URowIndex];
        for (j = LUMat->rowPointers[i]; j < LUMat->rowPointers[i + 1]; j++) {
            if (j > diagIndex[i]) {
                memcpy(UMat->nnzValues + (UMat->rowPointers[URowIndex + 1]) * bs * bs, LUMat->nnzValues + j * bs * bs, sizeof(double) * bs * bs);
                UMat->colIndices[UMat->rowPointers[URowIndex + 1]] = LUMat->colIndices[j];
                UMat->rowPointers[URowIndex + 1] = UMat->rowPointers[URowIndex + 1] + 1;
            }
        }
        URowIndex++;
    }

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "FPGABILU0 decomposition: " << t_decomposition.stop() << " s";
        OpmLog::info(out.str());
    }

    std::vector<int> URowsPerColor(numColors);
    rowSize = block_size * rMat->Nb;
    LRowSize = block_size * LMat->Nb;
    URowSize = block_size * UMat->Nb;
    LNumColors = numColors;
    UNumColors = numColors;
    for (int c = 0; c < numColors; c++) {
        URowsPerColor[numColors - c - 1] = rowsPerColor[c];
    }
    int err;
    err = rMat->toRDF(numColors, rowsPerColor.data(), /*isUMatrix:*/ false,
                      colIndicesInColor, maxNNZsPerRow, nnzValsSizes.data(),
                      nnzValues, colIndices.data(), NROffsets.data(), colorSizes.data(), &valSize);
    if (err != 0) {
        return false;
    }
    err = LMat->toRDF(LNumColors, rowsPerColor.data(), /*isUMatrix:*/ false,
                      LColIndicesInColor, maxNNZsPerRow, LnnzValsSizes.data(),
                      LnnzValues, LColIndices.data(), LNROffsets.data(), LColorSizes.data(), &LValSize);
    if (err != 0) {
        return false;
    }
    err = UMat->toRDF(UNumColors, URowsPerColor.data(), /*isUMatrix:*/ true,
                      UColIndicesInColor, maxNNZsPerRow, UnnzValsSizes.data(),
                      UnnzValues, UColIndices.data(), UNROffsets.data(), UColorSizes.data(), &UValSize);
    if (err != 0) {
        return false;
    }
    blockedDiagtoRDF(invDiagVals, rMat->Nb, numColors, URowsPerColor, blockDiag.data());
    // resultPointers are set in the init method
    resultSizes[0] = rowSize;
    resultSizes[1] = colorSizes[3]; // zeropadded valSize;
    resultSizes[2] = numColors;
    resultSizes[3] = colorSizes[2]; //totalCols
    resultSizes[4] = colorSizes[5]; //NRFlagSize
    resultSizes[5] = colorSizes[6]; //diagValsSize
    resultSizes[6] = LRowSize;
    resultSizes[7] = LColorSizes[3]; // zeropadded LValSize;
    resultSizes[8] = LNumColors;
    resultSizes[9] = LColorSizes[2]; //LTotalCols
    resultSizes[10] = LColorSizes[5]; //LNRFlagSize
    resultSizes[11] = LColorSizes[6]; //LDiagValsSize
    resultSizes[12] = URowSize;
    resultSizes[13] = UColorSizes[3]; // zeropadded UValSize;
    resultSizes[14] = UNumColors;
    resultSizes[15] = UColorSizes[2]; //UTotalCols
    resultSizes[16] = UColorSizes[5]; //UNRFlagSize
    resultSizes[17] = UColorSizes[6]; //UDiagValsSize
    return true;
} // end create_preconditioner()


#define INSTANTIATE_BDA_FUNCTIONS(n)                                    \
template FPGABILU0<n>::FPGABILU0(ILUReorder, int, int, int, int, int);  \
template FPGABILU0<n>::~FPGABILU0();                                    \
template bool FPGABILU0<n>::init(BlockedMatrix*);                       \
template bool FPGABILU0<n>::create_preconditioner(BlockedMatrix *);

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
