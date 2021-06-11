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

#include <cstring>
#include <cmath>

#include <config.h>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/FPGAUtils.hpp>

namespace bda
{

using Opm::OpmLog;

/*Sort a row of matrix elements from a blocked CSR-format.*/

template <unsigned int block_size>
void sortBlockedRow(int *colIndices, double *data, int left, int right) {
    const unsigned int bs = block_size;
    int l = left;
    int r = right;
    int middle = colIndices[(l + r) >> 1];
    double lDatum[bs * bs];
    do {
        while (colIndices[l] < middle)
            l++;
        while (colIndices[r] > middle)
            r--;
        if (l <= r) {
            int lColIndex = colIndices[l];
            colIndices[l] = colIndices[r];
            colIndices[r] = lColIndex;
            memcpy(lDatum, data + l * bs * bs, sizeof(double) * bs * bs);
            memcpy(data + l * bs * bs, data + r * bs * bs, sizeof(double) * bs * bs);
            memcpy(data + r * bs * bs, lDatum, sizeof(double) * bs * bs);

            l++;
            r--;
        }
    } while (l < r);

    if (left < r)
        sortBlockedRow<bs>(colIndices, data, left, r);

    if (right > l)
        sortBlockedRow<bs>(colIndices, data, l, right);
}


// LUMat->nnzValues[ik] = LUMat->nnzValues[ik] - (pivot * LUMat->nnzValues[jk]) in ilu decomposition
// a = a - (b * c)
template <unsigned int block_size>
void blockMultSub(double *a, double *b, double *c)
{
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            double temp = 0.0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += b[block_size * row + k] * c[block_size * k + col];
            }
            a[block_size * row + col] -= temp;
        }
    }
}

/*Perform a 3x3 matrix-matrix multiplicationj on two blocks*/

template <unsigned int block_size>
void blockMult(double *mat1, double *mat2, double *resMat) {
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            double temp = 0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += mat1[block_size * row + k] * mat2[block_size * k + col];
            }
            resMat[block_size * row + col] = temp;
        }
    }
}

#if HAVE_FPGA

/*Subtract two blocks from one another element by element*/
template <unsigned int block_size>
void blockSub(double *mat1, double *mat2, double *resMat) {
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            resMat[row * block_size + col] = mat1[row * block_size + col] - mat2[row * block_size + col];
        }
    }
}

/*Multiply a block with a vector block, and add the result, scaled by a constant, to the result vector*/
template <unsigned int block_size>
void blockVectMult(double *mat, double *vect, double scale, double *resVect, bool resetRes) {
    for (unsigned int row = 0; row < block_size; row++) {
        if (resetRes) {
            resVect[row] = 0.0;
        }
        for (unsigned int col = 0; col < block_size; col++) {
            resVect[row] += scale * mat[row * block_size + col] * vect[col];
        }
    }
}



template <unsigned int block_size>
int BlockedMatrix<block_size>::countUnblockedNnzs() {
    int numNnzsOverThreshold = 0;
    int totalNnzs = rowPointers[Nb];
    for (unsigned int idx = 0; idx < totalNnzs * block_size * block_size; idx++) {
        if (fabs(nnzValues[idx]) > nnzThreshold) {
            numNnzsOverThreshold++;
        }
    }
    return numNnzsOverThreshold;
}

/*
 * Unblock the blocked matrix. Input the blocked matrix and output a CSR matrix without blocks.
 * If unblocking the U matrix, the rows in all blocks need to written to the new matrix in reverse order.
*/
template <unsigned int block_size>
void BlockedMatrix<block_size>::unblock(Matrix *mat, bool isUMatrix) {
    const unsigned int bs = block_size;
    int valIndex = 0, nnzsPerRow;

    mat->rowPointers[0] = 0;
    // go through the blocked matrix row-by row of blocks, and then row-by-row inside the block, and
    // write all non-zero values and corresponding column indices that belong to the same row into the new matrix.
    for (int row = 0; row < Nb; row++) {
        for (unsigned int bRow = 0; bRow < bs; bRow++) {
            nnzsPerRow = 0;
            for (int col = rowPointers[row]; col < rowPointers[row + 1]; col++) {
                for (unsigned int bCol = 0; bCol < bs; bCol++) {
                    int idx = 0;
                    // If the matrix is the U matrix, store the rows inside a block in reverse order.
                    if (isUMatrix) {
                        idx = col * bs * bs + (bs - bRow - 1) * bs + bCol;
                    } else {
                        idx = col * bs * bs + bRow * bs + bCol;
                    }

                    if (fabs(nnzValues[idx]) > nnzThreshold) {
                        mat->nnzValues[valIndex] = nnzValues[idx];
                        mat->colIndices[valIndex] = colIndices[col] * bs + bCol;
                        valIndex++;
                        nnzsPerRow++;
                    }
                }
            }
            // Update the rowpointers of the new matrix
            mat->rowPointers[row * bs + bRow + 1] = mat->rowPointers[row * bs + bRow] + nnzsPerRow;
        }
    }
}



/*Optimized version*/
// ub* prefixes indicate unblocked data
template <unsigned int block_size>
int BlockedMatrix<block_size>::toRDF(int numColors, int *nodesPerColor, bool isUMatrix,
                                     std::vector<std::vector<int> >& colIndicesInColor, int nnzsPerRowLimit, int *nnzValsSizes,
                                     std::vector<std::vector<double> >& ubNnzValues, short int *ubColIndices, unsigned char *NROffsets, int *colorSizes, int *valSize)
{
    int res;
    int numUnblockedNnzs = countUnblockedNnzs();

    // initialize the non-blocked matrix with the obtained size.
    std::unique_ptr<Matrix> ubMat = std::make_unique<Matrix>(Nb * block_size, numUnblockedNnzs);

    unblock(ubMat.get(), isUMatrix);

    std::vector<int> ubNodesPerColor(numColors);
    for (int i = 0; i < numColors; i++) {
        ubNodesPerColor[i] = nodesPerColor[i] * block_size;
    }

    *valSize = ubMat->nnzs;

    res = ubMat->toRDF(numColors, ubNodesPerColor,
                       colIndicesInColor, nnzsPerRowLimit,
                       ubNnzValues, ubColIndices, nnzValsSizes,
                       NROffsets, colorSizes);
    return res;
}


// coloring is already done, numColors and nodesPerColor are set
// [rows|columns]PerColorLimit are already queried from the FPGA
// colIndicesInColor, PIndicesAddr and colorSizes are written here
// There are 3 matrices analysed: the full matrix for spmv, L and U for ILU
// node == row
// color == partition
// colorSizes: contains meta info about a color/partition, like number of rows and number of nnzs
// colIndicesInColor: for each color: mapping of colIdx to colValue, unblocked. Used in Matrix::toRDF().
//                    due to partitioning, lots of columns are removed, this matrix keeps track of the mapping
// PIndicesAddr: contiguously for each color: indices of x in global x vector, unblocked
//               if color 0 has A unique colAccesses, PIndicesAddr[0 - A] are for color 0
//               then PIndicesAddr[A - A+B] are for color 1. Directly copied to FPGA
template <unsigned int block_size>
int BlockedMatrix<block_size>::findPartitionColumns(int numColors, int *nodesPerColor,
        int rowsPerColorLimit, int columnsPerColorLimit,
        std::vector<std::vector<int> >& colIndicesInColor, int *PIndicesAddr, int *colorSizes,
        std::vector<std::vector<int> >& LColIndicesInColor, int *LPIndicesAddr, int *LColorSizes,
        std::vector<std::vector<int> >& UColIndicesInColor, int *UPIndicesAddr, int *UColorSizes)
{
    // Data related to column indices per partition
    int doneRows = 0;
    std::vector<bool> isColAccessed(Nb); // std::vector<bool> might have some different optimized implementation, initialize in a loop
    std::vector<bool> isLColAccessed(Nb);
    int totalCols = 0;    // sum of numColAccesses for each color, blocked
    int LTotalCols = 0, UTotalCols = 0;
    int maxCols = 0;         // max value of numColAccesses for any color
    int maxRowsPerColor = 0; // max value of numRows for any color
    int maxColsPerRow = 0;   // max value of colsPerRow for any color
    // colsInColor holds all (blocked) columnIndices that are accessed by that color without duplicates
    // colsInColor[c][i] contains the ith column that color c accesses
    // initial size allows for each color to access all columns, with space for padding
    std::vector<std::vector<int> > colsInColor(numColors, std::vector<int>(roundUpTo(Nb, 16)));
    std::vector<std::vector<int> > LColsInColor(numColors, std::vector<int>(roundUpTo(Nb, 16)));
    std::vector<std::vector<int> > UColsInColor(numColors, std::vector<int>(roundUpTo(Nb, 16)));

    // find which columns are accessed in each color, as well as how many non-zeroes there are per color.
    for (int c = 0; c < numColors; c++) {
        int numRows = 0;
        // initialize
        for (int row = 0; row < Nb; row++) {
            isColAccessed[row] = false;
            isLColAccessed[row] = false;
        }
        if (c > 0) {
            for (int i = doneRows - nodesPerColor[c - 1]; i < doneRows; i++) {
                isLColAccessed[i] = true;
            }
        }
        int numColAccesses = 0, LNumColAccesses = 0, UNumColAccesses = 0;   // number of unique accesses, blocked
        // for every row in this color
        for (int row = doneRows; row < doneRows + nodesPerColor[c]; row++) {
            int colsPerRow = 0;    // number of blocks for this row
            bool rowIsEmpty = (rowPointers[row] == rowPointers[row + 1]);
            for (int idx = rowPointers[row]; idx < rowPointers[row + 1]; idx++) {
                // for every column in the current row, check if that column was accessed before this color
                int col = colIndices[idx];
                if (isColAccessed[col] == false) {
                    colsInColor[c][numColAccesses] = col;
                    isColAccessed[col] = true;
                    numColAccesses++;
                    if (col > row) {
                        UColsInColor[numColors - c - 1][UNumColAccesses] = col;
                        UNumColAccesses++;
                    }
                }
                if (isLColAccessed[col] == false) {
                    if (col < row) {
                        LColsInColor[c][LNumColAccesses] = col;
                        LNumColAccesses++;
                        isLColAccessed[col] = true;
                    }
                }
                colsPerRow++;
            }
            if (rowIsEmpty != true) {
                numRows++;
            }
            maxColsPerRow = std::max(maxColsPerRow, colsPerRow);
        }

        // add columns from previous color into L partition to simplify data forwarding
        if (c > 0) {
            for (int i = doneRows - nodesPerColor[c - 1]; i < doneRows; i++) {
                LColsInColor[c][LNumColAccesses] = i;
                LNumColAccesses++;
            }
        }

        colorSizes[c * 4 + 10] = numColAccesses * block_size;
        LColorSizes[c * 4 + 10] = LNumColAccesses * block_size;
        UColorSizes[(numColors - c - 1) * 4 + 10] = UNumColAccesses * block_size;

        // store mapping
        for (int col = 0; col < numColAccesses; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                colIndicesInColor[c][colsInColor[c][col]*block_size + i] = col * block_size + i;
            }
        }
        for (int col = 0; col < LNumColAccesses; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                LColIndicesInColor[c][LColsInColor[c][col]*block_size + i] = col * block_size + i;
            }
        }
        for (int col = 0; col < UNumColAccesses; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                UColIndicesInColor[numColors - c - 1][UColsInColor[numColors - c - 1][col]*block_size + i] = col * block_size + i;
            }
        }

        // zeropad the colsInColor number to the nearest multiple of 16, because there are 16 32-bit color_col_index values per cacheline
        while (numColAccesses % 16 != 0) {
            colsInColor[c][numColAccesses] = colsInColor[c][numColAccesses - 1];
            numColAccesses++;
        }
        while (LNumColAccesses % 16 != 0) {
            LColsInColor[c][LNumColAccesses] = LColsInColor[c][LNumColAccesses - 1];
            LNumColAccesses++;
        }
        while (UNumColAccesses % 16 != 0) {
            UColsInColor[numColors - c - 1][UNumColAccesses] = UColsInColor[numColors - c - 1][UNumColAccesses - 1];
            UNumColAccesses++;
        }
        maxCols = std::max(numColAccesses, maxCols);
        totalCols += numColAccesses;
        LTotalCols += LNumColAccesses;
        UTotalCols += UNumColAccesses;
        doneRows = doneRows + nodesPerColor[c];
        maxRowsPerColor = std::max(numRows, maxRowsPerColor);
    }

    if (maxCols * static_cast<int>(block_size) > columnsPerColorLimit) {
        std::ostringstream errorstring;
        errorstring << "ERROR: Current reordering exceeds maximum number of columns per color limit: " << maxCols * block_size << " > " << columnsPerColorLimit;
        OPM_THROW(std::logic_error, errorstring.str());
    }

    doneRows = 0;
    int diagValsSize = 0;
    int maxRows = 0;

    for (int c = 0; c < numColors; c++) {
        // calculate sizes that include zeropadding
        diagValsSize += roundUpTo(nodesPerColor[c] * block_size * 4, 8);
        doneRows += nodesPerColor[c];
        if (nodesPerColor[c] * static_cast<int>(block_size) > maxRows)
            maxRows = nodesPerColor[c];
        colorSizes[c * 4 + 9] = nodesPerColor[c] * block_size;
        LColorSizes[c * 4 + 9] = nodesPerColor[c] * block_size;
        UColorSizes[c * 4 + 9] = nodesPerColor[numColors - c - 1] * block_size;
    }

    if (maxRows * static_cast<int>(block_size) > rowsPerColorLimit) {
        std::ostringstream errorstring;
        errorstring << "ERROR: Current reordering exceeds maximum number of columns per color limit: " << maxRows * block_size << " > " << rowsPerColorLimit;
        OPM_THROW(std::logic_error, errorstring.str());
    }

    // create and fill sizes array as far as already possible
    colorSizes[0] = Nb * block_size;
    LColorSizes[0] = Nb * block_size;
    UColorSizes[0] = Nb * block_size;
    // col_sizes (but the matrix is square)
    colorSizes[1] = Nb * block_size;
    LColorSizes[1] = Nb * block_size;
    UColorSizes[1] = Nb * block_size;
    colorSizes[2] = totalCols * block_size;
    LColorSizes[2] = LTotalCols * block_size;
    UColorSizes[2] = UTotalCols * block_size;
    // missing val_size, written in Matrix::toRDF()
    colorSizes[4] = numColors;
    LColorSizes[4] = numColors;
    UColorSizes[4] = numColors;
    // missing NRFlagsSize, written in Matrix::toRDF()
    colorSizes[6] = diagValsSize;
    LColorSizes[6] = diagValsSize;
    UColorSizes[6] = diagValsSize;

    int paddingIdx = numColors;
    while (paddingIdx % 4 != 0) {
        for (unsigned int i = 0; i < 4; i++) {
            colorSizes[paddingIdx * 4 + 8 + i] = 0;
            LColorSizes[paddingIdx * 4 + 8 + i] = 0;
            UColorSizes[paddingIdx * 4 + 8 + i] = 0;
        }
        paddingIdx++;
    }

    int index = 0, Lindex = 0, Uindex = 0;
    for (int c = 0; c < numColors; c++) {
        // for each unique col access
        for (int col = 0; col < colorSizes[c * 4 + 10] / static_cast<int>(block_size) ; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                PIndicesAddr[index] = colsInColor[c][col] * block_size + i;
                index++;
            }
        }
        // add padding
        while (index % 16 != 0) {
            PIndicesAddr[index] = PIndicesAddr[index - 1];
            index++;
        }
        for (int col = 0; col < LColorSizes[c * 4 + 10] / static_cast<int>(block_size) ; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                LPIndicesAddr[Lindex] = LColsInColor[c][col] * block_size + i;
                Lindex++;
            }
        }
        while (Lindex % 16 != 0) {
            LPIndicesAddr[Lindex] = LPIndicesAddr[Lindex - 1];
            Lindex++;
        }
        for (int col = 0; col < UColorSizes[c * 4 + 10] / static_cast<int>(block_size) ; col++) {
            for (unsigned int i = 0; i < block_size; i++) {
                UPIndicesAddr[Uindex] = UColsInColor[c][col] * block_size + i;
                Uindex++;
            }
        }
        while (Uindex % 16 != 0) {
            UPIndicesAddr[Uindex] = UPIndicesAddr[Uindex - 1];
            Uindex++;
        }
    }
    return 0;
}


void blockedDiagtoRDF(double *blockedDiagVals, int rowSize, int numColors, std::vector<int>& rowsPerColor, double *RDFDiag) {
    const unsigned int block_size = 3;
    int doneRows = rowSize - 1;    // since the rows of U are reversed, the rows of the diag are also reversed
    int RDFIndex = 0;
    for (int c = 0; c < numColors; c++) {
        for (int r = 0; r < rowsPerColor[c]; r++) {

            // the rows in the block are reversed
            for (int i = static_cast<int>(block_size) - 1; i >= 0; i--) {
                for (unsigned int j = 0; j < block_size; j++) {
                    RDFDiag[RDFIndex] = blockedDiagVals[(doneRows - r) * block_size * block_size + i * block_size + j];
                    RDFIndex++;
                }
                // add 4th column, zeropadding
                RDFDiag[RDFIndex] = 0.0;
                RDFIndex++;
            }
        }
        doneRows -= rowsPerColor[c];

        // make sure the color completely fills a cacheline
        // a color with 3 blocks would otherwise leave space
        while (RDFIndex % 8 != 0) {
            RDFDiag[RDFIndex] = 0.0;
            RDFIndex++;
        }
    }
    assert(RDFIndex % 8 == 0);
}

#endif // HAVE_FPGA



#define INSTANTIATE_BDA_FUNCTIONS(n)                                        \
template void sortBlockedRow<n>(int *, double *, int, int);                 \
template void blockMultSub<n>(double *, double *, double *);                \
template void blockMult<n>(double *, double *, double *);                   \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

#if HAVE_FPGA
#define INSTANTIATE_BDA_FPGA_FUNCTIONS(n)                                             \
template void blockSub<n>(double *, double *, double *);                              \
template void blockVectMult<n>(double *, double *, double, double *, bool);           \
template int BlockedMatrix<n>::toRDF(int, int *, bool,                                \
    std::vector<std::vector<int> >& , int, int *,                                     \
    std::vector<std::vector<double> >&, short int *, unsigned char *, int *,  int *); \
template int BlockedMatrix<n>::findPartitionColumns(int, int *,                       \
        int, int,                                                                     \
        std::vector<std::vector<int> >& , int *, int *,                               \
        std::vector<std::vector<int> >& , int *, int *,                               \
        std::vector<std::vector<int> >& , int *, int *);

INSTANTIATE_BDA_FPGA_FUNCTIONS(1);
INSTANTIATE_BDA_FPGA_FUNCTIONS(2);
INSTANTIATE_BDA_FPGA_FUNCTIONS(3);
INSTANTIATE_BDA_FPGA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FPGA_FUNCTIONS
#endif // HAVE_FPGA


} // end namespace bda
