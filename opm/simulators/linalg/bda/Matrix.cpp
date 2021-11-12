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

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/Matrix.hpp>
#include <opm/simulators/linalg/bda/FPGAUtils.hpp>

namespace Opm
{
namespace Accelerator
{

template <unsigned int block_size>
void OpenclMatrix<block_size>::upload(cl::CommandQueue *queue, double *vals, int *cols, int *rows) {
    std::vector<cl::Event> events(3);

    cl_int err = queue->enqueueWriteBuffer(nnzValues, CL_FALSE, 0, sizeof(double) * block_size * block_size * nnzbs, vals, nullptr, &events[0]);
    err |= queue->enqueueWriteBuffer(colIndices, CL_FALSE, 0, sizeof(int) * nnzbs, cols, nullptr, &events[1]);
    err |= queue->enqueueWriteBuffer(rowPointers, CL_FALSE, 0, sizeof(int) * (Nb + 1), rows, nullptr, &events[2]);

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "OpenclMatrix OpenCL enqueueWriteBuffer error");
    }
}

template <unsigned int block_size>
void OpenclMatrix<block_size>::upload(cl::CommandQueue *queue, Matrix *matrix) {
    upload(queue, matrix->nnzValues.data(), matrix->colIndices.data(), matrix->rowPointers.data());
}

template <unsigned int block_size>
void OpenclMatrix<block_size>::upload(cl::CommandQueue *queue, BlockedMatrix *matrix) {
    upload(queue, matrix->nnzValues, matrix->colIndices, matrix->rowPointers);
}

/*Sort a row of matrix elements from a CSR-format.*/
void sortRow(int *colIndices, double *data, int left, int right) {
    int l = left;
    int r = right;
    int middle = colIndices[(l + r) >> 1];
    do {
        while (colIndices[l] < middle)
            l++;
        while (colIndices[r] > middle)
            r--;
        if (l <= r) {
            int lColIndex = colIndices[l];
            colIndices[l] = colIndices[r];
            colIndices[r] = lColIndex;
            double lDatum = data[l];
            data[l] = data[r];
            data[r] = lDatum;

            l++;
            r--;
        }
    } while (l < r);
    if (left < r)
        sortRow(colIndices, data, left, r);
    if (right > l)
        sortRow(colIndices, data, l, right);

}

#if HAVE_FPGA
/*
 * Write all data used by the VHDL testbenches to raw data arrays. The arrays are as follows:
 * - The "colorSizes" array, which first contains the number of rows, columns, non-zero values
 *   and colors, and the size, in elements, of the NROffsets array, followed by:
 *   the number of rows (rounded to the nearest 32), the number of rows (not rounded),
 *   the number of columns (not rounded) and the number of non-zero values
 *   (rounded to the nearest 32) for every partition.
 *   This array is zero padded up to the nearest 64-byte cacheline.
 * - The "colIndicesInColor" array, which contains for every partition, from which elements
 *   in the global X vector the elements of that X vector partition came.
 *   For example, if a matrix partition only has non-zero values in columns 1, 3 and 6, then
 *   that X vector partition will only have three elements, and the color_col_indices array
 *   will contain 1, 3 and 6 for that partition.
 *   This array is zero padded up to the nearest 64-byte cacheline for every partition.
 * - The "nnzValues" array contains all non-zero values of each partition of the matrix.
 *   This array is zero-padded so that each color has a multiple of 32 elements (to have the
 *   same number of elements per partition as the column indices array).
 * - The "colIndices" array contains all column indices of each partition of the matrix.
 *   These column indices are the local indices for that partition, so to be used, first a
 *   local X vector partition needs to be loaded into some local memory (this is done using
 *   data from the _color_col_indices array), before these column indices can be used as
 *   addresses to that local memory to read the desired X vector values.
 *   This array is zero-padded so that data for every partition fills up a number of complete
 *   cachelines (this means every color has a multiple of 32 elements).
 * - "NROffsets" is the name of the array that contains the new row offsets for
 *   all elements of every partition of the matrix. New row offsets are 8-bit values which
 *   are 0 if that element is not the first element in a row, or which, if that element is
 *   the first element of a row) is equal to the amount of empty rows between that new row
 *   and the row before it plus 1. This array is zero-padded so that data for every partition
 *   fills up a number of complete cachelines (this means every color has a multiple of 64 elements).
 */
int Matrix::toRDF(int numColors, std::vector<int>& nodesPerColor,
                  std::vector<std::vector<int> >& colIndicesInColor, int nnzsThisRowLimit,
                  std::vector<std::vector<double> >& ubNnzValues, short int *ubColIndices, int *nnzValsSizes, unsigned char *NROffsets, int *colorSizes)
{
    auto mat = this;

    int doneRows = 0;
    int totalRowNum = 0;  // total number of non-empty rows
    int nnzsPerColor = 0; // total number of nnzs in current color, padded to multiple of 32 for each color
    int maxNNZsPerColor = 0; // max of nnzsPerColor

    int totalValSize = 0; // sum of nnzsPerColor, padded

    std::vector<int> nnzRowsPerColor(numColors);

    // find number of nnzs per color and number of non-empty rows
    for (int c = 0; c < numColors; c++) {
        int numRows = 0;
        nnzRowsPerColor[c] = 0;
        int firstNnzOfColor = mat->rowPointers[doneRows];
        int lastNnzOfColor = mat->rowPointers[doneRows + nodesPerColor[c]];
        nnzsPerColor = roundUpTo(lastNnzOfColor - firstNnzOfColor, 32); // round up to nearest 16 for short ints of column indices
        totalValSize += nnzsPerColor;
        maxNNZsPerColor = std::max(nnzsPerColor, maxNNZsPerColor);
        int row = doneRows;
        for (; row < doneRows + nodesPerColor[c]; row++) {
            if ( mat->rowPointers[row] != mat->rowPointers[row + 1]) {
                numRows++;
                nnzRowsPerColor[c] = nnzRowsPerColor[c] + 1;
            }
        }

        doneRows = row;
        totalRowNum += numRows;
    }

    int conseqZeroRows = 0;       // number of consecutive empty rows
    int maxConseqZeroRows = 0;
    int numEmptyRows = 0;         // total number of empty rows
    std::vector<int> rowOffsets(totalRowNum);
    std::vector<int> nnzRowPointers(totalRowNum + 1, 0); // rowPointers, but only for non empty rows
    std::vector<int> colorValPointers(numColors + 1);    // points to first nnz of first row of each color
    std::vector<int> colorValZeroPointers(numColors);    // points to first padded zero for each color

    int nonEmptyRowIdx = 0;   // read all rows, but only keep non empty rows, this idx keeps track of how many non empty rows where seen
    doneRows = 0;

    int totalPaddingSize = 0;  // number of padded zeros from previous colors
    int NROffsetSize = 0;  // number of NROffsets entries, padded to multiple of 64 for each color
    int maxRows = 0;
    int maxNNZsPerRow = 0;

    // determine the row offset of each row (amount of zero rows between it and the previous non-zero row)
    // this is later converted to rowOffset for each nnz
    for (int c = 0; c < numColors; c++) {
        conseqZeroRows = 0;
        for (int row = doneRows; row < doneRows + nodesPerColor[c]; row++) {
            int nnzsThisRow = mat->rowPointers[row + 1] - mat->rowPointers[row];
            if (nnzsThisRow == 0) {
                conseqZeroRows++;
                numEmptyRows++;
            } else {
                maxNNZsPerRow = std::max(nnzsThisRow, maxNNZsPerRow);
                nnzRowPointers[nonEmptyRowIdx + 1] = mat->rowPointers[row + 1];
                rowOffsets[nonEmptyRowIdx] = conseqZeroRows;
                maxConseqZeroRows = std::max(conseqZeroRows, maxConseqZeroRows);
                conseqZeroRows = 0;
                nonEmptyRowIdx++;
            }
        }
        // calculate sizes that include zeropadding
        colorValZeroPointers[c] = nnzRowPointers[nonEmptyRowIdx] + totalPaddingSize;
        colorValPointers[c + 1] = roundUpTo(colorValZeroPointers[c], 32);
        totalPaddingSize += colorValPointers[c + 1] - colorValZeroPointers[c];
        NROffsetSize += roundUpTo(colorValPointers[c + 1] - colorValPointers[c], 64);

        doneRows += nodesPerColor[c];
        maxRows = std::max(nodesPerColor[c], maxRows);
    }

    if (maxNNZsPerRow > nnzsThisRowLimit) {
        std::ostringstream errorstring;
        errorstring << "ERROR: Current reordering exceeds maximum number of non-zero values per row limit: " << maxNNZsPerRow << " > " << nnzsThisRowLimit;
        OPM_THROW(std::logic_error, errorstring.str());
    }

    // create and fill RDF arrays
    colorSizes[3] = colorValPointers[numColors];  // total number of nnzs the FPGA has to process, including zeropadding
    colorSizes[5] = NROffsetSize;

    for (int c = 0; c < numColors; c++) {
        colorSizes[c * 4 + 8] = nnzRowsPerColor[c];
        colorSizes[c * 4 + 11] = colorValPointers[c + 1] - colorValPointers[c];
    }

    int rowIndex = 0;  // keep track of where to read/write
    int valIndex = 0;
    int NRIndex = 0;
    int halfwayPoint = colorValPointers[numColors] / 2;
    nnzValsSizes[0] = colorValPointers[numColors];

    colorSizes[7] = halfwayPoint;

    for (int c = 0; c < numColors; c++) {
        int nnzsThisRow;
        // make sure 32 values are written in batches (pad with zeros if needed)
        for (int v = colorValPointers[c]; v < colorValPointers[c + 1]; v += 32) {
            for (int vb = 0; vb < 32; vb++) {

                // if there are enough values for the whole cacheline
                if (v + vb < colorValZeroPointers[c]) {
                    ubNnzValues[0][v + vb] = mat->nnzValues[valIndex];
                    ubColIndices[v + vb] = static_cast<short int>(colIndicesInColor[c][mat->colIndices[valIndex]]);

                    // if this val is the first of a row
                    if (nnzRowPointers[rowIndex] == valIndex) {

                        if (rowOffsets[rowIndex] + 1 >= 255) {
                            std::ostringstream errorstring;
                            errorstring << "ERROR: row offset size exceeded in row " << rowIndex << " with an offset of " << rowOffsets[rowIndex] + 1;
                            OPM_THROW(std::logic_error, errorstring.str());
                        }

                        NROffsets[NRIndex] =  static_cast<unsigned char>(rowOffsets[rowIndex] + 1);

                        // skip all empty rows
                        while (rowIndex < mat->N && nnzRowPointers[rowIndex] == valIndex) {
                            rowIndex++;
                            nnzsThisRow = 0;
                        }
                        nnzsThisRow++;
                    }
                    else
                    {
                        NROffsets[NRIndex] = (unsigned char) 0;
                        nnzsThisRow++;
                    }
                    valIndex++;
                }
                else // zeropadding is needed
                {
                    ubNnzValues[0][v + vb] = 0.0;
                    ubColIndices[v + vb] = static_cast<short int>(colIndicesInColor[c][mat->colIndices[valIndex - 1]]);
                    NROffsets[NRIndex] = 0;
                }
                NRIndex++;

            }
        }

        // zeropad the NROffsets file
        while (NRIndex % 64 != 0) {
            NROffsets[NRIndex] = 0;
            NRIndex++;
        }
    }

    return 0;
}
#endif

#define INSTANTIATE_BDA_FUNCTIONS(n)  \
template class OpenclMatrix<n>;


INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
