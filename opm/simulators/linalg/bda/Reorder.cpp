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
#include <vector>
#include <cstring>
#include <algorithm> // for fill()
#include <sys/time.h>

#include <opm/simulators/linalg/bda/Reorder.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

namespace bda
{


/* Give every node in the matrix (of which only the sparsity pattern in the
 * form of row pointers and column indices arrays are in the input), a color
 * in the colors array. Also return the amount of colors in the return integer. */

int colorBlockedNodes(int rows, const int *rowPointers, const int *colIndices, std::vector<int>& colors, int maxRowsPerColor, int maxColsPerColor)
{
    int left, c, t, i, j, k;
    const int max_tries = 100;            // since coloring is random, it is possible that a coloring fails. In that case, try again.
    std::vector<int> randoms;
    randoms.reserve(rows);

    std::vector<bool> visitedColumns;
    std::fill(visitedColumns.begin(), visitedColumns.end(), false);
    int colsInColor;
    int additionalColsInRow;

    for (t = 0; t < max_tries; t++) {
        struct timeval tm;
        gettimeofday(&tm, nullptr);
        srand(tm.tv_sec + tm.tv_usec * 1000000ul);
        for (i = 0; i  < rows; i++) {
            randoms[i] = rand();
            colors[i] = -1;
        }

        for (c = 0; c < MAX_COLORS; c++) {
            int rowsInColor = 0;
            colsInColor = 0;
            for (i = 0; i < rows; i++)
            {
                char f = 1; // true iff you have max random

                // ignore nodes colored earlier
                if ((colors[i] != -1))
                    continue;

                int ir = randoms[i];

                // look at neighbors to check their random number
                for (k = rowPointers[i]; k < rowPointers[i + 1]; k++) {

                    // ignore nodes colored earlier (and yourself)
                    j = colIndices[k];
                    int jc = colors[j];
                    if (((jc != -1) && (jc != c)) || (i == j)) {
                        continue;
                    }
                    // The if statement below makes it both true graph coloring and no longer guaranteed to converge
                    if (jc == c) {
                        f = 0;
                        break;
                    }
                    int jr = randoms[j];
                    if (ir <= jr) f = 0;
                }

                // assign color if you have the maximum random number
                if (f == 1) {
                    additionalColsInRow = 0;
                    for (k = rowPointers[i]; k < rowPointers[i + 1]; k++) {
                        j = colIndices[k];
                        if (!visitedColumns[j]) {
                            visitedColumns[j] = true;
                            additionalColsInRow += 3;
                        }
                    }
                    if ((colsInColor + additionalColsInRow) > maxColsPerColor) {
                        break;
                    }
                    colsInColor += additionalColsInRow;
                    colors[i] = c;
                    rowsInColor += 3;
                    if ((rowsInColor + 2) >= maxRowsPerColor) {
                        break;
                    }

                }

            }
            // Check if graph coloring is done.
            left = 0;
            for (k = 0; k < rows; k++) {
                if (colors[k] == -1) {
                    left++;
                }
            }
            if (left == 0) {
                return c + 1;
            }
        }
    }

    printf("Could not find a graph coloring with %d colors after %d tries.\nAmount of colorless nodes: %d.\n", c, t, left);
    return -1;
}


/* Reorder a matrix by a specified input order.
 * Both a to order array, which contains for every node from the old matrix where it will move in the new matrix,
 * and the from order, which contains for every node in the new matrix where it came from in the old matrix.*/

template <unsigned int block_size>
void blocked_reorder_matrix_by_pattern(BlockedMatrix *mat, int *toOrder, int *fromOrder, BlockedMatrix *rMat) {
    const unsigned int bs = block_size;
    int rIndex = 0;
    int i, k;
    unsigned int j;

    rMat->rowPointers[0] = 0;
    for (i = 0; i < mat->Nb; i++) {
        int thisRow = fromOrder[i];
        // put thisRow from the old matrix into row i of the new matrix
        rMat->rowPointers[i + 1] = rMat->rowPointers[i] + mat->rowPointers[thisRow + 1] - mat->rowPointers[thisRow];
        for (k = mat->rowPointers[thisRow]; k < mat->rowPointers[thisRow + 1]; k++) {
            for (j = 0; j < bs * bs; j++){
                rMat->nnzValues[rIndex * bs * bs + j] = mat->nnzValues[k * bs * bs + j];
            }
            rMat->colIndices[rIndex] = mat->colIndices[k];
            rIndex++;
        }
    }
    // re-assign column indices according to the new positions of the nodes referenced by the column indices
    for (i = 0; i < mat->nnzbs; i++) {
        rMat->colIndices[i] = toOrder[rMat->colIndices[i]];
    }
    // re-sort the column indices of every row.
    for (i = 0; i < mat->Nb; i++) {
        sortBlockedRow<bs>(rMat->colIndices, rMat->nnzValues, rMat->rowPointers[i], rMat->rowPointers[i + 1] - 1);
    }
}

/* Reorder a matrix according to the colors that every node of the matrix has received*/

void colorsToReordering(int Nb, std::vector<int>& colors, int *toOrder, int *fromOrder, int *iters) {
    int reordered = 0;
    int i, c;
    for (i = 0; i < MAX_COLORS; i++) {
        iters[i] = 0;
    }

    // Find reordering patterns
    for (c = 0; c < MAX_COLORS; c++) {
        for (i = 0; i < Nb; i++) {
            if (colors[i] == c) {
                iters[c]++;
                toOrder[i] = reordered;

                fromOrder[reordered] = i;
                reordered++;
            }
        }
    }
}

// Reorder a matrix according to a reordering pattern

template <unsigned int block_size>
void blocked_reorder_vector_by_pattern(int Nb, double *vector, int *fromOrder, double *rVector) {
    for (int i = 0; i < Nb; i++) {
        for (unsigned int j = 0; j < block_size; j++) {
            rVector[block_size * i + j] = vector[block_size * fromOrder[i] + j];
        }
    }
}


/* Check is operations on a node in the matrix can be started
 * A node can only be started if all nodes that it depends on during sequential execution have already completed.*/

bool canBeStarted(int rowIndex, int *rowPointers, int *colIndices, std::vector<bool>& doneRows) {
    bool canStart = !doneRows[rowIndex];
    int i, thisDependency;
    if (canStart) {
        for (i = rowPointers[rowIndex]; i < rowPointers[rowIndex + 1]; i++) {
            thisDependency = colIndices[i];
            // Only dependencies on rows that should execute before the current one are relevant
            if (thisDependency >= rowIndex)
                break;
            // Check if dependency has been resolved
            if (!doneRows[thisDependency]) {
                return false;
            }
        }
    }
    return canStart;
}

/*
 * The level scheduling of a non-symmetric, blocked matrix requires access to a CSC encoding and a CSR encoding of the same matrix.
*/

int *findLevelScheduling(int *CSRColIndices, int *CSRRowPointers, int *CSCColIndices, int *CSCRowPointers, int Nb, int *iters, int *toOrder, int* fromOrder) {
    int activeRowIndex = 0, iterEnd, nextActiveRowIndex = 0;
    int thisRow;
    std::vector<bool> doneRows(Nb, false);
    std::vector<int> rowsPerIter;
    rowsPerIter.reserve(Nb);
    int *resRowsPerIter;

    std::vector <int> rowsToStart;

    // find starting rows: rows that are independent from all rows that come before them.
    for (thisRow = 0; thisRow < Nb; thisRow++) {
        if (canBeStarted(thisRow, CSRRowPointers, CSCColIndices, doneRows)) {
            fromOrder[nextActiveRowIndex] = thisRow;
            toOrder[thisRow] = nextActiveRowIndex;
            nextActiveRowIndex++;
        }
    }
    // 'do' compute on all active rows
    for (iterEnd = 0; iterEnd < nextActiveRowIndex; iterEnd++) {
        doneRows[fromOrder[iterEnd]] = true;
    }

    rowsPerIter.emplace_back(nextActiveRowIndex - activeRowIndex);

    while (iterEnd < Nb) {
        // Go over all rows active from the last iteration, and check which of their neighbours can be activated this iteration
        for (; activeRowIndex < iterEnd; activeRowIndex++) {
            thisRow = fromOrder[activeRowIndex];

            for (int i = CSCRowPointers[thisRow]; i < CSCRowPointers[thisRow + 1]; i++) {
                int thatRow = CSCColIndices[i];

                if (canBeStarted(thatRow, CSRRowPointers, CSRColIndices, doneRows)) {
                    rowsToStart.emplace_back(thatRow);
                }
            }
        }
        // 'do' compute on all active rows
        for (unsigned int i = 0; i < rowsToStart.size(); i++) {
            thisRow = rowsToStart[i];
            if (!doneRows[thisRow]) {
                doneRows[thisRow] = true;
                fromOrder[nextActiveRowIndex] = thisRow;
                toOrder[thisRow] = nextActiveRowIndex;
                nextActiveRowIndex++;
            }
        }
        iterEnd = nextActiveRowIndex;
        rowsPerIter.emplace_back(nextActiveRowIndex - activeRowIndex);
    }
    // Crop the rowsPerIter array to it minimum size.
    int numColors = rowsPerIter.size();
    resRowsPerIter = new int[numColors];
    for (int i = 0; i < numColors; i++) {
        resRowsPerIter[i] = rowsPerIter[i];
    }

    *iters = rowsPerIter.size();

    return resRowsPerIter;
}

/* Perform the complete graph coloring algorithm on a matrix. Return an array with the amount of nodes per color.*/

int* findGraphColoring(int *colIndices, int *rowPointers, int Nb, int maxRowsPerColor, int maxColsPerColor, int *numColors, int *toOrder, int* fromOrder) {
    std::vector<int> rowColor;
    rowColor.reserve(Nb);
    int *rowsPerColor = new int[MAX_COLORS];

    if (colorBlockedNodes(Nb, rowPointers, colIndices, rowColor, maxRowsPerColor, maxColsPerColor) == -1) {
        return nullptr;
    }

    colorsToReordering(Nb, rowColor, toOrder, fromOrder, rowsPerColor);

    *numColors = MAX_COLORS;
    while (rowsPerColor[*numColors - 1] == 0) {
        *numColors = *numColors - 1;
    }

    return rowsPerColor;
}

// based on the scipy package from python, scipy/sparse/sparsetools/csr.h on github
// input : matrix A via Avals, Acols, Arows, Nb
// output: matrix B via Bvals, Bcols, Brows
// arrays for B must be preallocated
template <unsigned int block_size>
void bcsr_to_bcsc(double *Avals, int *Acols, int *Arows, double *Bvals, int *Bcols, int *Brows, int Nb) {

    int nnz = Arows[Nb];

    // compute number of nnzs per column
    std::fill(Brows, Brows + Nb, 0);

    for (int n = 0; n < nnz; ++n) {
        Brows[Acols[n]]++;
    }

    // cumsum the nnz per col to get Brows
    for (int col = 0, cumsum = 0; col < Nb; ++col) {
        int temp = Brows[col];
        Brows[col] = cumsum;
        cumsum += temp;
    }
    Brows[Nb] = nnz;

    for (int row = 0; row < Nb; ++row) {
        for (int j = Arows[row]; j < Arows[row + 1]; ++j) {
            int col = Acols[j];
            int dest = Brows[col];
            Bcols[dest] = row;
            memcpy(Bvals + dest, Avals + dest, sizeof(double) * block_size * block_size);
            Brows[col]++;
        }
    }

    for (int col = 0, last = 0; col <= Nb; ++col) {
        int temp = Brows[col];
        Brows[col] = last;
        last = temp;
    }
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                 \
template void blocked_reorder_matrix_by_pattern<n>(BlockedMatrix *, int *, int *, BlockedMatrix *);  \
template void blocked_reorder_vector_by_pattern<n>(int, double*, int*, double*);                     \
template void bcsr_to_bcsc<n>(double *, int *, int *, double *, int *, int *, int );                 \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} //namespace bda