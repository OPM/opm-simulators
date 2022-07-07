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

#include <opm/simulators/linalg/bda/Reorder.hpp>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <random>
#include <sstream>
#include <vector>
#include <cstring>

namespace {

    std::mt19937 make_urng()
    {
        std::random_device rd;
        std::array<unsigned int, std::mt19937::state_size> seed_data{};

        std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
        std::seed_seq seq(seed_data.begin(), seed_data.end());

        return std::mt19937{ seq };
    }

}

namespace Opm
{
namespace Accelerator
{


/* Give every node in the matrix (of which only the sparsity pattern in the
 * form of row pointers and column indices arrays are in the input), a color
 * in the colors array. Also return the amount of colors in the return integer.
 * This graph-coloring algorithm is based on the Jones-Plassmann algorithm, proposed in:
 * "A Parallel Graph Coloring Heuristic" by M.T. Jones and P.E. Plassmann in SIAM Journal of Scientific Computing 14 (1993) */

template <unsigned int block_size>
int colorBlockedNodes(int rows, const int *CSRRowPointers, const int *CSRColIndices, const int *CSCColPointers, const int *CSCRowIndices, std::vector<int>& colors, int maxRowsPerColor, int maxColsPerColor)
{
    auto left = static_cast<std::vector<int>::difference_type>(colors.size());
    int c = -1;
    const int max_tries = 100;            // since coloring is random, it is possible that a coloring fails. In that case, try again.

    std::vector<bool> visitedColumns(rows, false);

    auto gen = make_urng();

    std::vector<int> randoms(rows);
    for (unsigned int t = 0; t < max_tries; t++) {
        // (re)initialize data for coloring process
        std::uniform_int_distribution<int> uniform{}; // 0 .. INT_MAX

        std::generate(randoms.begin(), randoms.end(),
            [&uniform, &gen]()
        {
            return uniform(gen);
        });

        std::fill(colors.begin(), colors.end(), -1);

        // actually perform coloring
        for (c = 0; c < MAX_COLORS; c++) {
            unsigned int rowsInColor = 0u;
            unsigned int colsInColor = 0u;
            for (int i = 0; i < rows; i++)
            {
                bool iMax = true; // true iff you have max random

                // ignore nodes colored earlier
                if ((colors[i] != -1))
                    continue;

                int ir = randoms[i];

                // look at all nodex that node i is connected to
                for (int k = CSRRowPointers[i]; k < CSRRowPointers[i + 1]; k++) {
                    // ignore nodes colored earlier (and yourself)
                    int j = CSRColIndices[k];
                    int jc = colors[j];
                    if (((jc != -1) && (jc != c)) || (i == j)) {
                        continue;
                    }
                    // node i is not in the current color if one of its neighbours shares this color,
                    if (jc == c) {
                        iMax = false;
                        break;
                    }
                    // or if one of its uncolored neighbours has a higher random value
                    int jr = randoms[j];
                    if (ir <= jr) {
                        iMax = false;
                        break;
                    }
                }
                // look at all nodes that have a connection to node i
                for (int k = CSCColPointers[i]; k < CSCColPointers[i + 1]; k++) {
                    // ignore nodes colored earlier (and yourself)
                    int j = CSCRowIndices[k];
                    int jc = colors[j];
                    if (((jc != -1) && (jc != c)) || (i == j)) {
                        continue;
                    }
                    // node i is not in the current color if one of its neighbours shares this color,
                    if (jc == c) {
                        iMax = false;
                        break;
                    }
                    // or if one of its uncolored neighbours has a higher random value
                    int jr = randoms[j];
                    if (ir <= jr) {
                        iMax = false;
                        break;
                    }
                }

                // assign color if you have the maximum random number
                if (iMax) {
                    unsigned int additionalColsInRow = 0u;
                    for (int k = CSRRowPointers[i]; k < CSRRowPointers[i + 1]; k++) {
                        int j = CSRColIndices[k];
                        if (!visitedColumns[j]) {
                            visitedColumns[j] = true;
                            additionalColsInRow += block_size;
                        }
                    }
                    if ((colsInColor + additionalColsInRow) > static_cast<unsigned int>(maxColsPerColor)) {
                        break;
                    }
                    colsInColor += additionalColsInRow;
                    colors[i] = c;
                    rowsInColor += block_size;
                    if ((rowsInColor + block_size - 1) >= static_cast<unsigned int>(maxRowsPerColor)) {
                        break;
                    }
                }
            }

            // Check if graph coloring is done.
            left = std::count_if(colors.begin(), colors.end(),
                                 [](const int color) { return color == -1; });
            if (left == 0) {
                return c + 1;
            }
        }
    }

    std::ostringstream oss;
    oss << "Error could not find a graph coloring with " << c << " colors after " << max_tries << " tries.\nNumber of colorless nodes: " << left;
    OPM_THROW(std::logic_error, oss.str());
    return -1;
}


/* Reorder a matrix by a specified input order.
 * Both a to order array, which contains for every node from the old matrix where it will move in the new matrix,
 * and the from order, which contains for every node in the new matrix where it came from in the old matrix.
 * reordermapping_nonzeroes is filled with increasing indices, and reordered using the translated colIndices as keys,
 * this means the resulting reordermapping_nonzeroes array contains the mapping
 */
void reorderBlockedMatrixByPattern(BlockedMatrix *mat, std::vector<int>& reordermapping_nonzeroes, int *toOrder, int *fromOrder, BlockedMatrix *rmat){

    int rIndex = 0;
    std::vector<int> tmp(mat->nnzbs);

    reordermapping_nonzeroes.resize(mat->nnzbs);
    for(int i = 0; i < mat->nnzbs; ++i){
        reordermapping_nonzeroes[i] = i;
    }

    rmat->rowPointers[0] = 0;
    for(int i = 0; i < mat->Nb; i++){
        int thisRow = fromOrder[i];
        // put thisRow from the old matrix into row i of the new matrix
        rmat->rowPointers[i+1] = rmat->rowPointers[i] + mat->rowPointers[thisRow+1] - mat->rowPointers[thisRow];
        for(int k = mat->rowPointers[thisRow]; k < mat->rowPointers[thisRow+1]; k++){
            tmp[rIndex] = reordermapping_nonzeroes[k]; // only get 1 entry per block
            rmat->colIndices[rIndex] = mat->colIndices[k];
            rIndex++;
        }
    }
    // re-assign column indices according to the new positions of the nodes referenced by the column indices
    for(int i = 0; i < mat->nnzbs; i++){
        rmat->colIndices[i] = toOrder[rmat->colIndices[i]];
    }
    // re-sort the column indices of every row.
    for(int i = 0; i < mat->Nb; i++){
        sortRow(rmat->colIndices, tmp.data(), rmat->rowPointers[i], rmat->rowPointers[i+1]-1);
    }
    for(int i = 0; i < mat->nnzbs; i++){
        reordermapping_nonzeroes[i] = tmp[i];
    }
    // std::copy();
}

/* Reorder an array of nonzero blocks into another array, using a mapping */
void reorderNonzeroes(BlockedMatrix *mat, std::vector<int>& reordermapping_nonzeroes, BlockedMatrix *rmat){
    assert(mat->block_size == rmat->block_size);

    const unsigned int bs = mat->block_size;

    for(int i = 0; i < mat->nnzbs; i++){
        int old_idx = reordermapping_nonzeroes[i];
        memcpy(rmat->nnzValues+i*bs*bs, mat->nnzValues+old_idx*bs*bs, sizeof(double)*bs*bs); // copy nnz block
    }
}

/* Find a reorder mapping according to the colors that every node of the matrix has received */

void colorsToReordering(int Nb, std::vector<int>& colors, int numColors, int *toOrder, int *fromOrder, std::vector<int>& rowsPerColor) {
    int reordered = 0;

    // Find reordering patterns
    for (int c = 0; c < numColors; c++) {
        for (int i = 0; i < Nb; i++) {
            if (colors[i] == c) {
                rowsPerColor[c]++;
                toOrder[i] = reordered;
                fromOrder[reordered] = i;
                reordered++;
            }
        }
    }
}

// Reorder a vector according to a reordering pattern

template <unsigned int block_size>
void reorderBlockedVectorByPattern(int Nb, double *vector, int *fromOrder, double *rVector) {
    for (int i = 0; i < Nb; i++) {
        for (unsigned int j = 0; j < block_size; j++) {
            rVector[block_size * i + j] = vector[block_size * fromOrder[i] + j];
        }
    }
}


/* Check is operations on a node in the matrix can be started
 * A node can only be started if all nodes that it depends on during sequential execution have already completed.*/

bool canBeStarted(const int rowIndex, const int *rowPointers, const int *colIndices, const std::vector<bool>& doneRows) {
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
 * The level scheduling of a non-symmetric, blocked matrix requires access to a CSC encoding and a CSR encoding of the sparsity pattern of the input matrix.
 * This function is based on a standard level scheduling algorithm, like the one described in:
 * "Iterative methods for Sparse Linear Systems" by Yousef Saad in section 11.6.3
 */

void findLevelScheduling(int *CSRColIndices, int *CSRRowPointers, int *CSCRowIndices, int *CSCColPointers, int Nb, int *numColors, int *toOrder, int* fromOrder, std::vector<int>& rowsPerColor) {
    int activeRowIndex = 0, colorEnd, nextActiveRowIndex = 0;
    int thisRow;
    std::vector<bool> doneRows(Nb, false);
    std::vector <int> rowsToStart;

    // since emplace_back() is used to fill, the vector must be empty
    assert(rowsPerColor.empty());

    // find starting rows: rows that are independent from all rows that come before them.
    for (thisRow = 0; thisRow < Nb; thisRow++) {
        if (canBeStarted(thisRow, CSCColPointers, CSCRowIndices, doneRows)) {
            fromOrder[nextActiveRowIndex] = thisRow;
            toOrder[thisRow] = nextActiveRowIndex;
            nextActiveRowIndex++;
        }
    }
    // 'do' compute on all active rows
    for (colorEnd = 0; colorEnd < nextActiveRowIndex; colorEnd++) {
        doneRows[fromOrder[colorEnd]] = true;
    }

    rowsPerColor.emplace_back(nextActiveRowIndex - activeRowIndex);

    while (colorEnd < Nb) {
        // Go over all rows active from the last color, and check which of their neighbours can be activated this color
        for (; activeRowIndex < colorEnd; activeRowIndex++) {
            thisRow = fromOrder[activeRowIndex];

            for (int i = CSCColPointers[thisRow]; i < CSCColPointers[thisRow + 1]; i++) {
                int thatRow = CSCRowIndices[i];

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
        rowsToStart.clear();
        colorEnd = nextActiveRowIndex;
        rowsPerColor.emplace_back(nextActiveRowIndex - activeRowIndex);
    }

    *numColors = rowsPerColor.size();
}

/* Perform the complete graph coloring algorithm on a matrix. Return an array with the amount of nodes per color.*/

template <unsigned int block_size>
void findGraphColoring(const int *CSRColIndices, const int *CSRRowPointers, const int *CSCRowIndices, const int *CSCColPointers, int Nb, int maxRowsPerColor, int maxColsPerColor, int *numColors, int *toOrder, int *fromOrder, std::vector<int>& rowsPerColor) {
    std::vector<int> rowColor(Nb);

    *numColors = colorBlockedNodes<block_size>(Nb, CSRRowPointers, CSRColIndices, CSCColPointers, CSCRowIndices, rowColor, maxRowsPerColor, maxColsPerColor);

    rowsPerColor.resize(*numColors);
    colorsToReordering(Nb, rowColor, *numColors, toOrder, fromOrder, rowsPerColor);
}

// based on the scipy package from python, scipy/sparse/sparsetools/csr.h on github
void csrPatternToCsc(int *CSRColIndices, int *CSRRowPointers, int *CSCRowIndices, int *CSCColPointers, int Nb) {

    int nnz = CSRRowPointers[Nb];

    // compute number of nnzs per column
    std::fill(CSCColPointers, CSCColPointers + Nb, 0);

    for (int n = 0; n < nnz; ++n) {
        CSCColPointers[CSRColIndices[n]]++;
    }

    // cumsum the nnz per col to get CSCColPointers
    for (int col = 0, cumsum = 0; col < Nb; ++col) {
        int temp = CSCColPointers[col];
        CSCColPointers[col] = cumsum;
        cumsum += temp;
    }
    CSCColPointers[Nb] = nnz;

    for (int row = 0; row < Nb; ++row) {
        for (int j = CSRRowPointers[row]; j < CSRRowPointers[row + 1]; ++j) {
            int col = CSRColIndices[j];
            int dest = CSCColPointers[col];
            CSCRowIndices[dest] = row;
            CSCColPointers[col]++;
        }
    }

    for (int col = 0, last = 0; col <= Nb; ++col) {
        int temp = CSCColPointers[col];
        CSCColPointers[col] = last;
        last = temp;
    }
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                                                            \
template int colorBlockedNodes<n>(int, const int *, const int *, const int *, const int *, std::vector<int>&, int, int);                        \
template void reorderBlockedVectorByPattern<n>(int, double*, int*, double*);                                                                    \
template void findGraphColoring<n>(const int *, const int *, const int *, const int *, int, int, int, int *, int *, int *, std::vector<int>&);  \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
