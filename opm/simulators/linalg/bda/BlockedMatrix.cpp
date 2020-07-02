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

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

using bda::BlockedMatrix;

namespace bda
{

template <unsigned int block_size>
BlockedMatrix *allocateBlockedMatrix(int Nb, int nnzbs) {
    BlockedMatrix *mat = new BlockedMatrix();

    mat->nnzValues  = new double[nnzbs * block_size * block_size];
    mat->colIndices = new int[nnzbs];
    mat->rowPointers = new int[Nb + 1];
    mat->Nb = Nb;
    mat->nnzbs = nnzbs;

    return mat;
}

void freeBlockedMatrix(BlockedMatrix **mat) {
    if (*mat) {
        delete[] (*mat)->nnzValues;
        delete[] (*mat)->colIndices;
        delete[] (*mat)->rowPointers;
        delete (*mat);
        *mat = NULL;
    }
}

BlockedMatrix *softCopyBlockedMatrix(BlockedMatrix *mat) {
    BlockedMatrix *res = new BlockedMatrix();
    res->nnzValues = mat->nnzValues;
    res->colIndices = mat->colIndices;
    res->rowPointers = mat->rowPointers;
    res->Nb = mat->Nb;
    res->nnzbs = mat->nnzbs;
    return res;
}


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


#define INSTANTIATE_BDA_FUNCTIONS(n)                                 \
template BlockedMatrix *allocateBlockedMatrix<n>(int Nb, int nnzbs); \
template void sortBlockedRow<n>(int *, double *, int, int);          \
template void blockMultSub<n>(double *, double *, double *);         \
template void blockMult<n>(double *, double *, double *);            \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // end namespace bda
