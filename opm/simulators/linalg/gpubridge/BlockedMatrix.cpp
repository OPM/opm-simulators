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

#include <opm/simulators/linalg/gpubridge/BlockedMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/Matrix.hpp>
#include <opm/simulators/linalg/gpubridge/Matrix.hpp>

namespace Opm::Accelerator {

void sortRow(int *colIndices, int *data, int left, int right)
{
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
            int tmp = data[l];
            data[l] = data[r];
            data[r] = tmp;

            l++;
            r--;
        }
    } while (l < r);

    if (left < r)
        sortRow(colIndices, data, left, r);

    if (right > l)
        sortRow(colIndices, data, l, right);
}

// LUMat->nnzValues[ik] = LUMat->nnzValues[ik] - (pivot * LUMat->nnzValues[jk]) in ilu decomposition
// a = a - (b * c)
template<class Scalar>
void blockMultSub(Scalar* a, Scalar* b, Scalar* c, unsigned int block_size)
{
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            Scalar temp = 0.0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += b[block_size * row + k] * c[block_size * k + col];
            }
            a[block_size * row + col] -= temp;
        }
    }
}

/*Perform a 3x3 matrix-matrix multiplicationj on two blocks*/
template<class Scalar>
void blockMult(Scalar* mat1, Scalar* mat2, Scalar* resMat, unsigned int block_size)
{
    for (unsigned int row = 0; row < block_size; row++) {
        for (unsigned int col = 0; col < block_size; col++) {
            Scalar temp = 0;
            for (unsigned int k = 0; k < block_size; k++) {
                temp += mat1[block_size * row + k] * mat2[block_size * k + col];
            }
            resMat[block_size * row + col] = temp;
        }
    }
}

#define INSTANTIATE_TYPE(T)                               \
    template void blockMultSub(T*, T*, T*, unsigned int); \
    template void blockMult(T*, T*, T*, unsigned int);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::Accelerator
