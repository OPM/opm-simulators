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

namespace Opm
{
namespace Accelerator
{

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

} // namespace Accelerator
} // namespace Opm
