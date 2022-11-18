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

#ifndef OPM_MATRIX_HEADER_INCLUDED
#define OPM_MATRIX_HEADER_INCLUDED

#include <vector>

namespace Opm
{
namespace Accelerator
{

/// This struct resembles a csr matrix, only doubles are supported
/// The data is stored in contiguous memory, such that they can be copied to a device in one transfer.
class Matrix {

public:

    /// Allocate square Matrix and data arrays with given sizes
    /// \param[in] N               number of rows
    /// \param[in] nnzs            number of nonzeros
    Matrix(int N_, int nnzs_)
    : N(N_),
      M(N_),
      nnzs(nnzs_)
    {
        nnzValues.resize(nnzs);
        colIndices.resize(nnzs);
        rowPointers.resize(N+1);
    }

    /// Allocate rectangular Matrix and data arrays with given sizes
    /// \param[in] N               number of rows
    /// \param[in] M               number of columns
    /// \param[in] nnzs            number of nonzeros
    Matrix(int N_, int M_, int nnzs_)
    : Matrix(N_, nnzs_)
    {
        M = M_;
    }

    std::vector<double> nnzValues;
    std::vector<int> colIndices;
    std::vector<int> rowPointers;
    int N, M;
    int nnzs;
};

} // namespace Accelerator
} // namespace Opm

#endif // OPM_MATRIX_HEADER_INCLUDED
