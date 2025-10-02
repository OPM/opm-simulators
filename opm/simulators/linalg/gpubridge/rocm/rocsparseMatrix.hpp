/*
  Copyright 2024 Equinor ASA

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

#ifndef OPM_ROCMMATRIX_HEADER_INCLUDED
#define OPM_ROCMMATRIX_HEADER_INCLUDED

#include <hip/hip_runtime_api.h>

namespace Opm::Accelerator {

template<class Scalar> class Matrix;
template<class Scalar> class BlockedMatrix;

/// This struct resembles a csr matrix
template<class Scalar>
class RocmMatrix {
public:

    RocmMatrix(int Nb_, int Mb_, int nnzbs_, unsigned int block_size_);
    ~RocmMatrix();

    void upload(Scalar *vals,
                int *cols,
                int *rows,
                hipStream_t stream);

    void upload(Matrix<Scalar> *matrix,
                hipStream_t stream);

    void upload(BlockedMatrix<Scalar> *matrix,
                hipStream_t stream);

    Scalar* nnzValues;
    int* colIndices;
    int* rowPointers;
    int Nb, Mb;
    int nnzbs;
    unsigned int block_size;
};

template <typename Scalar>
class RocmVector {
public:

    RocmVector(int N);
    ~RocmVector();

    void upload(Scalar *vals,
                hipStream_t stream);

    void upload(Matrix<Scalar> *matrix,
                hipStream_t stream);

    Scalar* nnzValues;
    int size;
};
} // namespace Opm

#endif
