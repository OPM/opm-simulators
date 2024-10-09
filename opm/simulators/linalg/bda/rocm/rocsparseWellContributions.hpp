/*
  Copyright 2023 Equinor ASA

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

#ifndef WELLCONTRIBUTIONS_ROCSPARSE_HEADER_INCLUDED
#define WELLCONTRIBUTIONS_ROCSPARSE_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <hip/hip_runtime_api.h>

#include <vector>

namespace Opm {

template<class Scalar>
class WellContributionsRocsparse : public WellContributions<Scalar>
{
private:
    hipStream_t stream;

public:
    void apply_stdwells(Scalar* d_x, Scalar* d_y);
    void apply_mswells(Scalar* d_x, Scalar* d_y);
    void apply(Scalar* d_x, Scalar* d_y);
    void setStream(hipStream_t stream);

protected:
    /// Allocate memory for the StandardWells
    void APIalloc() override;

    using MatrixType = typename WellContributions<Scalar>::MatrixType;

    void APIaddMatrix(MatrixType type, int* colIndices,
                      Scalar* values, unsigned int val_size) override;

    Scalar *d_Cnnzs_hip, *d_Dnnzs_hip, *d_Bnnzs_hip;
    unsigned *d_Ccols_hip, *d_Bcols_hip;
    unsigned *d_val_pointers_hip;

    std::vector<Scalar> h_x;
    std::vector<Scalar> h_y;
};

} //namespace Opm

#endif
