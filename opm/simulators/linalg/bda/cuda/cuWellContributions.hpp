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

#ifndef WELLCONTRIBUTIONS_CUDA_HEADER_INCLUDED
#define WELLCONTRIBUTIONS_CUDA_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <cuda_runtime.h>


namespace Opm
{

class WellContributionsCuda : public WellContributions
{
public:
    ~WellContributionsCuda() override;

    /// Set a cudaStream to be used
    /// \param[in] stream           the cudaStream that is used to launch the kernel in
    void setCudaStream(cudaStream_t stream);

    /// Apply all Wells in this object
    /// performs y -= (C^T * (D^-1 * (B*x))) for all Wells
    /// \param[in] d_x        vector x, must be on GPU
    /// \param[inout] d_y     vector y, must be on GPU
    void apply(double *d_x, double *d_y);

protected:
    /// Allocate memory for the StandardWells
    void APIalloc() override;

    /// Store a matrix in this object, in blocked csr format, can only be called after alloc() is called
    /// \param[in] type        indicate if C, D or B is sent
    /// \param[in] colIndices  columnindices of blocks in C or B, ignored for D
    /// \param[in] values      array of nonzeroes
    /// \param[in] val_size    number of blocks in C or B, ignored for D
    void APIaddMatrix(MatrixType type, int *colIndices, double *values, unsigned int val_size) override;

    cudaStream_t stream;

    // data for StandardWells, could remain nullptrs if not used
    double *d_Cnnzs = nullptr;
    double *d_Dnnzs = nullptr;
    double *d_Bnnzs = nullptr;
    int *d_Ccols = nullptr;
    int *d_Bcols = nullptr;
    double *d_z1 = nullptr;
    double *d_z2 = nullptr;
    unsigned int *d_val_pointers = nullptr;
    double* h_x = nullptr;
    double* h_y = nullptr;

};

} //namespace Opm

#endif
