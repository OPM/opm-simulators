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

#ifndef OPM_HIPKERNELS_HPP
#define OPM_HIPKERNELS_HPP

#include <string>
#include <memory>
#include <cstddef>

#include <hip/hip_runtime_api.h>
#include <hip/hip_version.h>

namespace Opm {

template<class Scalar>
class HipKernels
{
private:
    static int  verbosity;
    static bool initialized;
    
    HipKernels();
    
public:
    /// Initialize verbosity level for the HIP kernels
    /// \param[in] verbosity   verbosity level
    static void init(int verbosity);
    
    /// Transform blocked vector to scalar vector using pressure-weights, where every workitem handles one blockrow
    /// \param[in]  fine_y     Input y vector
    /// \param[in]  weights    Weights used to combine cells
    /// \param[out] coarse_y   Output y vector
    /// \param[in]  Nb         Number of blocks in the original matrix
    /// \param[in]  stream     Hip stream to use for the computations
    static void full_to_pressure_restriction(const Scalar* fine_y,
                                             Scalar* weights,
                                             Scalar* coarse_y,
                                             int Nb,
                                             hipStream_t stream);

    /// Add the coarse pressure solution back to the finer, complete solution; every workitem handles one blockrow
    /// \param[in]  coarse_x     Input scalar x vector
    /// \param[out] fine_x       Output blocked x vector 
    /// \param[in]  pressure_idx Pressure index
    /// \param[in]  Nb           Number of blocks in the original matrix
    /// \param[in]  stream       Hip stream to use for the computations
    static void add_coarse_pressure_correction(Scalar* coarse_x,
                                               Scalar* fine_x,
                                               int pressure_idx,
                                               int Nb,
                                               hipStream_t stream);


    /// Function to multiply vector with another vector and a scalar, element-wise and add the result to a third vector (out = alpha * in1 + in2) 
    /// \param[in]  alpha   Input scalar 
    /// \param[in]  in1     First input vector 
    /// \param[in]  in2     Second input vector
    /// \param[out] out     Output vector
    /// \param[in]  N       Size of the vector
    /// \param[in]  stream  Hip stream to use for the computations
    static void vmul(const Scalar alpha,
                     Scalar* in1,
                     Scalar* in2,
                     Scalar* out,
                     int N,
                     hipStream_t stream);
    
    /// Function to prolongate vector during amg cycle, every workitem handles one row
    /// \param[in]  in      Input fine-grained vector
    /// \param[out] out     Output course-graned vector 
    /// \param[in]  cols    Column indexes
    /// \param[in]  N       Size of the vector
    /// \param[in]  stream  Hip stream to use for the computations
    static void prolongate_vector(const Scalar* in,
                                  Scalar* out,
                                  const int* cols,
                                  int N,
                                  hipStream_t stream);

    /// Function to perform res = rhs - mat * x
    /// \param[in]  vals       Matrix values
    /// \param[in]  cols       Column indexes
    /// \param[in]  rows       Row pointers
    /// \param[in]  x          X vector
    /// \param[in]  rhs        Rhs vector
    /// \param[out] out        Output res vector
    /// \param[in]  Nb         Number of non-zero blocks in the original matrix
    /// \param[in]  block_size Block size
    /// \param[in]  stream     Hip stream to use for the computations
    static void residual(Scalar* vals,
                         int* cols,
                         int* rows,
                         Scalar* x,
                         const Scalar* rhs,
                         Scalar* out,
                         int Nb,
                         unsigned int block_size,
                         hipStream_t stream);

    /// Function to perform sparse matrix vector multipliation
    /// \param[in]  vals       Matrix values
    /// \param[in]  cols       Column indexes
    /// \param[in]  rows       Row pointers
    /// \param[in]  x          Input x vector
    /// \param[out] y          Output y vector
    /// \param[in]  Nb         Number of non-zero blocks in the original matrix
    /// \param[in]  block_size Block size
    /// \param[in]  stream     Hip stream to use for the computations
    static void spmv(Scalar* vals,
                     int* cols,
                     int* rows,
                     Scalar* x,
                     Scalar* y,
                     int Nb,
                     unsigned int block_size,
                     hipStream_t stream);
};

} // namespace Opm

#endif
