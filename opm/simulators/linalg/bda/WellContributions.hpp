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

#ifndef WELLCONTRIBUTIONS_HEADER_INCLUDED
#define WELLCONTRIBUTIONS_HEADER_INCLUDED

#include <config.h>

#if ! HAVE_CUDA
  #error "This file should only be included if CUDA is found"
#endif

#include <vector>

#include <cuda_runtime.h>

namespace Opm
{

    /// This class serves to eliminate the need to include the WellContributions into the matrix (with --matrix-add-well-contributions=true) for the cusparseSolver
    /// If the --matrix-add-well-contributions commandline parameter is true, this class should not be used
    class WellContributions
    {

    private:
        unsigned int num_blocks = 0;    // total number of blocks in all wells
        unsigned int dim;
        unsigned int dim_wells;
        unsigned int num_wells = 0;
        unsigned int num_blocks_so_far = 0;
        unsigned int num_wells_so_far = 0;
        unsigned int *val_pointers = nullptr;     // val_pointers[wellID] == index of first block for this well in Ccols and Bcols
        bool allocated = false;

        double *d_Cnnzs = nullptr;
        double *d_Dnnzs = nullptr;
        double *d_Bnnzs = nullptr;
        int *d_Ccols = nullptr;
        int *d_Bcols = nullptr;
        double *d_z1 = nullptr;
        double *d_z2 = nullptr;
        unsigned int *d_val_pointers = nullptr;
        cudaStream_t stream;
    public:
        /// Set a cudaStream to be used
        /// \param[in] stream           the cudaStream that is used to launch the kernel in
        void setCudaStream(cudaStream_t stream);

        /// Create a new WellContributions, implementation is empty
        WellContributions(){};

        /// Destroy a WellContributions, and free memory
        ~WellContributions();

        /// Apply all wellcontributions in this object
        void apply(double *x, double *y);

        /// Allocate memory for the wellcontributions
        void alloc();

        /// Indicate how large the next wellcontributions are, this function cannot be called after alloc_all() is called
        void addSizes(unsigned int nnz, unsigned int numEq, unsigned int numWellEq);

        /// Store a matrix in this object, in blocked csr format
        void addMatrix(int idx, int *colIndices, double *values, unsigned int val_size);

        /// Return the number of wells added to this object
        unsigned int getNumWells(){
            return num_wells;
        }

        /// WellContributions can be applied on CPU or GPU
        /// This function sets the static variable, so each WellContributions is applied on the correct hardware
        static void setMode(bool use_gpu);

    };

} //namespace Opm

#endif
