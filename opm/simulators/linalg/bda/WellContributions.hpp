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
    /// A StandardWell uses C, D and B and performs y -= (C^T * (D^-1 * (B*x)))
    /// B and C are vectors, disguised as matrices and contain blocks of StandardWell::numEq by StandardWell::numStaticWellEq
    /// D is a block, disguised as matrix, the square block has size StandardWell::numStaticWellEq. D is actually stored as D^-1
    /// B*x and D*B*x are a vector with numStaticWellEq doubles
    /// C*D*B*x is a blocked matrix with a symmetric sparsity pattern, contains square blocks with size numEq. For every columnindex i, j in StandardWell::duneB_, there is a block on (i, j) in C*D*B*x.
    ///
    /// This class is used in 3 phases:
    /// - get total size of all wellcontributions that must be stored here
    /// - allocate memory
    /// - copy data of wellcontributions
    class WellContributions
    {

    private:
        unsigned int num_blocks = 0;             // total number of blocks in all wells
        unsigned int dim;                        // number of columns of blocks in B and C, equal to StandardWell::numEq
        unsigned int dim_wells;                  // number of rows of blocks in B and C, equal to StandardWell::numStaticWellEq
        unsigned int num_wells = 0;              // number of wellcontributions in this object
        unsigned int num_blocks_so_far = 0;      // keep track of where next data is written
        unsigned int num_wells_so_far = 0;       // keep track of where next data is written
        unsigned int *val_pointers = nullptr;    // val_pointers[wellID] == index of first block for this well in Ccols and Bcols
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
        /// performs y -= (C^T * (D^-1 * (B*x))) for StandardWell
        /// \param[in] x          vector x
        /// \param[inout] y       vector y
        void apply(double *x, double *y);

        /// Allocate memory for the wellcontributions
        void alloc();

        /// Indicate how large the blocks of the wellcontributions (C and B) are
        /// \param[in] dim         number of columns
        /// \param[in] dim_wells   number of rows
        void setBlockSize(unsigned int dim, unsigned int dim_wells);

        /// Indicate how large the next wellcontribution is, this function cannot be called after alloc() is called
        /// \param[in] numBlocks   number of blocks in C and B of next wellcontribution
        void addNumBlocks(unsigned int numBlocks);

        /// Store a matrix in this object, in blocked csr format, can only be called after alloc() is called
        /// \param[in] idx         indicate if C, D or B is sent
        /// \param[in] colIndices  columnindices of blocks in C or B, ignored for D
        /// \param[in] values      array of nonzeroes
        /// \param[in] val_size    number of blocks in C or B, ignored for D
        void addMatrix(int idx, int *colIndices, double *values, unsigned int val_size);

        /// Return the number of wells added to this object
        /// \return the number of wells added to this object
        unsigned int getNumWells(){
            return num_wells;
        }
    };

} //namespace Opm

#endif
