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

#ifndef MULTISEGMENTWELLCONTRIBUTION_HEADER_INCLUDED
#define MULTISEGMENTWELLCONTRIBUTION_HEADER_INCLUDED

#include <config.h>

#include <vector>

#include <cuda_runtime.h>

namespace Opm
{

    /// This class serves to duplicate the functionality of the MultisegmentWell
    /// A MultisegmentWell uses C, D and B and performs y -= (C^T * (D^-1 * (B*x)))
    /// B and C are matrices, with M rows and N columns, where N is the size of the matrix. They contain blocks of MultisegmentWell::numEq by MultisegmentWell::numWellEq.
    /// D is a MxM matrix, the square blocks have size MultisegmentWell::numWellEq. 
    /// B*x and D*B*x are a vector with M*numWellEq doubles
    /// C*D*B*x is a vector with N*numEq doubles.
    
    class MultisegmentWellContribution
    {

    private:
        unsigned int dim;                        // number of columns of blocks in B and C, equal to MultisegmentWell::numEq
        unsigned int dim_wells;                  // number of rows of blocks in B and C, equal to MultisegmentWell::numWellEq
        unsigned int N;                          // number of rows in vectors x and y, N == dim*Nb
        unsigned int Nb;                         // number of blockrows in x and y
        unsigned int M;                          // number of rows, M == dim_wells*Mb
        unsigned int Mb;                         // number of blockrows in C, D and B

        cudaStream_t stream;
        double *h_x = nullptr, *h_y = nullptr;  // CUDA pinned memory for GPU memcpy

        // C and B are stored in BCRS format, D is stored in CSC format (Dune::UMFPack).
        unsigned int DnumBlocks;
        unsigned int BnumBlocks;
        std::vector<double> Cvals;
        std::vector<double> Dvals;
        std::vector<double> Bvals;
        std::vector<unsigned int> Ccols;
        std::vector<int> Dcols;              // Columnpointers, contains M+1 entries
        std::vector<unsigned int> Bcols;
        std::vector<unsigned int> Crows;
        std::vector<int> Drows;              // Rowindicies, contains DnumBlocks*dim*dim_wells entries
        std::vector<unsigned int> Brows;
        std::vector<double> z1;          // z1 = B * x
        std::vector<double> z2;          // z2 = D^-1 * B * x
        void *UMFPACK_Symbolic, *UMFPACK_Numeric;


    public:

        /// Set a cudaStream to be used
        /// \param[in] stream           the cudaStream that is used
        void setCudaStream(cudaStream_t stream);

        /// Create a new MultisegmentWellContribution
        /// Matrices C and B are passed in Blocked CSR, matrix D in CSC
        /// \param[in] dim              size of blocks in vectors x and y, equal to MultisegmentWell::numEq
        /// \param[in] dim_wells        size of blocks of C, B and D, equal to MultisegmentWell::numWellEq
        /// \param[in] Nb               number of blocks in vectors x and y
        /// \param[in] Mb               number of blockrows in C, B and D
        /// \param[in] BnumBlocks       number of blocks in C and B
        /// \param[in] Bvalues          nonzero values of matrix B
        /// \param[in] BcolIndices      columnindices of blocks of matrix B
        /// \param[in] BrowPointers     rowpointers of matrix B
        /// \param[in] DnumBlocks       number of blocks in D
        /// \param[in] Dvalues          nonzero values of matrix D
        /// \param[in] DcolPointers     columnpointers of matrix D
        /// \param[in] DrowIndices      rowindices of matrix D
        /// \param[in] Cvalues          nonzero values of matrix C
        MultisegmentWellContribution(unsigned int dim, unsigned int dim_wells,
            unsigned int Nb, unsigned int Mb,
            unsigned int BnumBlocks, double *Bvalues, unsigned int *BcolIndices, unsigned int *BrowPointers,
            unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
            double *Cvalues);

        /// Destroy a MultisegmentWellContribution, and free memory
        ~MultisegmentWellContribution();

        /// Apply the MultisegmentWellContribution on GPU
        /// performs y -= (C^T * (D^-1 * (B*x))) for MultisegmentWell
        /// \param[in] d_x          vector x, must be on GPU
        /// \param[inout] d_y       vector y, must be on GPU
        void apply(double *d_x, double *d_y);

    };

} //namespace Opm

#endif
