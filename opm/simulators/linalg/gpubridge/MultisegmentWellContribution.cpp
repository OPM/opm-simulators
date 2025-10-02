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

#include <config.h> // CMake
#include <opm/simulators/linalg/gpubridge/MultisegmentWellContribution.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif // HAVE_UMFPACK

namespace Opm {

template<class Scalar>
MultisegmentWellContribution<Scalar>::
MultisegmentWellContribution(unsigned int dim_, unsigned int dim_wells_,
                             unsigned int Mb_,
                             std::vector<Scalar>& Bvalues,
                             std::vector<unsigned int>& BcolIndices,
                             std::vector<unsigned int>& BrowPointers,
                             unsigned int DnumBlocks_,
                             Scalar* Dvalues,
                             UMFPackIndex* DcolPointers,
                             UMFPackIndex* DrowIndices,
                             std::vector<Scalar>& Cvalues)
    : dim(dim_)                // size of blockvectors in vectors x and y, equal to MultisegmentWell::numEq
    , dim_wells(dim_wells_)    // size of blocks in C, B and D, equal to MultisegmentWell::numWellEq
    , M(Mb_ * dim_wells)       // number of rows, M == dim_wells*Mb
    , Mb(Mb_)                  // number of blockrows in C, D and B
    , DnumBlocks(DnumBlocks_)  // number of blocks in D
    // copy data for matrix D into vectors to prevent it going out of scope
    , Dvals(Dvalues, Dvalues + DnumBlocks * dim_wells * dim_wells)
    , Dcols(DcolPointers, DcolPointers + M + 1)
    , Drows(DrowIndices, DrowIndices + DnumBlocks * dim_wells * dim_wells)
{
    Cvals = std::move(Cvalues);
    Bvals = std::move(Bvalues);
    Bcols = std::move(BcolIndices);
    Brows = std::move(BrowPointers);

    z1.resize(Mb * dim_wells);
    z2.resize(Mb * dim_wells);

    if constexpr (std::is_same_v<Scalar,float>) {
        OPM_THROW(std::runtime_error, "Cannot use multisegment wells with float");
    } else {
        umfpack_di_symbolic(M, M, Dcols.data(), Drows.data(), Dvals.data(), &UMFPACK_Symbolic, nullptr, nullptr);
        umfpack_di_numeric(Dcols.data(), Drows.data(), Dvals.data(), UMFPACK_Symbolic, &UMFPACK_Numeric, nullptr, nullptr);
    }
}

template<class Scalar>
MultisegmentWellContribution<Scalar>::~MultisegmentWellContribution()
{
    if constexpr (std::is_same_v<Scalar,double>) {
        umfpack_di_free_symbolic(&UMFPACK_Symbolic);
        umfpack_di_free_numeric(&UMFPACK_Numeric);
    }
}

// Apply the MultisegmentWellContribution, similar to MultisegmentWell::apply()
// h_x and h_y reside on host
// y -= (C^T * (D^-1 * (B * x)))
template<class Scalar>
void MultisegmentWellContribution<Scalar>::apply(Scalar* h_x, Scalar* h_y)
{
    OPM_TIMEBLOCK(apply);
    // reset z1 and z2
    std::fill(z1.begin(), z1.end(), 0.0);
    std::fill(z2.begin(), z2.end(), 0.0);

    // z1 = B * x
    for (unsigned int row = 0; row < Mb; ++row) {
        // for every block in the row
        for (unsigned int blockID = Brows[row]; blockID < Brows[row + 1]; ++blockID) {
            unsigned int colIdx = Bcols[blockID];
            for (unsigned int j = 0; j < dim_wells; ++j) {
                Scalar temp = 0.0;
                for (unsigned int k = 0; k < dim; ++k) {
                    temp += Bvals[blockID * dim * dim_wells + j * dim + k] * h_x[colIdx * dim + k];
                }
                z1[row * dim_wells + j] += temp;
            }
        }
    }

    // z2 = D^-1 * (B * x)
    // umfpack
    if constexpr (std::is_same_v<Scalar,float>) {
        OPM_THROW(std::runtime_error, "Cannot use multisegment wells with float");
    } else {
        umfpack_di_solve(UMFPACK_A, Dcols.data(), Drows.data(), Dvals.data(), z2.data(), z1.data(), UMFPACK_Numeric, nullptr, nullptr);
    }

    // y -= (C^T * z2)
    // y -= (C^T * (D^-1 * (B * x)))
    for (unsigned int row = 0; row < Mb; ++row) {
        // for every block in the row
        for (unsigned int blockID = Brows[row]; blockID < Brows[row + 1]; ++blockID) {
            unsigned int colIdx = Bcols[blockID];
            for (unsigned int j = 0; j < dim; ++j) {
                Scalar temp = 0.0;
                for (unsigned int k = 0; k < dim_wells; ++k) {
                    temp += Cvals[blockID * dim * dim_wells + j + k * dim] * z2[row * dim_wells + k];
                }
                h_y[colIdx * dim + j] -= temp;
            }
        }
    }
}

#if HAVE_CUDA
template<class Scalar>
void MultisegmentWellContribution<Scalar>::setCudaStream(cudaStream_t stream_)
{
    stream = stream_;
}
#endif

template class MultisegmentWellContribution<double>;

#if FLOW_INSTANTIATE_FLOAT
template class MultisegmentWellContribution<float>;
#endif

} //namespace Opm
