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
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/MultisegmentWellContribution.hpp>

#ifdef HAVE_OPENCL
#include <opm/simulators/linalg/bda/opencl/openclWellContributions.hpp>
#endif

#ifdef HAVE_CUDA
#include <opm/simulators/linalg/bda/cuda/cuWellContributions.hpp>
#endif

#ifdef HAVE_ROCSPARSE
#include <opm/simulators/linalg/bda/rocm/rocsparseWellContributions.hpp>
#endif

namespace Opm {

template<class Scalar>
WellContributions<Scalar>::~WellContributions() = default;

template<class Scalar>
std::unique_ptr<WellContributions<Scalar>>
WellContributions<Scalar>::create(const std::string& accelerator_mode, bool useWellConn)
{
    if (accelerator_mode.compare("cusparse") == 0) {
#if HAVE_CUDA
        return std::make_unique<WellContributionsCuda<Scalar>>();
#else
        OPM_THROW(std::runtime_error, "Cannot initialize well contributions: CUDA is not enabled");
#endif
    }
    else if (accelerator_mode.compare("opencl") == 0) {
#if HAVE_OPENCL
        return std::make_unique<WellContributionsOCL<Scalar>>();
#else
        OPM_THROW(std::runtime_error, "Cannot initialize well contributions: OpenCL is not enabled");
#endif
    }
    else if (accelerator_mode.compare("rocsparse") == 0) {
        if (!useWellConn) {
#if HAVE_ROCSPARSE
            return std::make_unique<WellContributionsRocsparse<Scalar>>();
#else
        OPM_THROW(std::runtime_error, "Cannot initialize well contributions: rocsparse is not enabled");
#endif
        }
        return std::make_unique<WellContributions>();
    }
    else if(accelerator_mode.compare("amgcl") == 0){
        if (!useWellConn) {
            OPM_THROW(std::logic_error, "Error amgcl requires --matrix-add-well-contributions=true");
        }
        return std::make_unique<WellContributions>();
    }
    else if(accelerator_mode.compare("rocalution") == 0){
        if (!useWellConn) {
            OPM_THROW(std::logic_error, "Error rocalution requires --matrix-add-well-contributions=true");
        }
        return std::make_unique<WellContributions>();
    }
    else{
        OPM_THROW(std::logic_error, "Invalid accelerator mode");
    }
}

template<class Scalar>
void WellContributions<Scalar>::
addMatrix([[maybe_unused]] MatrixType type,
          [[maybe_unused]] int* colIndices,
          [[maybe_unused]] Scalar* values,
          [[maybe_unused]] unsigned int val_size)
{
#if !HAVE_CUDA && !HAVE_OPENCL
    OPM_THROW(std::logic_error, "Error cannot add StandardWell matrix on GPU because neither CUDA nor OpenCL were found by cmake");
#endif

    if (!allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }

    this->APIaddMatrix(type, colIndices, values, val_size);

    if(MatrixType::B == type) {
        num_blocks_so_far += val_size;
        num_std_wells_so_far++;
    }
}

template<class Scalar>
void WellContributions<Scalar>::setBlockSize(unsigned int dim_, unsigned int dim_wells_)
{
    dim = dim_;
    dim_wells = dim_wells_;

    if(dim != 3 || dim_wells != 4){
        OPM_THROW(std::logic_error,
                  "WellContributions::setBlockSize error: "
                  "dim and dim_wells must be equal to 3 and 4, "
                  "respectively, otherwise the add well contributions "
                  "kernel won't work.\n");
    }
}

template<class Scalar>
void WellContributions<Scalar>::setVectorSize(unsigned N_)
{
    N = N_;
}

template<class Scalar>
void WellContributions<Scalar>::addNumBlocks(unsigned int numBlocks)
{
    if (allocated) {
        OPM_THROW(std::logic_error, "Error cannot add more sizes after allocated in WellContributions");
    }
    num_blocks += numBlocks;
    num_std_wells++;
}

template<class Scalar>
void WellContributions<Scalar>::alloc()
{
    if (num_std_wells > 0) {
        val_pointers.resize(num_std_wells+1);

        this->APIalloc();
        allocated = true;
    }
}

template<class Scalar>
void WellContributions<Scalar>::
addMultisegmentWellContribution(unsigned int dim_,
                                unsigned int dim_wells_,
                                unsigned int Mb,
                                std::vector<Scalar>& Bvalues,
                                std::vector<unsigned int>& BcolIndices,
                                std::vector<unsigned int>& BrowPointers,
                                unsigned int DnumBlocks,
                                Scalar* Dvalues,
                                UMFPackIndex* DcolPointers,
                                UMFPackIndex* DrowIndices,
                                std::vector<Scalar>& Cvalues)
{
    assert(dim==dim_);
    using MSW = MultisegmentWellContribution<Scalar>;
    multisegments.push_back(std::make_unique<MSW>(dim_,
                                                  dim_wells_,
                                                  Mb,
                                                  Bvalues,
                                                  BcolIndices,
                                                  BrowPointers,
                                                  DnumBlocks,
                                                  Dvalues,
                                                  DcolPointers,
                                                  DrowIndices,
                                                  Cvalues));
    ++num_ms_wells;
}

template class WellContributions<double>;

} //namespace Opm
