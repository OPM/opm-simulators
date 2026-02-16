/*
  Copyright 2022-2023 SINTEF AS, 2025 Equinor ASA

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
#ifndef OPM_GPUISTL_GPUOWNEROVERLAPCOPY_HPP
#define OPM_GPUISTL_GPUOWNEROVERLAPCOPY_HPP

#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/gpuistl/GpuAwareMPISender.hpp>
#include <opm/simulators/linalg/gpuistl/GpuObliviousMPISender.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <mpi.h>

#include <memory>

#ifdef OPEN_MPI
#if OPEN_MPI
#include "mpi-ext.h"
#endif
#endif

namespace Opm::gpuistl
{

/**
 * @brief CUDA compatiable variant of Dune::OwnerOverlapCopyCommunication
 *
 * This class can essentially be seen as an adapter around Dune::OwnerOverlapCopyCommunication, and should work as
 * a Dune::OwnerOverlapCopyCommunication on GpuVectors
 *
 * @note This currently only has the functionality to parallelize the linear solve.
 *
 * @tparam field_type should be a field_type supported by GpuVector (double, float)
 * @tparam OwnerOverlapCopyCommunicationType should mimic Dune::OwnerOverlapCopyCommunication.
 */
template <class field_type, class OwnerOverlapCopyCommunicationType>
class GpuOwnerOverlapCopy
{
public:
    using X = GpuVector<field_type>;

    explicit GpuOwnerOverlapCopy(std::shared_ptr<GPUSender<field_type, OwnerOverlapCopyCommunicationType>> sender)
        : m_sender(sender)
    {
    }

    void copyOwnerToAll(const X& source, X& dest) const
    {
        m_sender->copyOwnerToAll(source, dest);
    }

    void dot(const X& x, const X& y, field_type& output) const
    {
        m_sender->dot(x, y, output);
    }

    field_type norm(const X& x) const
    {
        return m_sender->norm(x);
    }

    void project(X& x) const
    {
        m_sender->project(x);
    }

    /**
     * @brief communicator returns the MPI communicator used by this GpuOwnerOverlapCopy
     * @return the MPI communicator
     */
    const ::Dune::Communication<MPI_Comm>& communicator() const
    {
        return m_sender->communicator();
    }

private:
    std::shared_ptr<GPUSender<field_type, OwnerOverlapCopyCommunicationType>> m_sender;
};

template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
std::shared_ptr<GpuOwnerOverlapCopy<field_type, OwnerOverlapCopyCommunicationType>>
makeGpuOwnerOverlapCopy(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
{

    const auto useGPUAwareMPI = Opm::Parameters::Get<Opm::Parameters::GpuAwareMpi>();
    const auto verifyGPUAwareMPI = Opm::Parameters::Get<Opm::Parameters::VerifyGpuAwareMpi>();
    std::shared_ptr<Opm::gpuistl::GPUSender<field_type, OwnerOverlapCopyCommunicationType>> gpuComm;

    if (useGPUAwareMPI) {
        if (verifyGPUAwareMPI) {

            // Temporary solution use the GPU Direct communication solely based on these prepcrosessor statements
            bool mpiSupportsCudaAwareAtCompileTime = false;
            bool mpiSupportsCudaAwareAtRunTime = false;

#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
            mpiSupportsCudaAwareAtCompileTime = true;
#endif /* MPIX_CUDA_AWARE_SUPPORT */

#if defined(MPIX_CUDA_AWARE_SUPPORT)
            if (1 == MPIX_Query_cuda_support()) {
                mpiSupportsCudaAwareAtRunTime = true;
            }
#endif /* MPIX_CUDA_AWARE_SUPPORT */

            if (!mpiSupportsCudaAwareAtCompileTime || !mpiSupportsCudaAwareAtRunTime) {
                OPM_THROW(std::runtime_error,
                          fmt::format(fmt::runtime("The GPU-aware MPI support is not available. "
                                      "CUDA aware support at compile time: {}, "
                                      "CUDA aware support at runtime: {}. "
                                      "Please check your MPI installation and the OPM configuration "
                                      "or run with --gpu-aware-mpi=false. If you are sure that your MPI "
                                      "implementation supports GPU aware MPI, you can disable this check "
                                      "by setting --verify-gpu-aware-mpi=false."),
                                      mpiSupportsCudaAwareAtCompileTime ? "yes" : "no",
                                      mpiSupportsCudaAwareAtRunTime ? "yes" : "no"));
            }
        }
        gpuComm = std::make_shared<
            Opm::gpuistl::GPUAwareMPISender<field_type, block_size, OwnerOverlapCopyCommunicationType>>(
            cpuOwnerOverlapCopy);

    } else {
        gpuComm = std::make_shared<
            Opm::gpuistl::GPUObliviousMPISender<field_type, block_size, OwnerOverlapCopyCommunicationType>>(
            cpuOwnerOverlapCopy);
    }

    using CudaCommunication = GpuOwnerOverlapCopy<field_type, OwnerOverlapCopyCommunicationType>;

    return std::make_shared<CudaCommunication>(gpuComm);
}
} // namespace Opm::gpuistl

#endif
