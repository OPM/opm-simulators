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
#ifndef OPM_GPUISTL_GPUOBLIVIOUSMPISENDER_HPP
#define OPM_GPUISTL_GPUOBLIVIOUSMPISENDER_HPP

#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSender.hpp>

#include <mpi.h>

#include <memory>

namespace Opm::gpuistl
{

/**
 * @brief Derived class of GPUSender that handles MPI calls that should NOT use GPU direct communicatoin
 * The implementation moves data fromthe GPU to the CPU and then sends it using regular MPI
 * @tparam field_type is float or double
 * @tparam block_size is the blocksize of the blockelements in the matrix
 * @tparam OwnerOverlapCopyCommunicationType is typically a Dune::LinearOperator::communication_type
 */
template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
class GPUObliviousMPISender : public GPUSender<field_type, OwnerOverlapCopyCommunicationType>
{
public:
    using X = GpuVector<field_type>;

    explicit GPUObliviousMPISender(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : GPUSender<field_type, OwnerOverlapCopyCommunicationType>(cpuOwnerOverlapCopy)
    {
    }

    void copyOwnerToAll(const X& source, X& dest) const override
    {
        // TODO: [perf] Can we reduce copying from the GPU here?
        // TODO: [perf] Maybe create a global buffer instead?
        auto sourceAsDuneVector = source.template asDuneBlockVector<block_size>();
        auto destAsDuneVector = dest.template asDuneBlockVector<block_size>();
        this->m_cpuOwnerOverlapCopy.copyOwnerToAll(sourceAsDuneVector, destAsDuneVector);
        dest.copyFromHost(destAsDuneVector);
    }

private:
    void initIndexSet() const override
    {
        // We need indices that we we will use in the project, dot and norm calls.
        // TODO: [premature perf] Can this be run once per instance? Or do we need to rebuild every time?
        const auto& pis = this->m_cpuOwnerOverlapCopy.indexSet();
        std::vector<int> indicesCopyOnCPU;
        std::vector<int> indicesOwnerCPU;
        for (const auto& index : pis) {
            if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::copy) {
                for (int component = 0; component < block_size; ++component) {
                    indicesCopyOnCPU.push_back(index.local().local() * block_size + component);
                }
            }

            if (index.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                for (int component = 0; component < block_size; ++component) {
                    indicesOwnerCPU.push_back(index.local().local() * block_size + component);
                }
            }
        }

        this->m_indicesCopy = std::make_unique<GpuVector<int>>(indicesCopyOnCPU);
        this->m_indicesOwner = std::make_unique<GpuVector<int>>(indicesOwnerCPU);
    }
};

} // namespace Opm::gpuistl

#endif // OPM_GPUISTL_GPUOBLIVIOUSMPISENDER_HPP