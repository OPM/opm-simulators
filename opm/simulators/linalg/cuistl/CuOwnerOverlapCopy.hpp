/*
  Copyright SINTEF AS 2022

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
#ifndef OPM_CUISTL_CUOWNEROVERLAPCOPY
#define OPM_CUISTL_CUOWNEROVERLAPCOPY
#include <dune/istl/owneroverlapcopy.hh>
#include <memory>
#include <mutex>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>

namespace Opm::cuistl
{
//! @brief CUDA compatiable variant of Dune::OwnerOverlapCopyCommunication
//!
//! @note This currently has very limited functionality.
template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
class CuOwnerOverlapCopy
{
public:
    using X = CuVector<field_type>;

    void dot(const X& x, const X& y, field_type& output) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });

        const auto dotAtRank = x.dot(y, *m_indicesOwner);
        output = m_cpuOwnerOverlapCopy.communicator().sum(dotAtRank);
    }

    field_type norm(const X& x) const
    {
        auto xDotX = field_type(0);
        this->dot(x, x, xDotX);

        using std::sqrt;
        return sqrt(xDotX);
    }

    void project(X& x) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });
        x.setZeroAtIndexSet(*m_indicesCopy);
    }

    void copyOwnerToAll(const X& source, X& dest) const
    {
        // TODO: [perf] Can we reduce copying from the GPU here?
        // TODO: [perf] Maybe create a global buffer instead?
        auto sourceAsDuneVector = source.template asDuneBlockVector<block_size>();
        auto destAsDuneVector = dest.template asDuneBlockVector<block_size>();
        m_cpuOwnerOverlapCopy.copyOwnerToAll(sourceAsDuneVector, destAsDuneVector);
        dest.copyFromHost(destAsDuneVector);
    }

    CuOwnerOverlapCopy(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : m_cpuOwnerOverlapCopy(cpuOwnerOverlapCopy)
    {
    }

private:
    mutable std::once_flag m_initializedIndices;
    const OwnerOverlapCopyCommunicationType& m_cpuOwnerOverlapCopy;

    mutable std::unique_ptr<CuVector<int>> m_indicesCopy;
    mutable std::unique_ptr<CuVector<int>> m_indicesOwner;


    void initIndexSet() const
    {
        // We need indices that we we will use in the project, dot and norm calls.
        // TODO: [premature perf] Can this be run once per instance? Or do we need to rebuild every time?
        const auto& pis = m_cpuOwnerOverlapCopy.indexSet();
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

        m_indicesCopy.reset(new CuVector<int>(indicesCopyOnCPU));
        m_indicesOwner.reset(new CuVector<int>(indicesOwnerCPU));
    }
};
} // namespace Opm::cuistl
#endif
