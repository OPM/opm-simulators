/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUISTL_CUOWNEROVERLAPCOPY_HPP
#define OPM_CUISTL_CUOWNEROVERLAPCOPY_HPP
#include <dune/istl/owneroverlapcopy.hh>
#include <memory>
#include <mutex>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>

namespace Opm::cuistl
{
/**
 * @brief CUDA compatiable variant of Dune::OwnerOverlapCopyCommunication
 *
 * This class can essentially be seen as an adapter around Dune::OwnerOverlapCopyCommunication, and should work as
 * a Dune::OwnerOverlapCopyCommunication on CuVectors
 *
 *
 * @note This currently only has the functionality to parallelize the linear solve.
 *
 * @tparam field_type should be a field_type supported by CuVector (double, float)
 * @tparam block_size the block size used (this is relevant for say figuring out the correct indices)
 * @tparam OwnerOverlapCopyCommunicationType should mimic Dune::OwnerOverlapCopyCommunication.
 */
template <class field_type, int block_size, class OwnerOverlapCopyCommunicationType>
class CuOwnerOverlapCopy
{
public:
    using X = CuVector<field_type>;

    CuOwnerOverlapCopy(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : m_cpuOwnerOverlapCopy(cpuOwnerOverlapCopy)
    {
    }
    /**
     * @brief dot will carry out the dot product between x and y on the owned indices, then sum up the result across MPI
     * processes.
     * @param[out] output result will be stored here
     *
     * @note This uses the same interface as its DUNE equivalent.
     */
    void dot(const X& x, const X& y, field_type& output) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });

        const auto dotAtRank = x.dot(y, *m_indicesOwner);
        output = m_cpuOwnerOverlapCopy.communicator().sum(dotAtRank);
    }

    /**
     * @brief norm computes the l^2-norm of x across processes.
     *
     * This will compute the dot product of x with itself on owned indices, then
     * sum the result across process and return the square root of the sum.
     */
    field_type norm(const X& x) const
    {
        auto xDotX = field_type(0);
        this->dot(x, x, xDotX);

        using std::sqrt;
        return sqrt(xDotX);
    }

    /**
     * @brief project will project x to the owned subspace
     *
     * For each component i which is not owned, x_i will be set to 0
     *
     * @param[inout] x the vector to project
     */
    void project(X& x) const
    {
        std::call_once(m_initializedIndices, [&]() { initIndexSet(); });
        x.setZeroAtIndexSet(*m_indicesCopy);
    }

    /**
     * @brief copyOwnerToAll will copy source to the CPU, then call OwnerOverlapCopyCommunicationType::copyOwnerToAll on
     * the copied data, and copy the result back to the GPU
     * @param[in] source
     * @param[out] dest
     */
    void copyOwnerToAll(const X& source, X& dest) const
    {
        // TODO: [perf] Can we reduce copying from the GPU here?
        // TODO: [perf] Maybe create a global buffer instead?
        auto sourceAsDuneVector = source.template asDuneBlockVector<block_size>();
        auto destAsDuneVector = dest.template asDuneBlockVector<block_size>();
        m_cpuOwnerOverlapCopy.copyOwnerToAll(sourceAsDuneVector, destAsDuneVector);
        dest.copyFromHost(destAsDuneVector);
    }



private:
    const OwnerOverlapCopyCommunicationType& m_cpuOwnerOverlapCopy;

    // Used to call the initIndexSet. Note that this is kind of a
    // premature optimization, in the sense that we could just initialize these indices
    // always, but they are not always used.
    mutable std::once_flag m_initializedIndices;
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

        m_indicesCopy = std::make_unique<CuVector<int>>(indicesCopyOnCPU);
        m_indicesOwner = std::make_unique<CuVector<int>>(indicesOwnerCPU);
    }
};
} // namespace Opm::cuistl
#endif
