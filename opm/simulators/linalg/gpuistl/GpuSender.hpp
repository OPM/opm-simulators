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
#ifndef OPM_GPUISTL_GPUSENDER_HPP
#define OPM_GPUISTL_GPUSENDER_HPP

#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <mpi.h>

#include <memory>
#include <mutex>

namespace Opm::gpuistl
{

/**
 * @brief GPUSender is a wrapper class for classes which will implement copOwnerToAll
 * This is implemented with the intention of creating communicators with generic GPUSender
 * To hide implementation that will either use GPU aware MPI or not
 * @tparam field_type is float or double
 * @tparam OwnerOverlapCopyCommunicationType is typically a Dune::LinearOperator::communication_type
 */
template <class field_type, class OwnerOverlapCopyCommunicationType>
class GPUSender
{
public:
    using X = GpuVector<field_type>;

    GPUSender(const OwnerOverlapCopyCommunicationType& cpuOwnerOverlapCopy)
        : m_cpuOwnerOverlapCopy(cpuOwnerOverlapCopy)
    {
    }

    virtual ~GPUSender() = default;

    /**
     * @brief copyOwnerToAll will copy the data in source to all processes. 
     *
     * @note Depending on the implementation, this may or may not use GPU aware MPI.
     *       If it does not use GPU aware MPI, the data will be copied to the
     *       CPU before the communication.
     * @param[in] source
     * @param[out] dest
     */
    virtual void copyOwnerToAll(const X& source, X& dest) const = 0;
    virtual void initIndexSet() const = 0;

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
     * @brief dot will carry out the dot product between x and y on the owned indices, then sum up the result across MPI
     * processes.
     * @param[in] x First vector in dot product
     * @param[in] y Second vector in dot product
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
        dot(x, x, xDotX);

        // using std::sqrt;
        return std::sqrt(xDotX);
    }


    /**
     * @brief communicator returns the MPI communicator used by this GPUSender
     * @return the MPI communicator
     */
    const ::Dune::Communication<MPI_Comm>& communicator() const
    {
        return m_cpuOwnerOverlapCopy.communicator();
    }

protected:
    // Used to call the initIndexSet. Note that this is kind of a
    // premature optimization, in the sense that we could just initialize these indices
    // always, but they are not always used.
    mutable std::once_flag m_initializedIndices;
    mutable std::unique_ptr<GpuVector<int>> m_indicesOwner;
    mutable std::unique_ptr<GpuVector<int>> m_indicesCopy;
    const OwnerOverlapCopyCommunicationType& m_cpuOwnerOverlapCopy;
};

} // namespace Opm::gpuistl

#endif // OPM_GPUISTL_GPUSENDER_HPP
