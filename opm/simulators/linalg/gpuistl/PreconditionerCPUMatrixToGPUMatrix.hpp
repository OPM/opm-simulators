/*
  Copyright 2025 Equinor ASA

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
#ifndef OPM_PRECONDITIONERCPUMATRIXTOGPUMATRIX_HPP
#define OPM_PRECONDITIONERCPUMATRIXTOGPUMATRIX_HPP
#include "opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp"
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/PreconditionerHolder.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditioner_should_call_post_pre.hpp>


namespace Opm::gpuistl
{
//!\brief Convert a CPU matrix to a GPU matrix and use a CUDA preconditioner on the GPU
//!
//! The use case is having a CPU matrix created by a CPU simulator being solved by a GPU BICGSTAB solver.
//!
//! \tparam X the domain type (should be on the GPU). Typicall a GpuVector
//! \tparam Y the range type (should be on the GPU). Typicall a GpuVector
//! \tparam CudaPreconditionerType the preconditioner taking GpuVector<real_type>  as arguments to apply
//!         **and** expecting a GPU matrix as **the first argument** to its constructor.
//! \tparam CPUMatrixType the type of the CPU matrix. This is the type that will be copied to the GPU.
template <class X, class Y, class CudaPreconditionerType, class CPUMatrixType>
class PreconditionerCPUMatrixToGPUMatrix : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;


    template <typename... Args>
    explicit PreconditionerCPUMatrixToGPUMatrix(const CPUMatrixType& A, Args&&... args)
        : m_cpuMatrix(A)
        , m_gpuMatrix(GpuSparseMatrixWrapper<field_type>::fromMatrix(A))
        , m_underlyingPreconditioner(m_gpuMatrix, std::forward<Args>(args)...)
    {
    }


    //! \brief Prepare the preconditioner.
    //!
    //! Currently not supported.
    void pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b) override
    {
        static_assert(!detail::shouldCallPreconditionerPre<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::pre().");
    }


    //! \brief Apply the preconditoner.
    //!
    //! \copydoc Preconditioner::apply(X&,const Y&)
    void apply(X& v, const Y& d) override
    {
        m_underlyingPreconditioner.apply(v, d);
    }


    //! \brief Clean up.
    //!
    //! Currently not supported.
    void post([[maybe_unused]] X& x) override
    {
        static_assert(!detail::shouldCallPreconditionerPost<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::post().");
    }


    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return m_underlyingPreconditioner.category();
    }

    //! Copies the CPU matrix to the GPU matrix
    //! then calls the GPU preconditioner update function
    void update() override
    {
        m_gpuMatrix.updateNonzeroValues(m_cpuMatrix, true);
        m_underlyingPreconditioner.update();
    }

    static constexpr bool shouldCallPre()
    {
        return detail::shouldCallPreconditionerPost<CudaPreconditionerType>();
    }
    static constexpr bool shouldCallPost()
    {
        return detail::shouldCallPreconditionerPre<CudaPreconditionerType>();
    }

    virtual bool hasPerfectUpdate() const override
    {
        return m_underlyingPreconditioner.hasPerfectUpdate();
    }

private:
    const CPUMatrixType& m_cpuMatrix;
    GpuSparseMatrixWrapper<field_type> m_gpuMatrix;

    //! \brief the underlying preconditioner to use
    CudaPreconditionerType m_underlyingPreconditioner;
};
} // end namespace Opm::gpuistl

#endif // OPM_PRECONDITIONERCPUMATRIXTOGPUMATRIX_HPP
