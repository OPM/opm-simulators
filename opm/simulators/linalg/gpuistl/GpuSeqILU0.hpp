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
#ifndef OPM_GPUSEQILU0_HPP
#define OPM_GPUSEQILU0_HPP

#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuSparseResource.hpp>
#include <opm/simulators/linalg/gpuistl/detail/fix_zero_diagonal.hpp>
#include <opm/simulators/linalg/is_gpu_operator.hpp>



namespace Opm::gpuistl
{
//! \brief Sequential ILU0 preconditioner on the GPU through the CuSparse library.
//!
//! This implementation calls the CuSparse functions, which in turn essentially
//! does a level decomposition to get some parallelism.
//!
//! \note This is not expected to be a fast preconditioner.
//!
//! \tparam M The matrix type to operate on
//! \tparam X Type of the update
//! \tparam Y Type of the defect
//! \tparam l Ignored. Just there to have the same number of template arguments
//!    as other preconditioners.
//!
//! \note We assume X and Y are both GpuVector<real_type>, but we leave them as template
//! arguments in case of future additions.
template <class M, class X, class Y, int l = 1>
class GpuSeqILU0 : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::remove_const<M>::type ;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    //! \brief Constructor.
    //!
    //!  Constructor gets all parameters to operate the prec.
    //! \param A The matrix to operate on.
    //! \param w The relaxation factor.
    //!
    template <typename T = M, typename std::enable_if_t<is_gpu_matrix_v<T>, int> = 0>
    GpuSeqILU0(const M& A, field_type w);


    //! \brief Constructor.
    //!
    //!  Constructor gets all parameters to operate the prec.
    //! \param A The matrix to operate on.
    //! \param w The relaxation factor.
    //!
    template <typename T = M, typename std::enable_if_t<!is_gpu_matrix_v<T>, int> = 0>
    GpuSeqILU0(const M& A, field_type w);


    //! \brief Prepare the preconditioner.
    //! \note Does nothing at the time being.
    virtual void pre(X& x, Y& b) override;

    //! \brief Apply the preconditoner.
    virtual void apply(X& v, const Y& d) override;

    //! \brief Post processing
    //! \note Does nothing at the moment
    virtual void post(X& x) override;

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const override;

    //! \brief Updates the matrix data.
    virtual void update() override;


    //! \returns false
    static constexpr bool shouldCallPre()
    {
        return false;
    }

    //! \returns false
    static constexpr bool shouldCallPost()
    {
        return false;
    }

    virtual bool hasPerfectUpdate() const override {
        return true;
    }


private:
    //! \brief Reference to the underlying matrix
    const M& m_underlyingMatrix;
    //! \brief The relaxation factor to use.
    field_type m_w;

    //! This is the storage for the LU composition.
    //! Initially this will have the values of A, but will be
    //! modified in the constructor to be the proper LU decomposition.
    GpuSparseMatrix<field_type> m_LU;

    GpuVector<field_type> m_temporaryStorage;


    detail::GpuSparseMatrixDescriptionPtr m_descriptionL;
    detail::GpuSparseMatrixDescriptionPtr m_descriptionU;
    detail::CuSparseResource<bsrsv2Info_t> m_infoL;
    detail::CuSparseResource<bsrsv2Info_t> m_infoU;
    detail::CuSparseResource<bsrilu02Info_t> m_infoM;

    std::unique_ptr<GpuVector<field_type>> m_buffer;
    detail::CuSparseHandle& m_cuSparseHandle;

    bool m_analysisDone = false;

    void analyzeMatrix();
    size_t findBufferSize();

    void createILU();

    void updateILUConfiguration();
};


// Case where matrix is a CPU matrix as input
template <class M, class X, class Y, int l>
template <typename T, typename std::enable_if_t<!is_gpu_matrix_v<T>, int>>
GpuSeqILU0<M, X, Y, l>::GpuSeqILU0(const M& A, field_type w)
    : m_underlyingMatrix(A)
    , m_w(w)
    , m_LU(GpuSparseMatrix<field_type>::fromMatrix(detail::makeMatrixWithNonzeroDiagonal(A)))
    , m_temporaryStorage(m_LU.N() * m_LU.blockSize())
    , m_descriptionL(detail::createLowerDiagonalDescription())
    , m_descriptionU(detail::createUpperDiagonalDescription())
    , m_cuSparseHandle(detail::CuSparseHandle::getInstance())
{
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_LU.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_LU.N(), A.N()));
    OPM_ERROR_IF(
        A[0][0].N() != m_LU.blockSize(),
        fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.", m_LU.blockSize(), A[0][0].N()));
    OPM_ERROR_IF(
        A.N() * A[0][0].N() != m_LU.dim(),
        fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.", m_LU.dim(), A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != m_LU.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_LU.nonzeroes(),
                             A.nonzeroes()));


    updateILUConfiguration();
}

// Case where matrix is a GPU matrix as input
template <class M, class X, class Y, int l>
template <typename T, typename std::enable_if_t<is_gpu_matrix_v<T>, int>>
GpuSeqILU0<M, X, Y, l>::GpuSeqILU0(const M& A, field_type w)
    : m_underlyingMatrix(A)
    , m_w(w)
    , m_LU(A)
    , m_temporaryStorage(m_LU.N() * m_LU.blockSize())
    , m_descriptionL(detail::createLowerDiagonalDescription())
    , m_descriptionU(detail::createUpperDiagonalDescription())
    , m_cuSparseHandle(detail::CuSparseHandle::getInstance())
{
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_LU.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_LU.N(), A.N()));
    OPM_ERROR_IF(A.nonzeroes() != m_LU.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_LU.nonzeroes(),
                             A.nonzeroes()));


    updateILUConfiguration();
}

} // end namespace Opm::gpuistl

#endif
