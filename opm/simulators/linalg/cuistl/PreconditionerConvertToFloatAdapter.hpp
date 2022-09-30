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
#ifndef OPM_PRECONDITIONERCONVERTOFLOATADAPTER_HEADER_INCLUDED
#define OPM_PRECONDITIONERCONVERTOFLOATADAPTER_HEADER_INCLUDED
#include <cusparse.h>
#include <dune/common/simd.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/impl/preconditioner_should_call_post_pre.hpp>


namespace Opm::cuistl
{
//! \brief Converts double to single precision to benchmark single precision preconditioners
//!
//! \note This is not a fast conversion, it is simply meant to benchmark the potential of some
//!       preconditioners on consumer grade GPUs where the double precision performance is often
//!       artificially limited.
//!
//! To use this, use something like the following pseudo code
//! \code{.cpp}
//! auto adapter = PreconditionerFloatAdapter<args>(matrixDouble);
//! auto underlyingPreconditioner = SomePreconditioner(adapter.getConvertedMatrix());
//! adapter.setUnderlyingPreconditioner(underlyingPreconditioner);
//!
//! // Now call adapter.apply(...) with double vectors
//! \endcode
template <class CudaPreconditionerType, class M, class X, class Y, int l = 1>
class PreconditionerConvertToFloatAdapter : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename std::remove_const<M>::type matrix_type;

    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    typedef typename X::field_type block_typefield_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::SimdScalar<field_type> scalar_field_type;

    typedef typename CudaPreconditionerType::domain_type domain_type_to;
    //! \brief The range type of the preconditioner.
    typedef typename CudaPreconditionerType::range_type range_type_to;
    //! \brief The field type of the preconditioner.
    typedef typename domain_type_to::field_type field_type_to;
    using block_type = typename domain_type::block_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::SimdScalar<field_type_to> scalar_field_type_to;
    using XTo = Dune::BlockVector<Dune::FieldVector<scalar_field_type_to, block_type::dimension>>;
    using YTo = Dune::BlockVector<Dune::FieldVector<scalar_field_type_to, block_type::dimension>>;
    using matrix_type_to = typename Dune::BCRSMatrix<
        Dune::FieldMatrix<scalar_field_type_to, block_type::dimension, block_type::dimension>>;

    //! \brief Constructor.
    //!
    //! \param A The matrix to operate on.
    PreconditionerConvertToFloatAdapter(const M& matrix)
        : m_matrix(matrix)
        , m_convertedMatrix(createConvertedMatrix())
    {
    }

    //! \brief Not used at the moment
    virtual void pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b) override
    {
        static_assert(!impl::shouldCallPreconditionerPre<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::pre().");
    }

    //! \brief Apply the preconditoner.
    virtual void apply(X& v, const Y& d) override
    {
        XTo convertedV(v.N());
        for (size_t i = 0; i < v.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                // This is probably unnecessary, but doing it anyway:
                convertedV[i][j] = scalar_field_type_to(v[i][j]);
            }
        }
        YTo convertedD(d.N());
        for (size_t i = 0; i < d.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                convertedD[i][j] = scalar_field_type_to(d[i][j]);
            }
        }

        m_underlyingPreconditioner->apply(convertedV, convertedD);

        for (size_t i = 0; i < v.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                v[i][j] = scalar_field_type(convertedV[i][j]);
            }
        }
    }

    //! \brief Not used at the moment
    virtual void post([[maybe_unused]] X& x) override
    {
        static_assert(!impl::shouldCallPreconditionerPost<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::post().");
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return m_underlyingPreconditioner->category();
    }

    virtual void update() override
    {
        updateMatrix();
        m_underlyingPreconditioner->update();
    }

    const matrix_type_to& getConvertedMatrix() const
    {
        return m_convertedMatrix;
    }

    void setUnderlyingPreconditioner(const std::shared_ptr<CudaPreconditionerType>& conditioner)
    {
        m_underlyingPreconditioner = conditioner;
    }


private:
    void updateMatrix()
    {
        const auto nnz = m_matrix.nonzeroes() * m_matrix[0][0].N() * m_matrix[0][0].N();
        const auto dataPointerIn = static_cast<const scalar_field_type*>(&((m_matrix[0][0][0][0])));
        auto dataPointerOut = static_cast<scalar_field_type_to*>(&((m_convertedMatrix[0][0][0][0])));

        std::vector<scalar_field_type_to> buffer(nnz, 0);
        for (size_t i = 0; i < nnz; ++i) {
            dataPointerOut[i] = scalar_field_type_to(dataPointerIn[i]);
        }
    }
    matrix_type_to createConvertedMatrix()
    {
        // TODO: Check if this whole conversion can be done more efficiently.
        const auto N = m_matrix.N();
        matrix_type_to matrixBuilder(N, N, m_matrix.nonzeroes(), matrix_type_to::row_wise);
        {
            auto rowIn = m_matrix.begin();
            for (auto rowOut = matrixBuilder.createbegin(); rowOut != matrixBuilder.createend(); ++rowOut) {
                for (auto column = rowIn->begin(); column != rowIn->end(); ++column) {
                    rowOut.insert(column.index());
                }
                ++rowIn;
            }
        }

        for (auto row = m_matrix.begin(); row != m_matrix.end(); ++row) {
            for (auto column = row->begin(); column != row->end(); ++column) {
                for (size_t i = 0; i < block_type::dimension; ++i) {
                    for (size_t j = 0; j < block_type::dimension; ++j) {
                        matrixBuilder[row.index()][column.index()][i][j]
                            = scalar_field_type_to(m_matrix[row.index()][column.index()][i][j]);
                    }
                }
            }
        }

        return matrixBuilder;
    }
    const M& m_matrix;
    matrix_type_to m_convertedMatrix;
    //! \brief the underlying preconditioner to use
    std::shared_ptr<CudaPreconditionerType> m_underlyingPreconditioner;
};
} // end namespace Opm::cuistl

#endif
