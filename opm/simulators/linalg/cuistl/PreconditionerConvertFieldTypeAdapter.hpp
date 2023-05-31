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
#ifndef OPM_PRECONDITIONERCONVERTOFLOATADAPTER_HPP
#define OPM_PRECONDITIONERCONVERTOFLOATADAPTER_HPP
#include <cusparse.h>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/preconditioner_should_call_post_pre.hpp>


namespace Opm::cuistl
{
//! \brief Converts the field type (eg. double to float) to benchmark single precision preconditioners
//!
//! \note This is not a fast conversion, it is simply meant to benchmark the potential of some
//!       preconditioners on consumer grade GPUs where the double precision performance is often
//!       artificially limited.
//!
//! \note In theory this can do any field_type conversion that is meaningful, but it is only tested
//! on double to float conversion.
//!
//! \note Remember to set the underlying preconditioner with setUnderlyingPreconditioner (should use the matrix from
//! getConvertedMatrix())
//!
//! \note One could simply change the constructor design by accepting a creator function for the underlying
//! preconditioner. For the current use cases this is however not needed.
//!
//! To use this, use something like the following code:
//! \code{.cpp}
//! #include <opm/simulators/linalg/cuistl/PreconditionerConvertFieldTypeAdapter.hpp>
//! #include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
//!
//! using XDouble = Dune::BlockVector<Dune::FieldVector<double, 2>>;
//! using MDouble = Dune::FieldMatrix<double, 2, 2>;
//! using SpMatrixDouble = Dune::BCRSMatrix<MDouble>;
//! using XFloat = Dune::BlockVector<Dune::FieldVector<float, 2>>;
//! using MFloat = Dune::FieldMatrix<float, 2, 2>;
//! using SpMatrixFloat = Dune::BCRSMatrix<MFloat>;
//!
//! template<class ParallelInfo>
//! void applyILU0AsFloat(const MDouble& matrix, const XDouble& x, XDouble& y) {
//!
//!     using FloatILU0 = typename Opm::ParallelOverlappingILU0<MFloat, XFloat, XFloat, ParallelInfo>;
//!     using DoubleToFloatConverter = typename Opm::cuistl::PreconditionerConvertFieldTypeAdapter<FloatILU0, MDouble,
//!     XDouble, XDouble>;
//!
//!     // Note that we do not need to make a new instance for every invocation, this
//!     // is just done for this example
//!     auto doubleToFloatConverter = DoubleToFloatConverter(matrix);
//!     const auto& convertedMatrix = doubleToFloatConverter.getConvertedMatrix();
//!
//!     auto floatILU0 = std::make_shared<FloatILU0>(convertedMatrix, 0, 1.0, 0);
//!
//!     doubleToFloatConverter.setUnderlyingPreconditioner(floatILU0);
//!
//!     // This will convert x and y to float, then call floatILU0.apply on the converted arguments
//!     doubleToFloatConverter.apply(x, y);
//! }
//!
//! \endcode
template <class CudaPreconditionerType, class M, class X, class Y, int l = 1>
class PreconditionerConvertFieldTypeAdapter : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::remove_const<M>::type;

    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;


    using domain_type_to = typename CudaPreconditionerType::domain_type;
    //! \brief The range type of the preconditioner.
    using range_type_to = typename CudaPreconditionerType::range_type;
    //! \brief The field type of the preconditioner.
    using field_type_to = typename domain_type_to::field_type;


    using block_type = typename domain_type::block_type;

    using XTo = Dune::BlockVector<Dune::FieldVector<field_type_to, block_type::dimension>>;
    using YTo = Dune::BlockVector<Dune::FieldVector<field_type_to, block_type::dimension>>;
    using matrix_type_to =
        typename Dune::BCRSMatrix<Dune::FieldMatrix<field_type_to, block_type::dimension, block_type::dimension>>;

    //! \brief Constructor.
    //!
    //! \param A The matrix to operate on.
    //!
    //! \note After the PreconditionerConvertFieldTypeAdapter you can get the converted matrix
    //! by calling getConvertedMatrix(), which in turn can be used to create the underlying preconditioner.
    //! Once the underlying preconditioner has been called, this must be supplied to setUnderlyingPreconditioner.
    explicit PreconditionerConvertFieldTypeAdapter(const M& matrix)
        : m_matrix(matrix)
        , m_convertedMatrix(createConvertedMatrix())
    {
    }

    //! \brief Not used at the moment
    virtual void pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b) override
    {
        static_assert(!detail::shouldCallPreconditionerPre<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::pre().");
    }

    //! \brief Apply the preconditoner.
    virtual void apply(X& v, const Y& d) override
    {
        OPM_ERROR_IF(!m_underlyingPreconditioner,
                     "You need to set the underlying preconditioner with setUnderlyingPreconditioner.");
        XTo convertedV(v.N());
        for (size_t i = 0; i < v.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                // This is probably unnecessary, but doing it anyway:
                convertedV[i][j] = field_type_to(v[i][j]);
            }
        }
        YTo convertedD(d.N());
        for (size_t i = 0; i < d.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                convertedD[i][j] = field_type_to(d[i][j]);
            }
        }

        m_underlyingPreconditioner->apply(convertedV, convertedD);

        for (size_t i = 0; i < v.N(); ++i) {
            for (size_t j = 0; j < block_type::dimension; ++j) {
                v[i][j] = field_type(convertedV[i][j]);
            }
        }
    }

    //! \brief Not used at the moment
    virtual void post([[maybe_unused]] X& x) override
    {
        static_assert(!detail::shouldCallPreconditionerPost<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::post().");
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const override
    {
        return m_underlyingPreconditioner->category();
    }

    virtual void update() override
    {
        OPM_ERROR_IF(!m_underlyingPreconditioner,
                     "You need to set the underlying preconditioner with setUnderlyingPreconditioner.");
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
        const auto dataPointerIn = static_cast<const field_type*>(&((m_matrix[0][0][0][0])));
        auto dataPointerOut = static_cast<field_type_to*>(&((m_convertedMatrix[0][0][0][0])));

        std::vector<field_type_to> buffer(nnz, 0);
        for (size_t i = 0; i < nnz; ++i) {
            dataPointerOut[i] = field_type_to(dataPointerIn[i]);
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
                            = field_type_to(m_matrix[row.index()][column.index()][i][j]);
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
