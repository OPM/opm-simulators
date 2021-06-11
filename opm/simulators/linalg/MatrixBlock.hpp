/*
  Copyright 2016 IRIS AS
  Copyright 2019 NORCE

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

#ifndef OPM_MATRIX_BLOCK_HEADER_INCLUDED
#define OPM_MATRIX_BLOCK_HEADER_INCLUDED

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/istl/matrixutils.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/superlu.hh>

namespace Dune
{
namespace FMatrixHelp {
//! invert 4x4 Matrix without changing the original matrix
template <typename K>
static inline K invertMatrix(const FieldMatrix<K,4,4>& matrix, FieldMatrix<K,4,4>& inverse)
{
    inverse[0][0] = matrix[1][1] * matrix[2][2] * matrix[3][3] -
            matrix[1][1] * matrix[2][3] * matrix[3][2] -
            matrix[2][1] * matrix[1][2] * matrix[3][3] +
            matrix[2][1] * matrix[1][3] * matrix[3][2] +
            matrix[3][1] * matrix[1][2] * matrix[2][3] -
            matrix[3][1] * matrix[1][3] * matrix[2][2];

    inverse[1][0] = -matrix[1][0] * matrix[2][2] * matrix[3][3] +
            matrix[1][0] * matrix[2][3] * matrix[3][2] +
            matrix[2][0] * matrix[1][2] * matrix[3][3] -
            matrix[2][0] * matrix[1][3] * matrix[3][2] -
            matrix[3][0] * matrix[1][2] * matrix[2][3] +
            matrix[3][0] * matrix[1][3] * matrix[2][2];

    inverse[2][0] = matrix[1][0] * matrix[2][1] * matrix[3][3] -
            matrix[1][0] * matrix[2][3] * matrix[3][1] -
            matrix[2][0] * matrix[1][1] * matrix[3][3] +
            matrix[2][0] * matrix[1][3] * matrix[3][1] +
            matrix[3][0] * matrix[1][1] * matrix[2][3] -
            matrix[3][0] * matrix[1][3] * matrix[2][1];

    inverse[3][0] = -matrix[1][0] * matrix[2][1] * matrix[3][2] +
            matrix[1][0] * matrix[2][2] * matrix[3][1] +
            matrix[2][0] * matrix[1][1] * matrix[3][2] -
            matrix[2][0] * matrix[1][2] * matrix[3][1] -
            matrix[3][0] * matrix[1][1] * matrix[2][2] +
            matrix[3][0] * matrix[1][2] * matrix[2][1];

    inverse[0][1]= -matrix[0][1]  * matrix[2][2] * matrix[3][3] +
            matrix[0][1] * matrix[2][3] * matrix[3][2] +
            matrix[2][1] * matrix[0][2] * matrix[3][3] -
            matrix[2][1] * matrix[0][3] * matrix[3][2] -
            matrix[3][1] * matrix[0][2] * matrix[2][3] +
            matrix[3][1] * matrix[0][3] * matrix[2][2];

    inverse[1][1] = matrix[0][0] * matrix[2][2] * matrix[3][3] -
            matrix[0][0] * matrix[2][3] * matrix[3][2] -
            matrix[2][0] * matrix[0][2] * matrix[3][3] +
            matrix[2][0] * matrix[0][3] * matrix[3][2] +
            matrix[3][0] * matrix[0][2] * matrix[2][3] -
            matrix[3][0] * matrix[0][3] * matrix[2][2];

    inverse[2][1] = -matrix[0][0] * matrix[2][1] * matrix[3][3] +
            matrix[0][0] * matrix[2][3] * matrix[3][1] +
            matrix[2][0] * matrix[0][1] * matrix[3][3] -
            matrix[2][0] * matrix[0][3] * matrix[3][1] -
            matrix[3][0] * matrix[0][1] * matrix[2][3] +
            matrix[3][0] * matrix[0][3] * matrix[2][1];

    inverse[3][1] = matrix[0][0] * matrix[2][1] * matrix[3][2] -
            matrix[0][0] * matrix[2][2] * matrix[3][1] -
            matrix[2][0] * matrix[0][1] * matrix[3][2] +
            matrix[2][0] * matrix[0][2] * matrix[3][1] +
            matrix[3][0] * matrix[0][1] * matrix[2][2] -
            matrix[3][0] * matrix[0][2] * matrix[2][1];

    inverse[0][2] = matrix[0][1] * matrix[1][2] * matrix[3][3] -
            matrix[0][1] * matrix[1][3] * matrix[3][2] -
            matrix[1][1] * matrix[0][2] * matrix[3][3] +
            matrix[1][1] * matrix[0][3] * matrix[3][2] +
            matrix[3][1] * matrix[0][2] * matrix[1][3] -
            matrix[3][1] * matrix[0][3] * matrix[1][2];

    inverse[1][2] = -matrix[0][0]  * matrix[1][2] * matrix[3][3] +
            matrix[0][0] * matrix[1][3] * matrix[3][2] +
            matrix[1][0] * matrix[0][2] * matrix[3][3] -
            matrix[1][0] * matrix[0][3] * matrix[3][2] -
            matrix[3][0] * matrix[0][2] * matrix[1][3] +
            matrix[3][0] * matrix[0][3] * matrix[1][2];

    inverse[2][2] = matrix[0][0] * matrix[1][1] * matrix[3][3] -
            matrix[0][0] * matrix[1][3] * matrix[3][1] -
            matrix[1][0] * matrix[0][1] * matrix[3][3] +
            matrix[1][0] * matrix[0][3] * matrix[3][1] +
            matrix[3][0] * matrix[0][1] * matrix[1][3] -
            matrix[3][0] * matrix[0][3] * matrix[1][1];

    inverse[3][2] = -matrix[0][0] * matrix[1][1] * matrix[3][2] +
            matrix[0][0] * matrix[1][2] * matrix[3][1] +
            matrix[1][0] * matrix[0][1] * matrix[3][2] -
            matrix[1][0] * matrix[0][2] * matrix[3][1] -
            matrix[3][0] * matrix[0][1] * matrix[1][2] +
            matrix[3][0] * matrix[0][2] * matrix[1][1];

    inverse[0][3] = -matrix[0][1] * matrix[1][2] * matrix[2][3] +
            matrix[0][1] * matrix[1][3] * matrix[2][2] +
            matrix[1][1] * matrix[0][2] * matrix[2][3] -
            matrix[1][1] * matrix[0][3] * matrix[2][2] -
            matrix[2][1] * matrix[0][2] * matrix[1][3] +
            matrix[2][1] * matrix[0][3] * matrix[1][2];

    inverse[1][3] = matrix[0][0] * matrix[1][2] * matrix[2][3] -
            matrix[0][0] * matrix[1][3] * matrix[2][2] -
            matrix[1][0] * matrix[0][2] * matrix[2][3] +
            matrix[1][0] * matrix[0][3] * matrix[2][2] +
            matrix[2][0] * matrix[0][2] * matrix[1][3] -
            matrix[2][0] * matrix[0][3] * matrix[1][2];

    inverse[2][3] = -matrix[0][0] * matrix[1][1] * matrix[2][3] +
            matrix[0][0] * matrix[1][3] * matrix[2][1] +
            matrix[1][0] * matrix[0][1] * matrix[2][3] -
            matrix[1][0] * matrix[0][3] * matrix[2][1] -
            matrix[2][0] * matrix[0][1] * matrix[1][3] +
            matrix[2][0] * matrix[0][3] * matrix[1][1];

    inverse[3][3] = matrix[0][0] * matrix[1][1] * matrix[2][2] -
            matrix[0][0] * matrix[1][2] * matrix[2][1] -
            matrix[1][0] * matrix[0][1] * matrix[2][2] +
            matrix[1][0] * matrix[0][2] * matrix[2][1] +
            matrix[2][0] * matrix[0][1] * matrix[1][2] -
            matrix[2][0] * matrix[0][2] * matrix[1][1];

    K det = matrix[0][0] * inverse[0][0] + matrix[0][1] * inverse[1][0] +
            matrix[0][2] * inverse[2][0] + matrix[0][3] * inverse[3][0];

    // return identity for singular or nearly singular matrices.
    if (std::abs(det) < 1e-40) {
        for (int i = 0; i < 4; ++i){
            inverse[i][i] = 1.0;
        }
        return 1.0;
    }
    K inv_det = 1.0 / det;
    inverse *= inv_det;

    return det;
}

template <typename K>
static inline K invertMatrix(const DynamicMatrix<K>& matrix, DynamicMatrix<K>& inverse)
{
    // this function is only for 4 X 4 matrix
    assert (matrix.rows() == 4);

    inverse[0][0] = matrix[1][1] * matrix[2][2] * matrix[3][3] -
            matrix[1][1] * matrix[2][3] * matrix[3][2] -
            matrix[2][1] * matrix[1][2] * matrix[3][3] +
            matrix[2][1] * matrix[1][3] * matrix[3][2] +
            matrix[3][1] * matrix[1][2] * matrix[2][3] -
            matrix[3][1] * matrix[1][3] * matrix[2][2];

    inverse[1][0] = -matrix[1][0] * matrix[2][2] * matrix[3][3] +
            matrix[1][0] * matrix[2][3] * matrix[3][2] +
            matrix[2][0] * matrix[1][2] * matrix[3][3] -
            matrix[2][0] * matrix[1][3] * matrix[3][2] -
            matrix[3][0] * matrix[1][2] * matrix[2][3] +
            matrix[3][0] * matrix[1][3] * matrix[2][2];

    inverse[2][0] = matrix[1][0] * matrix[2][1] * matrix[3][3] -
            matrix[1][0] * matrix[2][3] * matrix[3][1] -
            matrix[2][0] * matrix[1][1] * matrix[3][3] +
            matrix[2][0] * matrix[1][3] * matrix[3][1] +
            matrix[3][0] * matrix[1][1] * matrix[2][3] -
            matrix[3][0] * matrix[1][3] * matrix[2][1];

    inverse[3][0] = -matrix[1][0] * matrix[2][1] * matrix[3][2] +
            matrix[1][0] * matrix[2][2] * matrix[3][1] +
            matrix[2][0] * matrix[1][1] * matrix[3][2] -
            matrix[2][0] * matrix[1][2] * matrix[3][1] -
            matrix[3][0] * matrix[1][1] * matrix[2][2] +
            matrix[3][0] * matrix[1][2] * matrix[2][1];

    inverse[0][1]= -matrix[0][1]  * matrix[2][2] * matrix[3][3] +
            matrix[0][1] * matrix[2][3] * matrix[3][2] +
            matrix[2][1] * matrix[0][2] * matrix[3][3] -
            matrix[2][1] * matrix[0][3] * matrix[3][2] -
            matrix[3][1] * matrix[0][2] * matrix[2][3] +
            matrix[3][1] * matrix[0][3] * matrix[2][2];

    inverse[1][1] = matrix[0][0] * matrix[2][2] * matrix[3][3] -
            matrix[0][0] * matrix[2][3] * matrix[3][2] -
            matrix[2][0] * matrix[0][2] * matrix[3][3] +
            matrix[2][0] * matrix[0][3] * matrix[3][2] +
            matrix[3][0] * matrix[0][2] * matrix[2][3] -
            matrix[3][0] * matrix[0][3] * matrix[2][2];

    inverse[2][1] = -matrix[0][0] * matrix[2][1] * matrix[3][3] +
            matrix[0][0] * matrix[2][3] * matrix[3][1] +
            matrix[2][0] * matrix[0][1] * matrix[3][3] -
            matrix[2][0] * matrix[0][3] * matrix[3][1] -
            matrix[3][0] * matrix[0][1] * matrix[2][3] +
            matrix[3][0] * matrix[0][3] * matrix[2][1];

    inverse[3][1] = matrix[0][0] * matrix[2][1] * matrix[3][2] -
            matrix[0][0] * matrix[2][2] * matrix[3][1] -
            matrix[2][0] * matrix[0][1] * matrix[3][2] +
            matrix[2][0] * matrix[0][2] * matrix[3][1] +
            matrix[3][0] * matrix[0][1] * matrix[2][2] -
            matrix[3][0] * matrix[0][2] * matrix[2][1];

    inverse[0][2] = matrix[0][1] * matrix[1][2] * matrix[3][3] -
            matrix[0][1] * matrix[1][3] * matrix[3][2] -
            matrix[1][1] * matrix[0][2] * matrix[3][3] +
            matrix[1][1] * matrix[0][3] * matrix[3][2] +
            matrix[3][1] * matrix[0][2] * matrix[1][3] -
            matrix[3][1] * matrix[0][3] * matrix[1][2];

    inverse[1][2] = -matrix[0][0]  * matrix[1][2] * matrix[3][3] +
            matrix[0][0] * matrix[1][3] * matrix[3][2] +
            matrix[1][0] * matrix[0][2] * matrix[3][3] -
            matrix[1][0] * matrix[0][3] * matrix[3][2] -
            matrix[3][0] * matrix[0][2] * matrix[1][3] +
            matrix[3][0] * matrix[0][3] * matrix[1][2];

    inverse[2][2] = matrix[0][0] * matrix[1][1] * matrix[3][3] -
            matrix[0][0] * matrix[1][3] * matrix[3][1] -
            matrix[1][0] * matrix[0][1] * matrix[3][3] +
            matrix[1][0] * matrix[0][3] * matrix[3][1] +
            matrix[3][0] * matrix[0][1] * matrix[1][3] -
            matrix[3][0] * matrix[0][3] * matrix[1][1];

    inverse[3][2] = -matrix[0][0] * matrix[1][1] * matrix[3][2] +
            matrix[0][0] * matrix[1][2] * matrix[3][1] +
            matrix[1][0] * matrix[0][1] * matrix[3][2] -
            matrix[1][0] * matrix[0][2] * matrix[3][1] -
            matrix[3][0] * matrix[0][1] * matrix[1][2] +
            matrix[3][0] * matrix[0][2] * matrix[1][1];

    inverse[0][3] = -matrix[0][1] * matrix[1][2] * matrix[2][3] +
            matrix[0][1] * matrix[1][3] * matrix[2][2] +
            matrix[1][1] * matrix[0][2] * matrix[2][3] -
            matrix[1][1] * matrix[0][3] * matrix[2][2] -
            matrix[2][1] * matrix[0][2] * matrix[1][3] +
            matrix[2][1] * matrix[0][3] * matrix[1][2];

    inverse[1][3] = matrix[0][0] * matrix[1][2] * matrix[2][3] -
            matrix[0][0] * matrix[1][3] * matrix[2][2] -
            matrix[1][0] * matrix[0][2] * matrix[2][3] +
            matrix[1][0] * matrix[0][3] * matrix[2][2] +
            matrix[2][0] * matrix[0][2] * matrix[1][3] -
            matrix[2][0] * matrix[0][3] * matrix[1][2];

    inverse[2][3] = -matrix[0][0] * matrix[1][1] * matrix[2][3] +
            matrix[0][0] * matrix[1][3] * matrix[2][1] +
            matrix[1][0] * matrix[0][1] * matrix[2][3] -
            matrix[1][0] * matrix[0][3] * matrix[2][1] -
            matrix[2][0] * matrix[0][1] * matrix[1][3] +
            matrix[2][0] * matrix[0][3] * matrix[1][1];

    inverse[3][3] = matrix[0][0] * matrix[1][1] * matrix[2][2] -
            matrix[0][0] * matrix[1][2] * matrix[2][1] -
            matrix[1][0] * matrix[0][1] * matrix[2][2] +
            matrix[1][0] * matrix[0][2] * matrix[2][1] +
            matrix[2][0] * matrix[0][1] * matrix[1][2] -
            matrix[2][0] * matrix[0][2] * matrix[1][1];

    K det = matrix[0][0] * inverse[0][0] + matrix[0][1] * inverse[1][0] +
            matrix[0][2] * inverse[2][0] + matrix[0][3] * inverse[3][0];

    // return identity for singular or nearly singular matrices.
    if (std::abs(det) < 1e-40) {
        for (int i = 0; i < 4; ++i){
            inverse[i][i] = 1.0;
        }
        return 1.0;
    }
    K inv_det = 1.0 / det;
    inverse *= inv_det;

    return det;
}
} // end FMatrixHelp

namespace ISTLUtility {

//! invert matrix by calling FMatrixHelp::invert
template <typename K>
static inline void invertMatrix(FieldMatrix<K,1,1>& matrix)
{
    FieldMatrix<K,1,1> A ( matrix );
    FMatrixHelp::invertMatrix(A, matrix );
}

//! invert matrix by calling FMatrixHelp::invert
template <typename K>
static inline void invertMatrix(FieldMatrix<K,2,2>& matrix)
{
    FieldMatrix<K,2,2> A ( matrix );
    FMatrixHelp::invertMatrix(A, matrix );
}

//! invert matrix by calling FMatrixHelp::invert
template <typename K>
static inline void invertMatrix(FieldMatrix<K,3,3>& matrix)
{
    FieldMatrix<K,3,3> A ( matrix );
    FMatrixHelp::invertMatrix(A, matrix );
}

//! invert matrix by calling FMatrixHelp::invert
template <typename K>
static inline void invertMatrix(FieldMatrix<K,4,4>& matrix)
{
    FieldMatrix<K,4,4> A ( matrix );
    FMatrixHelp::invertMatrix(A, matrix );
}

//! invert matrix by calling matrix.invert
template <typename K, int n>
static inline void invertMatrix(FieldMatrix<K,n,n>& matrix)
{
#if ! DUNE_VERSION_NEWER( DUNE_COMMON, 2, 7 )
    Dune::FMatrixPrecision<K>::set_singular_limit(1.e-20);
#endif
    matrix.invert();
}

//! invert matrix by calling matrix.invert
template <typename K>
static inline void invertMatrix(Dune::DynamicMatrix<K>& matrix)
{
    // for 4 X 4 matrix, using the invertMatrix() function above
    // it is for temporary usage, mainly to reduce the huge burden of testing
    // what algorithm should be used to invert 4 X 4 matrix will be handled
    // as a seperate issue
    if (matrix.rows() == 4) {
        Dune::DynamicMatrix<K> A = matrix;
        FMatrixHelp::invertMatrix(A, matrix);
        return;
    }

#if ! DUNE_VERSION_NEWER( DUNE_COMMON, 2, 7 )
    Dune::FMatrixPrecision<K>::set_singular_limit(1.e-30);
#endif
    matrix.invert();

}

} // end ISTLUtility

template <class Scalar, int n, int m>
class MatrixBlock : public Dune::FieldMatrix<Scalar, n, m>
{
public:
    typedef Dune::FieldMatrix<Scalar, n, m>  BaseType;

    using BaseType :: operator= ;
    using BaseType :: rows;
    using BaseType :: cols;
    explicit MatrixBlock( const Scalar scalar = 0 ) : BaseType( scalar ) {}
    void invert()
    {
        ISTLUtility::invertMatrix( *this );
    }
    const BaseType& asBase() const { return static_cast< const BaseType& > (*this); }
    BaseType& asBase() { return static_cast< BaseType& > (*this); }
};

template<class K, int n, int m>
void
print_row(std::ostream& s, const MatrixBlock<K,n,m>& A,
          typename FieldMatrix<K,n,m>::size_type I,
          typename FieldMatrix<K,n,m>::size_type J,
          typename FieldMatrix<K,n,m>::size_type therow, int width,
          int precision)
{
    print_row(s, A.asBase(), I, J, therow, width, precision);
}

template<class K, int n, int m>
K& firstmatrixelement(MatrixBlock<K,n,m>& A)
{
   return firstmatrixelement( A.asBase() );
}



template<typename Scalar, int n, int m>
struct MatrixDimension< MatrixBlock< Scalar, n, m > >
: public MatrixDimension< typename MatrixBlock< Scalar, n, m >::BaseType >
{
};


#if HAVE_UMFPACK

/// \brief UMFPack specialization for MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template<typename T, typename A, int n, int m>
class UMFPack<BCRSMatrix<MatrixBlock<T,n,m>, A> >
    : public UMFPack<BCRSMatrix<FieldMatrix<T,n,m>, A> >
{
    typedef UMFPack<BCRSMatrix<FieldMatrix<T,n,m>, A> > Base;
    typedef BCRSMatrix<FieldMatrix<T,n,m>, A> Matrix;

public:
    typedef BCRSMatrix<MatrixBlock<T,n,m>, A> RealMatrix;

    UMFPack(const RealMatrix& matrix, int verbose, bool)
        : Base(reinterpret_cast<const Matrix&>(matrix), verbose)
    {}
};
#endif

#if HAVE_SUPERLU

/// \brief SuperLU specialization for MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template<typename T, typename A, int n, int m>
class SuperLU<BCRSMatrix<MatrixBlock<T,n,m>, A> >
    : public SuperLU<BCRSMatrix<FieldMatrix<T,n,m>, A> >
{
    typedef SuperLU<BCRSMatrix<FieldMatrix<T,n,m>, A> > Base;
    typedef BCRSMatrix<FieldMatrix<T,n,m>, A> Matrix;

public:
    typedef BCRSMatrix<MatrixBlock<T,n,m>, A> RealMatrix;

    SuperLU(const RealMatrix& matrix, int verbose, bool reuse=true)
        : Base(reinterpret_cast<const Matrix&>(matrix), verbose, reuse)
    {}
};
#endif


} // end namespace Dune

namespace Opm
{
namespace Detail
{
    //! calculates ret = A * B
    template< class TA, class TB, class TC, class PositiveSign >
    static inline void multMatrixImpl( const TA &A, // n x m
                                       const TB &B, // n x p
                                       TC &ret,     // m x p
                                       const PositiveSign )
    {
        typedef typename TA :: size_type size_type;
        typedef typename TA :: field_type K;
        assert( A.N() == B.N() );
        assert( A.M() == ret.N() );
        assert( B.M() == ret.M() );

        const size_type n = A.N();
        const size_type m = ret.N();
        const size_type p = B.M();
        for( size_type i = 0; i < m; ++i )
        {
            for( size_type j = 0; j < p; ++j )
            {
                K sum = 0;
                for( size_type k = 0; k < n; ++k )
                {
                    sum += A[ i ][ k ] * B[ k ][ j ];
                }
                // set value depending on given sign
                ret[ i ][ j ] = PositiveSign::value ? sum : -sum;
            }
        }
    }

    //! calculates ret = sign * (A^T * B)
    //! TA, TB, and TC are not necessarily FieldMatrix, but those should
    //! follow the Dune::DenseMatrix interface.
    template< class TA, class TB, class TC, class PositiveSign >
    static inline void multMatrixTransposedImpl ( const TA &A, // n x m
                                                  const TB &B, // n x p
                                                  TC &ret,     // m x p
                                                  const PositiveSign )
    {
        typedef typename TA :: size_type size_type;
        typedef typename TA :: field_type K;
        assert( A.N() == B.N() );
        assert( A.M() == ret.N() );
        assert( B.M() == ret.M() );

        const size_type n = A.N();
        const size_type m = ret.N();
        const size_type p = B.M();
        for( size_type i = 0; i < m; ++i )
        {
            for( size_type j = 0; j < p; ++j )
            {
                K sum = 0;
                for( size_type k = 0; k < n; ++k )
                {
                    sum += A[ k ][ i ] * B[ k ][ j ];
                }
                // set value depending on given sign
                ret[ i ][ j ] = PositiveSign::value ? sum : -sum;
            }
        }
    }

    //! calculates ret = A^T * B
    template <class DenseMatrixA, class DenseMatrixB, class DenseMatrixC>
    static inline void multMatrixTransposed(const DenseMatrixA& A,
                                            const DenseMatrixB& B,
                                            DenseMatrixC& ret)
    {
        multMatrixTransposedImpl( A, B, ret, std::true_type() );
    }

    //! calculates ret = -A^T * B
    template <class DenseMatrixA, class DenseMatrixB, class DenseMatrixC>
    static inline void negativeMultMatrixTransposed(const DenseMatrixA& A,
                                                    const DenseMatrixB& B,
                                                    DenseMatrixC& ret)
    {
        multMatrixTransposedImpl( A, B, ret, std::false_type() );
    }

    //! calculates ret = A * B
    template< class K>
    static inline void multMatrix(const Dune::DynamicMatrix<K>& A,
                                  const Dune::DynamicMatrix<K>& B,
                                  Dune::DynamicMatrix<K>& ret )
    {
        typedef typename Dune::DynamicMatrix<K> :: size_type size_type;

        const size_type m = A.rows();
        const size_type n = A.cols();

        assert(n == B.rows() );

        const size_type p = B.cols();

        ret.resize(m, p);

        for( size_type i = 0; i < m; ++i )
        {
            for( size_type j = 0; j < p; ++j )
            {
                ret[ i ][ j ] = K( 0 );
                for( size_type k = 0; k < n; ++k )
                    ret[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
            }
        }
    }

    //! perform out of place matrix inversion on C-style arrays
    //! must have a specified block_size
    template <int block_size>
    struct Inverter
    {
        template <typename K>
        void operator()(const K *matrix [[maybe_unused]], K *inverse [[maybe_unused]])
        {
            throw std::logic_error("Not implemented");
        }
    };

    //! perform out of place matrix inversion on C-style arrays
    template <>
    struct Inverter<4>
    {
        template <typename K>
        void operator()(const K *matrix, K *inverse)
        {
            // based on Dune::FMatrixHelp::invertMatrix
            inverse[0] = matrix[5] * matrix[10] * matrix[15] -
                    matrix[5] * matrix[11] * matrix[14] -
                    matrix[9] * matrix[6] * matrix[15] +
                    matrix[9] * matrix[7] * matrix[14] +
                    matrix[13] * matrix[6] * matrix[11] -
                    matrix[13] * matrix[7] * matrix[10];

            inverse[4] = -matrix[4] * matrix[10] * matrix[15] +
                    matrix[4] * matrix[11] * matrix[14] +
                    matrix[8] * matrix[6] * matrix[15] -
                    matrix[8] * matrix[7] * matrix[14] -
                    matrix[12] * matrix[6] * matrix[11] +
                    matrix[12] * matrix[7] * matrix[10];

            inverse[8] = matrix[4] * matrix[9] * matrix[15] -
                    matrix[4] * matrix[11] * matrix[13] -
                    matrix[8] * matrix[5] * matrix[15] +
                    matrix[8] * matrix[7] * matrix[13] +
                    matrix[12] * matrix[5] * matrix[11] -
                    matrix[12] * matrix[7] * matrix[9];

            inverse[12] = -matrix[4] * matrix[9] * matrix[14] +
                    matrix[4] * matrix[10] * matrix[13] +
                    matrix[8] * matrix[5] * matrix[14] -
                    matrix[8] * matrix[6] * matrix[13] -
                    matrix[12] * matrix[5] * matrix[10] +
                    matrix[12] * matrix[6] * matrix[9];

            inverse[1]= -matrix[1]  * matrix[10] * matrix[15] +
                    matrix[1] * matrix[11] * matrix[14] +
                    matrix[9] * matrix[2] * matrix[15] -
                    matrix[9] * matrix[3] * matrix[14] -
                    matrix[13] * matrix[2] * matrix[11] +
                    matrix[13] * matrix[3] * matrix[10];

            inverse[5] = matrix[0] * matrix[10] * matrix[15] -
                    matrix[0] * matrix[11] * matrix[14] -
                    matrix[8] * matrix[2] * matrix[15] +
                    matrix[8] * matrix[3] * matrix[14] +
                    matrix[12] * matrix[2] * matrix[11] -
                    matrix[12] * matrix[3] * matrix[10];

            inverse[9] = -matrix[0] * matrix[9] * matrix[15] +
                    matrix[0] * matrix[11] * matrix[13] +
                    matrix[8] * matrix[1] * matrix[15] -
                    matrix[8] * matrix[3] * matrix[13] -
                    matrix[12] * matrix[1] * matrix[11] +
                    matrix[12] * matrix[3] * matrix[9];

            inverse[13] = matrix[0] * matrix[9] * matrix[14] -
                    matrix[0] * matrix[10] * matrix[13] -
                    matrix[8] * matrix[1] * matrix[14] +
                    matrix[8] * matrix[2] * matrix[13] +
                    matrix[12] * matrix[1] * matrix[10] -
                    matrix[12] * matrix[2] * matrix[9];

            inverse[2] = matrix[1] * matrix[6] * matrix[15] -
                    matrix[1] * matrix[7] * matrix[14] -
                    matrix[5] * matrix[2] * matrix[15] +
                    matrix[5] * matrix[3] * matrix[14] +
                    matrix[13] * matrix[2] * matrix[7] -
                    matrix[13] * matrix[3] * matrix[6];

            inverse[6] = -matrix[0]  * matrix[6] * matrix[15] +
                    matrix[0] * matrix[7] * matrix[14] +
                    matrix[4] * matrix[2] * matrix[15] -
                    matrix[4] * matrix[3] * matrix[14] -
                    matrix[12] * matrix[2] * matrix[7] +
                    matrix[12] * matrix[3] * matrix[6];

            inverse[10] = matrix[0] * matrix[5] * matrix[15] -
                    matrix[0] * matrix[7] * matrix[13] -
                    matrix[4] * matrix[1] * matrix[15] +
                    matrix[4] * matrix[3] * matrix[13] +
                    matrix[12] * matrix[1] * matrix[7] -
                    matrix[12] * matrix[3] * matrix[5];

            inverse[14] = -matrix[0] * matrix[5] * matrix[14] +
                    matrix[0] * matrix[6] * matrix[13] +
                    matrix[4] * matrix[1] * matrix[14] -
                    matrix[4] * matrix[2] * matrix[13] -
                    matrix[12] * matrix[1] * matrix[6] +
                    matrix[12] * matrix[2] * matrix[5];

            inverse[3] = -matrix[1] * matrix[6] * matrix[11] +
                    matrix[1] * matrix[7] * matrix[10] +
                    matrix[5] * matrix[2] * matrix[11] -
                    matrix[5] * matrix[3] * matrix[10] -
                    matrix[9] * matrix[2] * matrix[7] +
                    matrix[9] * matrix[3] * matrix[6];

            inverse[7] = matrix[0] * matrix[6] * matrix[11] -
                    matrix[0] * matrix[7] * matrix[10] -
                    matrix[4] * matrix[2] * matrix[11] +
                    matrix[4] * matrix[3] * matrix[10] +
                    matrix[8] * matrix[2] * matrix[7] -
                    matrix[8] * matrix[3] * matrix[6];

            inverse[11] = -matrix[0] * matrix[5] * matrix[11] +
                    matrix[0] * matrix[7] * matrix[9] +
                    matrix[4] * matrix[1] * matrix[11] -
                    matrix[4] * matrix[3] * matrix[9] -
                    matrix[8] * matrix[1] * matrix[7] +
                    matrix[8] * matrix[3] * matrix[5];

            inverse[15] = matrix[0] * matrix[5] * matrix[10] -
                    matrix[0] * matrix[6] * matrix[9] -
                    matrix[4] * matrix[1] * matrix[10] +
                    matrix[4] * matrix[2] * matrix[9] +
                    matrix[8] * matrix[1] * matrix[6] -
                    matrix[8] * matrix[2] * matrix[5];

            K det = matrix[0] * inverse[0] + matrix[1] * inverse[4] +
                    matrix[2] * inverse[8] + matrix[3] * inverse[12];

            // return identity for singular or nearly singular matrices.
            if (std::abs(det) < 1e-40) {
                for (int i = 0; i < 4; ++i){
                    inverse[4*i + i] = 1.0;
                }
            }
            K inv_det = 1.0 / det;

            for (unsigned int i = 0; i < 4 * 4; ++i) {
                inverse[i] *= inv_det;
            }
        }
    };

    //! perform out of place matrix inversion on C-style arrays
    template <>
    struct Inverter<3>
    {
        template <typename K>
        void operator()(const K *matrix, K *inverse)
        {
            // code generated by maple, copied from Dune::DenseMatrix
            K t4  = matrix[0] * matrix[4];
            K t6  = matrix[0] * matrix[5];
            K t8  = matrix[1] * matrix[3];
            K t10 = matrix[2] * matrix[3];
            K t12 = matrix[1] * matrix[6];
            K t14 = matrix[2] * matrix[6];

            K det = (t4 * matrix[8] - t6 * matrix[7] - t8 * matrix[8] +
                          t10 * matrix[7] + t12 * matrix[5] - t14 * matrix[4]);
            K t17 = 1.0 / det;

            inverse[0] =  (matrix[4] * matrix[8] - matrix[5] * matrix[7]) * t17;
            inverse[1] = -(matrix[1] * matrix[8] - matrix[2] * matrix[7]) * t17;
            inverse[2] =  (matrix[1] * matrix[5] - matrix[2] * matrix[4]) * t17;
            inverse[3] = -(matrix[3] * matrix[8] - matrix[5] * matrix[6]) * t17;
            inverse[4] =  (matrix[0] * matrix[8] - t14) * t17;
            inverse[5] = -(t6 - t10) * t17;
            inverse[6] =  (matrix[3] * matrix[7] - matrix[4] * matrix[6]) * t17;
            inverse[7] = -(matrix[0] * matrix[7] - t12) * t17;
            inverse[8] =  (t4 - t8) * t17;
        }
    };

    //! perform out of place matrix inversion on C-style arrays
    template <>
    struct Inverter<2>
    {
        template <typename K>
        void operator()(const K *matrix, K *inverse)
        {
            // code based on Dune::DenseMatrix
            K detinv = matrix[0] * matrix[3] - matrix[1] * matrix[2];
            detinv = 1 / detinv;
            inverse[0] =  matrix[3] * detinv;
            inverse[1] = -matrix[1] * detinv;
            inverse[2] = -matrix[2] * detinv;
            inverse[3] =  matrix[0] * detinv;
        }
    };

    //! perform out of place matrix inversion on C-style arrays
    template <>
    struct Inverter<1>
    {
        template <typename K>
        void operator()(const K *matrix, K *inverse)
        {
            inverse[0] = 1.0 / matrix[0];
        }
    };

} // namespace Detail
} // namespace Opm

#endif
