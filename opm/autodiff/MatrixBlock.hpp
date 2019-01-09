/*
  Copyright 2016 IRIS AS

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
    //! calculates ret = A^T * B
    template< class K, int m, int n, int p >
    static inline void multMatrixTransposed(const Dune::FieldMatrix< K, n, m >& A,
                                            const Dune::FieldMatrix< K, n, p >& B,
                                            Dune::FieldMatrix< K, m, p >& ret)
    {
        typedef typename Dune::FieldMatrix< K, m, p > :: size_type size_type;

        for( size_type i = 0; i < m; ++i )
        {
            for( size_type j = 0; j < p; ++j )
            {
                ret[ i ][ j ] = K( 0 );
                for( size_type k = 0; k < n; ++k )
                    ret[ i ][ j ] += A[ k ][ i ] * B[ k ][ j ];
            }
        }
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


    //! calculates ret = A^T * B
    template< class K, int m, int p >
    static inline void multMatrixTransposed(const Dune::DynamicMatrix<K>& A,
                                            const Dune::DynamicMatrix<K>& B,
                                            Dune::FieldMatrix< K, m, p>& ret )
    {
        typedef typename Dune::DynamicMatrix<K> :: size_type size_type;

        // A is a tranpose matrix
        const size_type n = A.rows();
        assert(m == A.cols() );

        assert(n == B.rows() );
        assert(p == B.cols() );

        for( size_type i = 0; i < m; ++i )
        {
            for( size_type j = 0; j < p; ++j )
            {
                ret[ i ][ j ] = K( 0 );
                for( size_type k = 0; k < n; ++k )
                    ret[ i ][ j ] += A[ k ][ i ] * B[ k ][ j ];
            }
        }
    }
} // namespace Detail
} // namespace Opm

#endif
