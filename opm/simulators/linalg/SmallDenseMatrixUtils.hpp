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

#ifndef OPM_SMALL_DENSE_MATRIX_UTILS_HEADER_INCLUDED
#define OPM_SMALL_DENSE_MATRIX_UTILS_HEADER_INCLUDED

#include <dune/common/dynmatrix.hh>

namespace Opm
{
namespace detail
{
    //! calculates ret = A * B
    template< class TA, class TB, class TC, class PositiveSign >
    static inline void multMatrixImpl( const TA &A, // n x m
                                       const TB &B, // n x p
                                       TC &ret,     // m x p
                                       const PositiveSign )
    {
        using size_type = typename TA :: size_type;
        using K = typename TA :: field_type;
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
        using size_type = typename TA :: size_type;
        using K = typename TA :: field_type;
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
        using size_type = typename Dune::DynamicMatrix<K> :: size_type;

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

} // namespace detail
} // namespace Opm

#endif
