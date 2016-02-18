// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2011 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file has been modified for use in the OPM project codebase.

#ifndef OPM_FASTSPARSEPRODUCT_HEADER_INCLUDED
#define OPM_FASTSPARSEPRODUCT_HEADER_INCLUDED

#include <Eigen/Sparse>

#include <algorithm>
#include <iterator>
#include <functional>
#include <limits>
#include <vector>

#include <Eigen/Core>

namespace Opm {

template < unsigned int depth >
struct QuickSort
{
  template <typename T>
  static inline void sort(T begin, T end)
  {
    if (begin != end)
    {
      T middle = std::partition (begin, end,
                                 std::bind2nd(std::less<typename std::iterator_traits<T>::value_type>(), *begin)
                                );
      QuickSort< depth-1 >::sort(begin, middle);

      // std::sort (max(begin + 1, middle), end);
      T new_middle = begin;
      QuickSort< depth-1 >::sort(++new_middle, end);
    }
  }
};

template <>
struct QuickSort< 0 >
{
  template <typename T>
  static inline void sort(T begin, T end)
  {
    // fall back to standard insertion sort
    std::sort( begin, end );
  }
};


template<typename Lhs, typename Rhs, typename ResultType>
void fastSparseProduct(const Lhs& lhs, const Rhs& rhs, ResultType& res)
{
  // initialize result
  res = ResultType(lhs.rows(), rhs.cols());

  // if one of the matrices does not contain non zero elements
  // the result will only contain an empty matrix
  if( lhs.nonZeros() == 0 || rhs.nonZeros() == 0 )
    return;

  typedef typename Eigen::internal::remove_all<Lhs>::type::Scalar Scalar;
  typedef typename Eigen::internal::remove_all<Lhs>::type::Index Index;

  // make sure to call innerSize/outerSize since we fake the storage order.
  Index rows = lhs.innerSize();
  Index cols = rhs.outerSize();
  eigen_assert(lhs.outerSize() == rhs.innerSize());

  std::vector<bool> mask(rows,false);
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> values(rows);
  Eigen::Matrix<Index, Eigen::Dynamic,1> indices(rows);

  // estimate the number of non zero entries
  // given a rhs column containing Y non zeros, we assume that the respective Y columns
  // of the lhs differs in average of one non zeros, thus the number of non zeros for
  // the product of a rhs column with the lhs is X+Y where X is the average number of non zero
  // per column of the lhs.
  // Therefore, we have nnz(lhs*rhs) = nnz(lhs) + nnz(rhs)
  Index estimated_nnz_prod = lhs.nonZeros() + rhs.nonZeros();

  res.setZero();
  res.reserve(Index(estimated_nnz_prod));

  //const Scalar epsilon = std::numeric_limits< Scalar >::epsilon();
  const Scalar epsilon = 0.0;

  // we compute each column of the result, one after the other
  for (Index j=0; j<cols; ++j)
  {
    Index nnz = 0;
    for (typename Rhs::InnerIterator rhsIt(rhs, j); rhsIt; ++rhsIt)
    {
      const Scalar y = rhsIt.value();
      for (typename Lhs::InnerIterator lhsIt(lhs, rhsIt.index()); lhsIt; ++lhsIt)
      {
        const Scalar val = lhsIt.value() * y;
        if( std::abs( val ) > epsilon )
        {
          const Index i = lhsIt.index();
          if(!mask[i])
          {
            mask[i] = true;
            values[i] = val;
            indices[nnz] = i;
            ++nnz;
          }
          else
            values[i] += val;
        }
      }
    }

    if( nnz > 1 )
    {
      // sort indices for sorted insertion to avoid later copying
      QuickSort< 1 >::sort( indices.data(), indices.data()+nnz );
    }

    res.startVec(j);
    // ordered insertion
    // still using insertBackByOuterInnerUnordered since we know what we are doing
    for(Index k=0; k<nnz; ++k)
    {
      const Index i = indices[k];
      res.insertBackByOuterInnerUnordered(j,i) = values[i];
      mask[i] = false;
    }

  }
  res.finalize();
}




inline void fastDiagSparseProduct(const std::vector<double>& lhs,
                                  const Eigen::SparseMatrix<double>& rhs,
                                  Eigen::SparseMatrix<double>& res)
{
    res = rhs;

    // Multiply rows by diagonal lhs.
    int n = res.cols();
    for (int col = 0; col < n; ++col) {
        typedef Eigen::SparseMatrix<double>::InnerIterator It;
        for (It it(res, col); it; ++it) {
            it.valueRef() *= lhs[it.row()];
        }
    }
}




inline void fastSparseDiagProduct(const Eigen::SparseMatrix<double>& lhs,
                                  const std::vector<double>& rhs,
                                  Eigen::SparseMatrix<double>& res)
{
    res = lhs;

    // Multiply columns by diagonal rhs.
    int n = res.cols();
    for (int col = 0; col < n; ++col) {
        typedef Eigen::SparseMatrix<double>::InnerIterator It;
        for (It it(res, col); it; ++it) {
            it.valueRef() *= rhs[col];
        }
    }
}

template<typename Lhs, typename Rhs>
inline bool
equalSparsityPattern(const Lhs& lhs, const Rhs& rhs)
{
    // if both matrices have equal storage and non zeros match, we can check sparsity pattern
    bool equal = (Lhs::IsRowMajor == Rhs::IsRowMajor) && (lhs.nonZeros() == rhs.nonZeros());

    // check complete sparsity pattern
    if( equal )
    {
        typedef typename Eigen::internal::remove_all<Lhs>::type::Index Index;
        const Index outerSize = lhs.outerSize();
        if( outerSize != rhs.outerSize() )
        {
            return false;
        }

        // outer indices
        const Index* rhsOuter = rhs.outerIndexPtr();
        const Index* lhsOuter = lhs.outerIndexPtr();
        for(Index i=0; i<=outerSize; ++i )
        {
            if( lhsOuter[ i ] != rhsOuter[ i ] ) {
                return false ;
            }
        }

        // inner indices
        const Index* rhsInner = rhs.innerIndexPtr();
        const Index* lhsInner = lhs.innerIndexPtr();

        const Index nnz = lhs.nonZeros();
        for( Index i=0; i<nnz; ++i)
        {
            if( lhsInner[ i ] != rhsInner[ i ] ) {
                return false;
            }
        }
    }

    return equal;
}

// this function substracts two sparse matrices
// if the sparsity pattern is the same a faster add/substract is performed
template<typename Lhs, typename Rhs>
inline void
fastSparseAdd(Lhs& lhs, const Rhs& rhs)
{
    if( equalSparsityPattern( lhs, rhs ) )
    {
        typedef typename Eigen::internal::remove_all<Lhs>::type::Scalar Scalar;
        typedef typename Eigen::internal::remove_all<Lhs>::type::Index Index;

        const Index nnz = lhs.nonZeros();

        // fast add using only the data pointers
        const Scalar* rhsV = rhs.valuePtr();
        Scalar* lhsV = lhs.valuePtr();

        for(Index i=0; i<nnz; ++i )
        {
            lhsV[ i ] += rhsV[ i ];
        }
    }
    else
    {
        // default Eigen operator+=
        lhs += rhs;
    }
}

// this function substracts two sparse matrices
// if the sparsity pattern is the same a faster add/substract is performed
template<typename Lhs, typename Rhs>
inline void
fastSparseSubstract(Lhs& lhs, const Rhs& rhs)
{
    if( equalSparsityPattern( lhs, rhs ) )
    {
        typedef typename Eigen::internal::remove_all<Lhs>::type::Scalar Scalar;
        typedef typename Eigen::internal::remove_all<Lhs>::type::Index Index;

        const Index nnz = lhs.nonZeros();

        // fast add using only the data pointers
        const Scalar* rhsV = rhs.valuePtr();
        Scalar* lhsV = lhs.valuePtr();

        for(Index i=0; i<nnz; ++i )
        {
            lhsV[ i ] -= rhsV[ i ];
        }
    }
    else
    {
        // default Eigen operator-=
        lhs -= rhs;
    }
}

} // end namespace Opm

#endif // OPM_FASTSPARSEPRODUCT_HEADER_INCLUDED
