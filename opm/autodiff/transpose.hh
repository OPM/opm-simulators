#ifndef DUNE_MATRIX_VECTOR_TRANSPOSE_HH
#define DUNE_MATRIX_VECTOR_TRANSPOSE_HH
// HMN copyed from dune-matrix-vector
#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/scaledidmatrix.hh>
#include <ewoms/linear/matrixblock.hh>
//#include <dune/matrix-vector/axpy.hh>

namespace Dune {
  namespace MatrixVector {

    template <class A>
    struct TransposeHelper;

    template <class MatrixType>
    using Transposed = typename TransposeHelper<MatrixType>::TransposedType;

    template <class A>
    struct TransposeHelper {
      typedef A TransposedType;

      static void transpose(const A& a, TransposedType& aT) {
        DUNE_UNUSED_PARAMETER(a);
        DUNE_UNUSED_PARAMETER(aT);
        DUNE_THROW(Dune::Exception,
                   "Not implemented for general matrix types!");
      }
    };

    //! Specialization for Dune::FieldMatrix
    template <class T, int n, int m>
    struct TransposeHelper<Dune::FieldMatrix<T, n, m>> {
      typedef Dune::FieldMatrix<T, n, m> MatrixType;
      typedef Dune::FieldMatrix<T, m, n> TransposedType;

      static void transpose(const MatrixType& a, TransposedType& aT) {
        for (int row = 0; row < m; ++row)
          for (int col = 0; col < n; ++col)
            aT[row][col] = a[col][row];
      }
    };

    //! Specialization for Dune::FieldMatrix
    template <class T, int n, int m>
    struct TransposeHelper<Ewoms::MatrixBlock<T, n, m>> {
      typedef Ewoms::MatrixBlock<T, n, m> MatrixType;
      typedef Ewoms::MatrixBlock<T, m, n> TransposedType;

      static void transpose(const MatrixType& a, TransposedType& aT) {
        for (int row = 0; row < m; ++row)
          for (int col = 0; col < n; ++col)
            aT[row][col] = a[col][row];
      }
    };

    template <class T, int n>
    struct TransposeHelper<Dune::DiagonalMatrix<T, n>> {
      typedef Dune::DiagonalMatrix<T, n> MatrixType;
      typedef Dune::DiagonalMatrix<T, n> TransposedType;

      static void transpose(const MatrixType& a, TransposedType& aT) { aT = a; }
    };

    template <class T, int n>
    struct TransposeHelper<Dune::ScaledIdentityMatrix<T, n>> {
      typedef Dune::ScaledIdentityMatrix<T, n> MatrixType;
      typedef Dune::ScaledIdentityMatrix<T, n> TransposedType;

      static void transpose(const MatrixType& a, TransposedType& aT) { aT = a; }
    };

    //! Specialization for Dune::BCRSMatrix Type
    template <class A>
    struct TransposeHelper<Dune::BCRSMatrix<A>> {
      typedef Dune::BCRSMatrix<Transposed<A>> TransposedType;

      static void transpose(const Dune::BCRSMatrix<A>& a, TransposedType& aT) {
        Dune::MatrixIndexSet idxSetaT(a.M(), a.N());

        typedef typename Dune::BCRSMatrix<A>::ConstColIterator ColIterator;

        // add indices into transposed matrix
        for (size_t row = 0; row < a.N(); ++row) {
          ColIterator col = a[row].begin();
          ColIterator end = a[row].end();

          for (; col != end; ++col)
            idxSetaT.add(col.index(), row);
        }

        idxSetaT.exportIdx(aT);

        for (size_t row = 0; row < a.N(); ++row) {
          ColIterator col = a[row].begin();
          ColIterator end = a[row].end();

          for (; col != end; ++col)
            TransposeHelper<A>::transpose(*col, aT[col.index()][row]);
        }
      }
    };

    //! Compute the transposed of a matrix
    template <class A>
    void transpose(const A& a, Transposed<A>& aT) {
      TransposeHelper<A>::transpose(a, aT);
    }
  }
}

#endif




