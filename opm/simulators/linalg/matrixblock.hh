// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef EWOMS_MATRIX_BLOCK_HH
#define EWOMS_MATRIX_BLOCK_HH

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <dune/common/fmatrix.hh>

namespace Opm {
namespace MatrixBlockHelp {

template <typename K, int m, int n>
static inline void invertMatrix(Dune::FieldMatrix<K, m, n>& matrix)
{ matrix.invert(); }

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 1, 1>& matrix)
{
    matrix[0][0] = 1.0/matrix[0][0];
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 2, 2>& matrix)
{
    Dune::FieldMatrix<K, 2, 2> tmp(matrix);
    Dune::FMatrixHelp::invertMatrix(tmp, matrix);
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 3, 3>& matrix)
{
    Dune::FieldMatrix<K, 3, 3> tmp(matrix);
    Dune::FMatrixHelp::invertMatrix(tmp, matrix);
}

template <typename K>
static inline void invertMatrix(Dune::FieldMatrix<K, 4, 4>& matrix)
{
    Dune::FieldMatrix<K, 4, 4> tmp(matrix);

    matrix[0][0] =
        tmp[1][1]*tmp[2][2]*tmp[3][3] -
        tmp[1][1]*tmp[2][3]*tmp[3][2] -
        tmp[2][1]*tmp[1][2]*tmp[3][3] +
        tmp[2][1]*tmp[1][3]*tmp[3][2] +
        tmp[3][1]*tmp[1][2]*tmp[2][3] -
        tmp[3][1]*tmp[1][3]*tmp[2][2];

    matrix[1][0] =
       -tmp[1][0]*tmp[2][2]*tmp[3][3] +
        tmp[1][0]*tmp[2][3]*tmp[3][2] +
        tmp[2][0]*tmp[1][2]*tmp[3][3] -
        tmp[2][0]*tmp[1][3]*tmp[3][2] -
        tmp[3][0]*tmp[1][2]*tmp[2][3] +
        tmp[3][0]*tmp[1][3]*tmp[2][2];

    matrix[2][0] =
        tmp[1][0]*tmp[2][1]*tmp[3][3] -
        tmp[1][0]*tmp[2][3]*tmp[3][1] -
        tmp[2][0]*tmp[1][1]*tmp[3][3] +
        tmp[2][0]*tmp[1][3]*tmp[3][1] +
        tmp[3][0]*tmp[1][1]*tmp[2][3] -
        tmp[3][0]*tmp[1][3]*tmp[2][1];

    matrix[3][0] =
       -tmp[1][0]*tmp[2][1]*tmp[3][2] +
        tmp[1][0]*tmp[2][2]*tmp[3][1] +
        tmp[2][0]*tmp[1][1]*tmp[3][2] -
        tmp[2][0]*tmp[1][2]*tmp[3][1] -
        tmp[3][0]*tmp[1][1]*tmp[2][2] +
        tmp[3][0]*tmp[1][2]*tmp[2][1];

    matrix[0][1] =
       -tmp[0][1]*tmp[2][2]*tmp[3][3] +
        tmp[0][1]*tmp[2][3]*tmp[3][2] +
        tmp[2][1]*tmp[0][2]*tmp[3][3] -
        tmp[2][1]*tmp[0][3]*tmp[3][2] -
        tmp[3][1]*tmp[0][2]*tmp[2][3] +
        tmp[3][1]*tmp[0][3]*tmp[2][2];

    matrix[1][1] =
        tmp[0][0]*tmp[2][2]*tmp[3][3] -
        tmp[0][0]*tmp[2][3]*tmp[3][2] -
        tmp[2][0]*tmp[0][2]*tmp[3][3] +
        tmp[2][0]*tmp[0][3]*tmp[3][2] +
        tmp[3][0]*tmp[0][2]*tmp[2][3] -
        tmp[3][0]*tmp[0][3]*tmp[2][2];

    matrix[2][1] =
       -tmp[0][0]*tmp[2][1]*tmp[3][3] +
        tmp[0][0]*tmp[2][3]*tmp[3][1] +
        tmp[2][0]*tmp[0][1]*tmp[3][3] -
        tmp[2][0]*tmp[0][3]*tmp[3][1] -
        tmp[3][0]*tmp[0][1]*tmp[2][3] +
        tmp[3][0]*tmp[0][3]*tmp[2][1];

    matrix[3][1] =
        tmp[0][0]*tmp[2][1]*tmp[3][2] -
        tmp[0][0]*tmp[2][2]*tmp[3][1] -
        tmp[2][0]*tmp[0][1]*tmp[3][2] +
        tmp[2][0]*tmp[0][2]*tmp[3][1] +
        tmp[3][0]*tmp[0][1]*tmp[2][2] -
        tmp[3][0]*tmp[0][2]*tmp[2][1];

    matrix[0][2] =
        tmp[0][1]*tmp[1][2]*tmp[3][3] -
        tmp[0][1]*tmp[1][3]*tmp[3][2] -
        tmp[1][1]*tmp[0][2]*tmp[3][3] +
        tmp[1][1]*tmp[0][3]*tmp[3][2] +
        tmp[3][1]*tmp[0][2]*tmp[1][3] -
        tmp[3][1]*tmp[0][3]*tmp[1][2];

    matrix[1][2] =
       -tmp[0][0] *tmp[1][2]*tmp[3][3] +
        tmp[0][0]*tmp[1][3]*tmp[3][2] +
        tmp[1][0]*tmp[0][2]*tmp[3][3] -
        tmp[1][0]*tmp[0][3]*tmp[3][2] -
        tmp[3][0]*tmp[0][2]*tmp[1][3] +
        tmp[3][0]*tmp[0][3]*tmp[1][2];

    matrix[2][2] =
        tmp[0][0]*tmp[1][1]*tmp[3][3] -
        tmp[0][0]*tmp[1][3]*tmp[3][1] -
        tmp[1][0]*tmp[0][1]*tmp[3][3] +
        tmp[1][0]*tmp[0][3]*tmp[3][1] +
        tmp[3][0]*tmp[0][1]*tmp[1][3] -
        tmp[3][0]*tmp[0][3]*tmp[1][1];

    matrix[3][2] =
       -tmp[0][0]*tmp[1][1]*tmp[3][2] +
        tmp[0][0]*tmp[1][2]*tmp[3][1] +
        tmp[1][0]*tmp[0][1]*tmp[3][2] -
        tmp[1][0]*tmp[0][2]*tmp[3][1] -
        tmp[3][0]*tmp[0][1]*tmp[1][2] +
        tmp[3][0]*tmp[0][2]*tmp[1][1];

    matrix[0][3] =
       -tmp[0][1]*tmp[1][2]*tmp[2][3] +
        tmp[0][1]*tmp[1][3]*tmp[2][2] +
        tmp[1][1]*tmp[0][2]*tmp[2][3] -
        tmp[1][1]*tmp[0][3]*tmp[2][2] -
        tmp[2][1]*tmp[0][2]*tmp[1][3] +
        tmp[2][1]*tmp[0][3]*tmp[1][2];

    matrix[1][3] =
        tmp[0][0]*tmp[1][2]*tmp[2][3] -
        tmp[0][0]*tmp[1][3]*tmp[2][2] -
        tmp[1][0]*tmp[0][2]*tmp[2][3] +
        tmp[1][0]*tmp[0][3]*tmp[2][2] +
        tmp[2][0]*tmp[0][2]*tmp[1][3] -
        tmp[2][0]*tmp[0][3]*tmp[1][2];

    matrix[2][3] =
       -tmp[0][0]*tmp[1][1]*tmp[2][3] +
        tmp[0][0]*tmp[1][3]*tmp[2][1] +
        tmp[1][0]*tmp[0][1]*tmp[2][3] -
        tmp[1][0]*tmp[0][3]*tmp[2][1] -
        tmp[2][0]*tmp[0][1]*tmp[1][3] +
        tmp[2][0]*tmp[0][3]*tmp[1][1];

    matrix[3][3] =
        tmp[0][0]*tmp[1][1]*tmp[2][2] -
        tmp[0][0]*tmp[1][2]*tmp[2][1] -
        tmp[1][0]*tmp[0][1]*tmp[2][2] +
        tmp[1][0]*tmp[0][2]*tmp[2][1] +
        tmp[2][0]*tmp[0][1]*tmp[1][2] -
        tmp[2][0]*tmp[0][2]*tmp[1][1];

    K det =
        tmp[0][0]*matrix[0][0] +
        tmp[0][1]*matrix[1][0] +
        tmp[0][2]*matrix[2][0] +
        tmp[0][3]*matrix[3][0];

    // return identity for singular or nearly singular matrices.
    if (std::abs(det) < 1e-40)
        matrix = std::numeric_limits<K>::quiet_NaN();
    else
        matrix *= 1.0/det;
}
} // namespace MatrixBlockHelp

template <class Scalar, int n, int m>
class MatrixBlock : public Dune::FieldMatrix<Scalar, n, m>
{
public:
    using BaseType = Dune::FieldMatrix<Scalar, n, m> ;

    using BaseType::operator= ;
    using BaseType::rows;
    using BaseType::cols;

    MatrixBlock()
        : BaseType(Scalar(0.0))
    {}

    explicit MatrixBlock(const Scalar value)
        : BaseType(value)
    {}

    void invert()
    { Opm::MatrixBlockHelp::invertMatrix(asBase()); }

    const BaseType& asBase() const
    { return static_cast<const BaseType&>(*this); }

    BaseType& asBase()
    { return static_cast<BaseType&>(*this); }
};

} // namespace Opm

namespace Dune {

template<class K, int n, int m>
void print_row(std::ostream& s, const Opm::MatrixBlock<K, n, m>& A,
               typename FieldMatrix<K, n, m>::size_type I,
               typename FieldMatrix<K, n, m>::size_type J,
               typename FieldMatrix<K, n, m>::size_type therow,
               int width,
               int precision)
{ print_row(s, A.asBase(), I, J, therow, width, precision); }

template<class K, int n, int m>
K& firstmatrixelement(Opm::MatrixBlock<K, n, m>& A)
{ return firstmatrixelement(A.asBase()); }

template <typename Scalar, int n, int m>
struct MatrixDimension<Opm::MatrixBlock<Scalar, n, m> >
    : public MatrixDimension<typename Opm::MatrixBlock<Scalar, n, m>::BaseType>
{ };


#if HAVE_UMFPACK
/// \brief UMFPack specialization for Opm::MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template <typename T, typename A, int n, int m>
class UMFPack<BCRSMatrix<Opm::MatrixBlock<T, n, m>, A> >
    : public UMFPack<BCRSMatrix<FieldMatrix<T, n, m>, A> >
{
    using Base = UMFPack<BCRSMatrix<FieldMatrix<T, n, m>, A> >;
    using Matrix = BCRSMatrix<FieldMatrix<T, n, m>, A>;

public:
    using RealMatrix = BCRSMatrix<Opm::MatrixBlock<T, n, m>, A>;

    UMFPack(const RealMatrix& matrix, int verbose, bool)
        : Base(reinterpret_cast<const Matrix&>(matrix), verbose)
    {}
};
#endif

#if HAVE_SUPERLU
/// \brief SuperLU specialization for Opm::MatrixBlock to make AMG happy
///
/// Without this the empty default implementation would be used.
template <typename T, typename A, int n, int m>
class SuperLU<BCRSMatrix<Opm::MatrixBlock<T, n, m>, A> >
    : public SuperLU<BCRSMatrix<FieldMatrix<T, n, m>, A> >
{
    using Base = SuperLU<BCRSMatrix<FieldMatrix<T, n, m>, A> >;
    using Matrix = BCRSMatrix<FieldMatrix<T, n, m>, A>;

public:
    using RealMatrix = BCRSMatrix<Opm::MatrixBlock<T, n, m>, A>;

    SuperLU(const RealMatrix& matrix, int verb, bool reuse=true)
        : Base(reinterpret_cast<const Matrix&>(matrix), verb, reuse)
    {}
};
#endif

} // end namespace Dune


#endif
