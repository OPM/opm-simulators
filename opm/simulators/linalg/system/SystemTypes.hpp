/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_SYSTEMTYPES_HEADER_INCLUDED
#define OPM_SYSTEMTYPES_HEADER_INCLUDED

#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

namespace Opm
{

// NOTE: These dimensions are hardcoded for standard 3-phase blackoil models
// (3 reservoir equations, 4 well equations). Models with a different number
// of conservation equations (e.g. EnablePolymerMW which adds an extra
// equation) are NOT supported by ISTLSolverSystem. A static_assert in
// ISTLSolverSystem guards against accidental misuse.
//
// To generalise, the types below (and the entire SystemPreconditioner /
// SystemPreconditionerFactory / WellMatrixMerger stack) would need to be
// templated on the dimension pair and the corresponding explicit
// instantiations updated.
inline constexpr int numResDofs = 3;
inline constexpr int numWellDofs = 4;

template<typename Scalar>
using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, numResDofs, numResDofs>>;
template<typename Scalar>
using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numResDofs, numWellDofs>>;
template<typename Scalar>
using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numWellDofs, numResDofs>>;
template<typename Scalar>
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numWellDofs, numWellDofs>>;

template<typename Scalar>
using ResVector = Dune::BlockVector<Dune::FieldVector<Scalar, numResDofs>>;
template<typename Scalar>
using WellVector = Dune::BlockVector<Dune::FieldVector<Scalar, numWellDofs>>;
template<typename Scalar>
using SystemVector = Dune::MultiTypeBlockVector<ResVector<Scalar>, WellVector<Scalar>>;

// --------------------------------------------------------------------------
// SystemMatrix: a lightweight read-only view over a 2×2 block-matrix
// structure.  All four sub-blocks are stored as const pointers; the actual
// data lives elsewhere (the reservoir block in ISTLSolver::matrix_, the
// well/coupling blocks in ISTLSolverSystem's merged-matrix members).
//
// Provides the operator interface required by Dune::MatrixAdapter and
// Dune::OverlappingSchwarzOperator (mv, usmv, N, M, field_type) as well as
// sub-block access via the index syntax  S[_0][_0]  used by
// SystemPreconditioner.
// --------------------------------------------------------------------------
template<typename Scalar> struct SystemMatrixRow0;  // forward
template<typename Scalar> struct SystemMatrixRow1;

template<typename Scalar>
class SystemMatrix
{
public:
    using size_type  = std::size_t;
    using field_type = Scalar;

    static constexpr size_type N() { return 2; }
    static constexpr size_type M() { return 2; }

    // Block pointers — set directly by the owning solver.
    const RRMatrix<Scalar>* A = nullptr;  // (0,0) reservoir
    const RWMatrix<Scalar>* C = nullptr;  // (0,1) reservoir–well coupling
    const WRMatrix<Scalar>* B = nullptr;  // (1,0) well–reservoir coupling
    const WWMatrix<Scalar>* D = nullptr;  // (1,1) well

    // Sub-block access: S[_0][_0], S[_0][_1], S[_1][_0], S[_1][_1]
    inline SystemMatrixRow0<Scalar> operator[](Dune::index_constant<0>) const;
    inline SystemMatrixRow1<Scalar> operator[](Dune::index_constant<1>) const;

    // Matrix-vector products required by Dune linear operators.
    // y = S * x  =  (A*x0 + C*x1;  B*x0 + D*x1)
    // Achieved by: y0 = A*x0 (mv), then y0 += C*x1 (umv); similarly for y1.
    void mv(const SystemVector<Scalar>& x, SystemVector<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->mv (x[_0], y[_0]);   C->umv(x[_1], y[_0]);
        B->mv (x[_0], y[_1]);   D->umv(x[_1], y[_1]);
    }

    // y += S * x  =  (y0 += A*x0 + C*x1;  y1 += B*x0 + D*x1)
    void umv(const SystemVector<Scalar>& x, SystemVector<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->umv(x[_0], y[_0]);   C->umv(x[_1], y[_0]);
        B->umv(x[_0], y[_1]);   D->umv(x[_1], y[_1]);
    }

    // y += alpha * S * x  =  (y0 += alpha*(A*x0 + C*x1);  y1 += alpha*(B*x0 + D*x1))
    void usmv(field_type alpha, const SystemVector<Scalar>& x, SystemVector<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->usmv(alpha, x[_0], y[_0]);   C->usmv(alpha, x[_1], y[_0]);
        B->usmv(alpha, x[_0], y[_1]);   D->usmv(alpha, x[_1], y[_1]);
    }
};

// Row proxies for  S[row][col]  — simple aggregates, no back-pointers.
template<typename Scalar>
struct SystemMatrixRow0
{
    const RRMatrix<Scalar>& A;
    const RWMatrix<Scalar>& C;
    const RRMatrix<Scalar>& operator[](Dune::index_constant<0>) const { return A; }
    const RWMatrix<Scalar>& operator[](Dune::index_constant<1>) const { return C; }
};

template<typename Scalar>
struct SystemMatrixRow1
{
    const WRMatrix<Scalar>& B;
    const WWMatrix<Scalar>& D;
    const WRMatrix<Scalar>& operator[](Dune::index_constant<0>) const { return B; }
    const WWMatrix<Scalar>& operator[](Dune::index_constant<1>) const { return D; }
};

template<typename Scalar>
SystemMatrixRow0<Scalar> SystemMatrix<Scalar>::operator[](Dune::index_constant<0>) const
{ return {*A, *C}; }

template<typename Scalar>
SystemMatrixRow1<Scalar> SystemMatrix<Scalar>::operator[](Dune::index_constant<1>) const
{ return {*B, *D}; }

} // namespace Opm

#endif // OPM_SYSTEMTYPES_HEADER_INCLUDED