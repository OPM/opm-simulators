#pragma once

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <opm/simulators/linalg/matrixblock.hh>

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
using RRMatrixT = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, numResDofs, numResDofs>>;
template<typename Scalar>
using RWMatrixT = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numResDofs, numWellDofs>>;
template<typename Scalar>
using WRMatrixT = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numWellDofs, numResDofs>>;
template<typename Scalar>
using WWMatrixT = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numWellDofs, numWellDofs>>;

template<typename Scalar>
using ResVectorT = Dune::BlockVector<Dune::FieldVector<Scalar, numResDofs>>;
template<typename Scalar>
using WellVectorT = Dune::BlockVector<Dune::FieldVector<Scalar, numWellDofs>>;
template<typename Scalar>
using SystemVectorT = Dune::MultiTypeBlockVector<ResVectorT<Scalar>, WellVectorT<Scalar>>;

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
template<typename Scalar> struct SystemMatrixRow0T;  // forward
template<typename Scalar> struct SystemMatrixRow1T;

template<typename Scalar>
class SystemMatrixT
{
public:
    using size_type  = std::size_t;
    using field_type = Scalar;

    static constexpr size_type N() { return 2; }
    static constexpr size_type M() { return 2; }

    // Block pointers — set directly by the owning solver.
    const RRMatrixT<Scalar>* A = nullptr;  // (0,0) reservoir
    const RWMatrixT<Scalar>* C = nullptr;  // (0,1) reservoir–well coupling
    const WRMatrixT<Scalar>* B = nullptr;  // (1,0) well–reservoir coupling
    const WWMatrixT<Scalar>* D = nullptr;  // (1,1) well

    // Sub-block access: S[_0][_0], S[_0][_1], S[_1][_0], S[_1][_1]
    inline SystemMatrixRow0T<Scalar> operator[](Dune::index_constant<0>) const;
    inline SystemMatrixRow1T<Scalar> operator[](Dune::index_constant<1>) const;

    // Matrix-vector products required by Dune linear operators.
    void mv(const SystemVectorT<Scalar>& x, SystemVectorT<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->mv (x[_0], y[_0]);   C->umv(x[_1], y[_0]);
        B->mv (x[_0], y[_1]);   D->umv(x[_1], y[_1]);
    }

    void umv(const SystemVectorT<Scalar>& x, SystemVectorT<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->umv(x[_0], y[_0]);   C->umv(x[_1], y[_0]);
        B->umv(x[_0], y[_1]);   D->umv(x[_1], y[_1]);
    }

    void usmv(field_type alpha, const SystemVectorT<Scalar>& x, SystemVectorT<Scalar>& y) const
    {
        using namespace Dune::Indices;
        A->usmv(alpha, x[_0], y[_0]);   C->usmv(alpha, x[_1], y[_0]);
        B->usmv(alpha, x[_0], y[_1]);   D->usmv(alpha, x[_1], y[_1]);
    }
};

// Row proxies for  S[row][col]  — simple aggregates, no back-pointers.
template<typename Scalar>
struct SystemMatrixRow0T
{
    const RRMatrixT<Scalar>& A;
    const RWMatrixT<Scalar>& C;
    const RRMatrixT<Scalar>& operator[](Dune::index_constant<0>) const { return A; }
    const RWMatrixT<Scalar>& operator[](Dune::index_constant<1>) const { return C; }
};

template<typename Scalar>
struct SystemMatrixRow1T
{
    const WRMatrixT<Scalar>& B;
    const WWMatrixT<Scalar>& D;
    const WRMatrixT<Scalar>& operator[](Dune::index_constant<0>) const { return B; }
    const WWMatrixT<Scalar>& operator[](Dune::index_constant<1>) const { return D; }
};

template<typename Scalar>
SystemMatrixRow0T<Scalar> SystemMatrixT<Scalar>::operator[](Dune::index_constant<0>) const
{ return {*A, *C}; }

template<typename Scalar>
SystemMatrixRow1T<Scalar> SystemMatrixT<Scalar>::operator[](Dune::index_constant<1>) const
{ return {*B, *D}; }

} // namespace Opm
