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
#ifndef OPM_CPUSYSTEMBACKEND_HEADER_INCLUDED
#define OPM_CPUSYSTEMBACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include <dune/common/indices.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <cstddef>
#include <type_traits>

namespace Opm
{

// Backend traits for the CPU instantiation of SystemPreconditioner.
//
// Collects the types and the few operations that differ between the CPU and
// GPU paths, so that a single SystemPreconditioner serves both.  The
// GPU counterpart is gpuistl::GpuSystemBackend.
template <typename ScalarT, class ResOp, class ResCommT = Dune::Amg::SequentialInformation>
struct CpuSystemBackend
{
    static constexpr bool isGpu = false;
    static constexpr bool isParallel = !std::is_same_v<ResCommT, Dune::Amg::SequentialInformation>;
    static constexpr bool hasPerfectUpdate = true;

    using Scalar = ScalarT;
    using SysMatrix = SystemMatrix<Scalar>;
    using SysVector = SystemVector<Scalar>;
    using ResWork = ResVector<Scalar>;
    using WellWork = WellVector<Scalar>;
    using ResOperator = ResOp;
    using ResComm = ResCommT;
    using ResSolver = Dune::FlexibleSolver<ResOp>;

    static ResWork& res(SysVector& v) { return v[Dune::Indices::_0]; }
    static const ResWork& res(const SysVector& v) { return v[Dune::Indices::_0]; }
    static WellWork& well(SysVector& v) { return v[Dune::Indices::_1]; }
    static const WellWork& well(const SysVector& v) { return v[Dune::Indices::_1]; }

    // Work vectors are block vectors; sizes are given in blocks.
    template <class Vec>
    static void resize(Vec& v, std::size_t numBlocks, int /*blockSize*/) { v.resize(numBlocks); }

    template <class Vec>
    static void setZero(Vec& v) { v = 0.0; }

    // y += x
    template <class Vec>
    static void add(Vec& y, const Vec& x) { y += x; }

    // y -= M * x
    template <class Mat, class X, class Y>
    static void mmv(const Mat& M, const X& x, Y& y) { M.mmv(x, y); }

    // The well sub-solver always runs on the CPU, on this matrix.
    static const WWMatrix<Scalar>& cpuWellMatrix(const SysMatrix& S) { return *S.D; }

    static void updateResSolver(ResSolver& s) { s.preconditioner().update(); }
};

} // namespace Opm

#endif // OPM_CPUSYSTEMBACKEND_HEADER_INCLUDED
