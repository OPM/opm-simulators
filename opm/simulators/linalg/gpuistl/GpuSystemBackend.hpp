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
#ifndef OPM_GPUSYSTEMBACKEND_HEADER_INCLUDED
#define OPM_GPUSYSTEMBACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/gpuistl/GpuSystemTypes.hpp>
#include <opm/simulators/linalg/gpuistl/detail/FlexibleSolverWrapper.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>

#include <dune/common/parallel/communication.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <cstddef>

namespace Opm::gpuistl
{

// Backend traits for the GPU instantiation of SystemPreconditioner; see
// CpuSystemBackend for the CPU counterpart.
//
// The reservoir sub-solvers run on the GPU through FlexibleSolverWrapper.
// The well sub-solver uses the host copy of the well matrix and runs on the
// CPU.  Sequential runs only.
template <typename ScalarT>
struct GpuSystemBackend
{
    static constexpr bool isGpu = true;
    static constexpr bool isParallel = false;
    static constexpr bool hasPerfectUpdate = false;

    using Scalar = ScalarT;
    using SysMatrix = GpuSystemMatrix<Scalar>;
    using SysVector = GpuSystemVector<Scalar>;
    using ResWork = GpuVector<Scalar>;
    using WellWork = GpuVector<Scalar>;
    // FlexibleSolverWrapper creates its own operator, so none is held here.
    struct ResOperator {};
    using ResComm = Dune::Amg::SequentialInformation;
    using ResSolver = detail::FlexibleSolverWrapper<GpuRRMatrix<Scalar>,
                                                    GpuVector<Scalar>,
                                                    Dune::Communication<int>>;

    static ResWork& res(SysVector& v) { return v.res; }
    static const ResWork& res(const SysVector& v) { return v.res; }
    static WellWork& well(SysVector& v) { return v.well; }
    static const WellWork& well(const SysVector& v) { return v.well; }

    // Work vectors are flat; sizes are given in blocks.
    static void resize(ResWork& v, std::size_t numBlocks, int blockSize)
    { v.resize(numBlocks * blockSize); }

    template <class Vec>
    static void setZero(Vec& v)
    {
        // A run without open wells has empty well vectors.
        if (v.dim() > 0) {
            v = 0.0;
        }
    }

    // y += x
    template <class Vec>
    static void add(Vec& y, const Vec& x)
    {
        if (y.dim() > 0) {
            y += x;
        }
    }

    // y -= M * x
    template <class Mat>
    static void mmv(const Mat& M, const GpuVector<Scalar>& x, GpuVector<Scalar>& y)
    {
        if (x.dim() > 0 && y.dim() > 0) {
            M.usmv(Scalar(-1), x, y);
        }
    }

    static const WWMatrix<Scalar>& cpuWellMatrix(const SysMatrix& S) { return **S.hostD; }

    static void updateResSolver(ResSolver& s) { s.update(); }
};

// SystemPreconditioner on the GPU backend.  The reservoir operator type is a
// placeholder: FlexibleSolverWrapper builds its own.
template <typename Scalar>
using GpuSystemPreconditioner = SystemPreconditioner<Scalar,
                                                     typename GpuSystemBackend<Scalar>::ResOperator,
                                                     Dune::Amg::SequentialInformation,
                                                     GpuSystemBackend<Scalar>>;

} // namespace Opm::gpuistl

#endif // OPM_GPUSYSTEMBACKEND_HEADER_INCLUDED
