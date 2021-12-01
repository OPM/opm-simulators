/*
  Copyright 2019 Equinor ASA

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

#ifndef BILU0_HPP
#define BILU0_HPP

#include <mutex>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/ILUReorder.hpp>

#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/openclKernels.hpp>
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>
#include <opm/simulators/linalg/bda/ChowPatelIlu.hpp>


namespace Opm
{
namespace Accelerator
{

/// This class implements a Blocked ILU0 preconditioner
/// The decomposition is done on CPU, and reorders the rows of the matrix
template <unsigned int block_size>
class BILU0 : public Preconditioner<block_size>
{
    typedef Preconditioner<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;

private:
    std::unique_ptr<BlockedMatrix> LUmat = nullptr;
    std::shared_ptr<BlockedMatrix> rmat = nullptr; // only used with PAR_SIM
#if CHOW_PATEL
    std::unique_ptr<BlockedMatrix> Lmat = nullptr, Umat = nullptr;
#endif
    std::vector<double> invDiagVals;
    std::vector<int> diagIndex;
    std::vector<int> rowsPerColor;  // color i contains rowsPerColor[i] rows, which are processed in parallel
    std::vector<int> rowsPerColorPrefix;  // the prefix sum of rowsPerColor
    std::vector<int> toOrder, fromOrder;
    int numColors;
    std::once_flag pattern_uploaded;

    ILUReorder opencl_ilu_reorder;

    typedef struct {
        cl::Buffer invDiagVals;
        cl::Buffer diagIndex;
        cl::Buffer rowsPerColor;
#if CHOW_PATEL
        cl::Buffer Lvals, Lcols, Lrows;
        cl::Buffer Uvals, Ucols, Urows;
#else
        cl::Buffer LUvals, LUcols, LUrows;
#endif
    } GPU_storage;

    GPU_storage s;
    cl::Context *context;
    cl::CommandQueue *queue;
    std::vector<cl::Event> events;
    cl_int err;

#if CHOW_PATEL
    ChowPatelIlu<block_size> chowPatelIlu;
#endif

public:

    BILU0(ILUReorder opencl_ilu_reorder, int verbosity);

    void init(int Nb, int nnzb, std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue) override;

    // analysis, find reordering if specified
    bool analyze_matrix(BlockedMatrix *mat) override;

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix *mat) override;

    // apply preconditioner, x = prec(y)
    void apply(const cl::Buffer& y, cl::Buffer& x) override;

    int* getToOrder() override
    {
        return toOrder.data();
    }

    int* getFromOrder() override
    {
        return fromOrder.data();
    }

    BlockedMatrix* getRMat() override
    {
        return rmat.get();
    }

};

} // namespace Accelerator
} // namespace Opm

#endif

