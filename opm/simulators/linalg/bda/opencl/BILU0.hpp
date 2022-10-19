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

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>
#include <opm/simulators/linalg/bda/opencl/ChowPatelIlu.hpp>


namespace Opm
{
namespace Accelerator
{

/// This class implements a Blocked ILU0 preconditioner
/// The decomposition is done on GPU, using exact decomposition, or ChowPatel decomposition
/// The preconditioner is applied via two exact triangular solves
template <unsigned int block_size>
class BILU0 : public Preconditioner<block_size>
{
    typedef Preconditioner<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::context;
    using Base::queue;
    using Base::events;
    using Base::err;

private:
    std::unique_ptr<BlockedMatrix> LUmat = nullptr;
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

    bool opencl_ilu_parallel;

    typedef struct {
        cl::Buffer invDiagVals;    // nnz values of diagonal blocks of the matrix, inverted
        cl::Buffer diagIndex;      // index of diagonal block of each row, used to differentiate between lower and upper triangular part
        cl::Buffer rowsPerColor;   // number of rows for every color
        cl::Buffer rowIndices;     // mapping every row to another index
                                   // after mapping, all rows that are processed in parallel are contiguous
                                   // equal to the contents of fromOrder
#if CHOW_PATEL
        cl::Buffer Lvals, Lcols, Lrows;
        cl::Buffer Uvals, Ucols, Urows;
#else
        cl::Buffer LUvals, LUcols, LUrows;
#endif
    } GPU_storage;

    GPU_storage s;

#if CHOW_PATEL
    ChowPatelIlu<block_size> chowPatelIlu;
#endif

public:

    BILU0(bool opencl_ilu_parallel, int verbosity);

    // analysis, extract parallelism if specified
    bool analyze_matrix(BlockedMatrix *mat) override;
    bool analyze_matrix(BlockedMatrix *mat, BlockedMatrix *jacMat) override;

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix *mat) override;
    bool create_preconditioner(BlockedMatrix *mat, BlockedMatrix *jacMat) override;

    // apply preconditioner, x = prec(y)
    // via Lz = y
    // and Ux = z
    void apply(const cl::Buffer& y, cl::Buffer& x) override;

    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> get_preconditioner_structure()
    {
        return {{LUmat->rowPointers, LUmat->rowPointers + (Nb + 1)}, {LUmat->colIndices, LUmat->colIndices + nnzb}, diagIndex};
    }

    std::pair<cl::Buffer, cl::Buffer> get_preconditioner_data()
    {
#if CHOW_PATEL
        return std::make_pair(s.Lvals, s.invDiagVals); // send dummy, BISAI is disabled when ChowPatel is selected
#else
        return std::make_pair(s.LUvals, s.invDiagVals);
#endif
    }
};

} // namespace Accelerator
} // namespace Opm

#endif

