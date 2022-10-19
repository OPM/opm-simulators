/*
  Copyright 2022 Equinor ASA

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

#ifndef BISAI_HPP
#define BISAI_HPP

#include <mutex>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

namespace Opm
{
namespace Accelerator
{

class BlockedMatrix;

/// This class implements a Blocked version of the Incomplete Sparse Approximate Inverse (ISAI) preconditioner.
/// Inspired by the paper "Incomplete Sparse Approximate Inverses for Parallel Preconditioning" by Anzt et. al.
template <unsigned int block_size>
class BISAI : public Preconditioner<block_size>
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
    std::once_flag initialize;

    std::vector<int> colPointers;
    std::vector<int> rowIndices;
    std::vector<int> diagIndex;
    std::vector<int> csrToCscOffsetMap;
    std::vector<double> invLvals;
    std::vector<double> invUvals;

    cl::Buffer d_colPointers;
    cl::Buffer d_rowIndices;
    cl::Buffer d_csrToCscOffsetMap;
    cl::Buffer d_diagIndex;
    cl::Buffer d_LUvals;
    cl::Buffer d_invDiagVals;
    cl::Buffer d_invLvals;
    cl::Buffer d_invUvals;
    cl::Buffer d_invL_x;

    bool opencl_ilu_parallel;
    std::unique_ptr<BILU0<block_size> > bilu0;

    /// Struct that holds the structure of the small subsystems for each column
    typedef struct{
        /// This vector holds the cumulative sum for the number of non-zero blocks for each subsystem.
        /// Works similarly to row and column pointers for the CSR and CSC matrix representations.
        std::vector<int> subsystemPointers;
        /// This vector holds the indices of the non-zero blocks for the target subsystem. These blocks are
        /// the ones that are present in the shadow set of the non-zero blocks of column j of the main matrix,
        /// as described in section 2.3 of the paper. The amount of non-zero blocks for j-th subsystem is
        /// given by subsystemPointers[j+1] - subsystemPointers[j].
        std::vector<int> nzIndices;
        /// This vector holds the indices of the already known values of the right hand sides of the subsystems.
        /// Its purpose is to aid in the parallel solution of the subsystems.
        std::vector<int> knownRhsIndices;
        /// This vector holds the indices of the unknown values of the right hand sides of the subsystems.
        std::vector<int> unknownRhsIndices;
    } subsystemStructure;

    /// GPU version of subsystemStructure
    typedef struct{
        cl::Buffer subsystemPointers;
        cl::Buffer nzIndices;
        cl::Buffer knownRhsIndices;
        cl::Buffer unknownRhsIndices;
    } subsystemStructureGPU;

    subsystemStructure lower, upper;
    subsystemStructureGPU d_lower, d_upper;

    /// An approximate inverse for L is computed by solving a small lower triangular system for each column of the main matrix.
    /// This function finds the structure of each of these subsystems and fills the 'lower' struct.
    void buildLowerSubsystemsStructures();

    /// An approximate inverse for U is computed by solving a small upper triangular system for each column of the main matrix.
    /// This function finds the structure of each of theses subsystems and fills the 'upper' struct.
    void buildUpperSubsystemsStructures();

public:
    BISAI(bool opencl_ilu_parallel, int verbosity);

    // set own Opencl variables, but also that of the bilu0 preconditioner
    void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue) override;

    // analysis, extract parallelism
    bool analyze_matrix(BlockedMatrix *mat) override;
    bool analyze_matrix(BlockedMatrix *mat, BlockedMatrix *jacMat) override;

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix *mat) override;
    bool create_preconditioner(BlockedMatrix *mat, BlockedMatrix *jacMat) override;

    // apply preconditioner, x = prec(y)
    void apply(const cl::Buffer& y, cl::Buffer& x) override;
};

/// Similar function to csrPatternToCsc. It gives an offset map from CSR to CSC instead of the full CSR to CSC conversion.
/// The map works as follows: if an element 'e' of the matrix is in the i-th position in the CSR representation, it will be
/// in the csrToCscOffsetMap[i]-th position in the CSC representation.
std::vector<int> buildCsrToCscOffsetMap(std::vector<int> colPointers, std::vector<int> rowIndices);

} // namespace Accelerator
} // namespace Opm

#endif
