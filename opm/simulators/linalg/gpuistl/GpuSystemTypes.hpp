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
#ifndef OPM_GPUSYSTEMTYPES_HEADER_INCLUDED
#define OPM_GPUSYSTEMTYPES_HEADER_INCLUDED

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixGeneric.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

namespace Opm::gpuistl
{

// The reservoir block A is stored in the block-sparse GPU format.  The well
// blocks are expanded to scalar CSR matrices: the coupling blocks C (3x4
// blocks) and B (4x3 blocks) are not square, which the block-sparse format
// does not support, and all three are empty in a run without open wells,
// which only the scalar format handles.  D is only used for a matrix-vector
// product on the GPU; the well sub-solver works on its host copy.
template <typename Scalar>
using GpuRRMatrix = GpuSparseMatrixWrapper<Scalar>;
template <typename Scalar>
using GpuRWMatrix = GpuSparseMatrixGeneric<Scalar>;
template <typename Scalar>
using GpuWRMatrix = GpuSparseMatrixGeneric<Scalar>;
template <typename Scalar>
using GpuWWMatrix = GpuSparseMatrixGeneric<Scalar>;

// GPU counterpart of SystemVector.  SystemVector is a
// Dune::MultiTypeBlockVector, which cannot hold GPU vectors, so this is a
// plain pair of flat vectors with the vector operations the Dune Krylov
// solvers need.  The layout follows the CPU one: res holds numResDofs scalars
// per reservoir block, well numWellDofs per well block.  v[_0] and v[_1]
// access the two parts as for SystemVector, which lets SystemMatrix work on
// both.
template <typename Scalar>
struct GpuSystemVector {
    using field_type = Scalar;

    GpuVector<Scalar> res;
    GpuVector<Scalar> well;

    GpuSystemVector() = default;

    // Sizes are in scalars, not blocks.
    GpuSystemVector(std::size_t numResScalars, std::size_t numWellScalars)
        : res(numResScalars)
        , well(numWellScalars)
    {
    }

    GpuVector<Scalar>& operator[](Dune::index_constant<0>)
    {
        return res;
    }
    const GpuVector<Scalar>& operator[](Dune::index_constant<0>) const
    {
        return res;
    }
    GpuVector<Scalar>& operator[](Dune::index_constant<1>)
    {
        return well;
    }
    const GpuVector<Scalar>& operator[](Dune::index_constant<1>) const
    {
        return well;
    }

    // A run without open wells has an empty well part, on which the
    // GpuVector operations would throw.
    bool hasWells() const
    {
        return well.dim() > 0;
    }

    GpuSystemVector& operator=(Scalar s)
    {
        res = s;
        if (hasWells()) {
            well = s;
        }
        return *this;
    }

    GpuSystemVector& operator+=(const GpuSystemVector& other)
    {
        res += other.res;
        if (hasWells()) {
            well += other.well;
        }
        return *this;
    }

    GpuSystemVector& operator-=(const GpuSystemVector& other)
    {
        res -= other.res;
        if (hasWells()) {
            well -= other.well;
        }
        return *this;
    }

    GpuSystemVector& operator*=(Scalar alpha)
    {
        res *= alpha;
        if (hasWells()) {
            well *= alpha;
        }
        return *this;
    }

    // this += alpha * other
    void axpy(Scalar alpha, const GpuSystemVector& other)
    {
        res.axpy(alpha, other.res);
        if (hasWells()) {
            well.axpy(alpha, other.well);
        }
    }

    Scalar dot(const GpuSystemVector& other) const
    {
        const Scalar resDot = res.dot(other.res);
        return hasWells() ? resDot + well.dot(other.well) : resDot;
    }

    Scalar two_norm() const
    {
        return std::sqrt(dot(*this));
    }

    std::size_t dim() const
    {
        return res.dim() + well.dim();
    }
};

// Prints sizes only; there is no device-to-host copy.  Dune's GMRes solvers
// stream vectors on breakdown.
template <typename Scalar>
std::ostream&
operator<<(std::ostream& os, const GpuSystemVector<Scalar>& v)
{
    return os << "GpuSystemVector(res=" << v.res.dim() << ", well=" << v.well.dim() << ")";
}

// Block types of the GPU system matrix, for SystemMatrix.
template <typename Scalar>
struct GpuSystemBlockTypes {
    static constexpr bool isGpu = true;
    using RR = GpuRRMatrix<Scalar>;
    using RW = GpuRWMatrix<Scalar>;
    using WR = GpuWRMatrix<Scalar>;
    using WW = GpuWWMatrix<Scalar>;
    using Vector = GpuSystemVector<Scalar>;
};

// The system matrix over device blocks.  Its hostD member points to the host
// copy of D, which the CPU well sub-solver uses.
template <typename Scalar>
using GpuSystemMatrix = SystemMatrix<Scalar, GpuSystemBlockTypes<Scalar>>;

namespace detail
{

    // Row pointers and column indices of the scalar CSR expansion of a matrix
    // with R x C blocks.  Each block row expands to R scalar rows, and scalar row
    // i*R + r holds row r of every block in block row i, in block column order.
    template <typename Scalar, int R, int C>
    void scalarCsrStructure(const Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, R, C>>& blockMat,
                            std::vector<int>& rowPtrs,
                            std::vector<int>& colIndices)
    {
        const std::size_t numBlockRows = blockMat.N();
        rowPtrs.assign(numBlockRows * R + 1, 0);
        colIndices.clear();
        colIndices.reserve(blockMat.nonzeroes() * R * C);
        for (std::size_t bi = 0; bi < numBlockRows; ++bi) {
            for (int r = 0; r < R; ++r) {
                for (auto it = blockMat[bi].begin(); it != blockMat[bi].end(); ++it) {
                    for (int c = 0; c < C; ++c) {
                        colIndices.push_back(static_cast<int>(it.index()) * C + c);
                    }
                }
                rowPtrs[bi * R + r + 1] = static_cast<int>(colIndices.size());
            }
        }
    }

    // Values of the scalar CSR expansion, in the order of scalarCsrStructure().
    template <typename Scalar, int R, int C>
    void scalarCsrValues(const Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, R, C>>& blockMat,
                         std::vector<Scalar>& values)
    {
        values.clear();
        values.reserve(blockMat.nonzeroes() * R * C);
        for (std::size_t bi = 0; bi < blockMat.N(); ++bi) {
            for (int r = 0; r < R; ++r) {
                for (auto it = blockMat[bi].begin(); it != blockMat[bi].end(); ++it) {
                    for (int c = 0; c < C; ++c) {
                        values.push_back((*it)[r][c]);
                    }
                }
            }
        }
    }

} // namespace detail

// Copy a matrix with non-square R x C blocks to the GPU as a scalar CSR matrix
// with N()*R rows and M()*C columns.
template <typename Scalar, int R, int C>
GpuSparseMatrixGeneric<Scalar>
makeScalarGpuMatrix(const Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, R, C>>& blockMat)
{
    std::vector<int> rowPtrs;
    std::vector<int> colIndices;
    std::vector<Scalar> values;
    detail::scalarCsrStructure(blockMat, rowPtrs, colIndices);
    detail::scalarCsrValues(blockMat, values);
    return GpuSparseMatrixGeneric<Scalar>(values.data(),
                                          rowPtrs.data(),
                                          colIndices.data(),
                                          values.size(),
                                          /*blockSize=*/1,
                                          blockMat.N() * R,
                                          blockMat.M() * C);
}

// Refresh the values of a matrix made by makeScalarGpuMatrix() from a block
// matrix with the same sparsity pattern.  values is scratch storage.
template <typename Scalar, int R, int C>
void
updateScalarGpuMatrix(const Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, R, C>>& blockMat,
                      GpuSparseMatrixGeneric<Scalar>& gpuMat,
                      std::vector<Scalar>& values)
{
    detail::scalarCsrValues(blockMat, values);
    gpuMat.getNonZeroValues().copyFromHost(values);
}

} // namespace Opm::gpuistl

#endif // OPM_GPUSYSTEMTYPES_HEADER_INCLUDED
