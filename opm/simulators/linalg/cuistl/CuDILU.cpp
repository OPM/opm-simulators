/*
  Copyright 2022-2023 SINTEF AS

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
#include <cuda.h>
#include <cuda_runtime.h>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuDILU.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <vector>

namespace
{
std::vector<int>
createReorderedToNatural(Opm::SparseTable<size_t> levelSets)
{
    auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[globCnt++] = static_cast<int>(col);
        }
    }
    return res;
}

std::vector<int>
createNaturalToReordered(Opm::SparseTable<size_t> levelSets)
{
    auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[col] = globCnt++;
        }
    }
    return res;
}

// TODO: When this function is called we already have the natural ordered matrix on the GPU
// TODO: could it be possible to create the reordered one in a kernel to speed up the constructor?
template <class M, class field_type>
Opm::cuistl::CuSparseMatrix<field_type>
createReorderedMatrix(const M& naturalMatrix, std::vector<int> reorderedToNatural)
{
    M reorderedMatrix(naturalMatrix.N(), naturalMatrix.N(), naturalMatrix.nonzeroes(), M::row_wise);
    for (auto dstRowIt = reorderedMatrix.createbegin(); dstRowIt != reorderedMatrix.createend(); ++dstRowIt) {
        auto srcRow = naturalMatrix.begin() + reorderedToNatural[dstRowIt.index()];
        // For elements in A
        for (auto elem = srcRow->begin(); elem != srcRow->end(); elem++) {
            dstRowIt.insert(elem.index());
        }
    }

    // TODO: There is probably a faster way to copy by copying whole rows at a time
    for (auto dstRowIt = reorderedMatrix.begin(); dstRowIt != reorderedMatrix.end(); ++dstRowIt) {
        auto srcRow = naturalMatrix.begin() + reorderedToNatural[dstRowIt.index()];
        for (auto elem = srcRow->begin(); elem != srcRow->end(); elem++) {
            reorderedMatrix[dstRowIt.index()][elem.index()] = *elem;
        }
    }

    return Opm::cuistl::CuSparseMatrix<field_type>::fromMatrix(reorderedMatrix, true);
}

} // NAMESPACE

namespace Opm::cuistl
{

template <class M, class X, class Y, int l>
CuDILU<M, X, Y, l>::CuDILU(const M& A)
    : m_cpuMatrix(A)
    , m_levelSets(Opm::getMatrixRowColoring(m_cpuMatrix, Opm::ColoringType::LOWER))
    , m_reorderedToNatural(createReorderedToNatural(m_levelSets))
    , m_naturalToReordered(createNaturalToReordered(m_levelSets))
    , m_gpuMatrix(CuSparseMatrix<field_type>::fromMatrix(m_cpuMatrix, true))
    , m_gpuMatrixReordered(createReorderedMatrix<M, field_type>(m_cpuMatrix, m_reorderedToNatural))
    , m_gpuNaturalToReorder(m_naturalToReordered)
    , m_gpuReorderToNatural(m_reorderedToNatural)
    , m_gpuDInv(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize())

{
    // TODO: Should in some way verify that this matrix is symmetric, only do it debug mode?
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_gpuMatrix.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_gpuMatrix.N(), A.N()));
    OPM_ERROR_IF(A[0][0].N() != m_gpuMatrix.blockSize(),
                 fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.blockSize(),
                             A[0][0].N()));
    OPM_ERROR_IF(A.N() * A[0][0].N() != m_gpuMatrix.dim(),
                 fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.dim(),
                             A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != m_gpuMatrix.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_gpuMatrix.nonzeroes(),
                             A.nonzeroes()));
    update();
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::apply(X& v, const Y& d)
{
    OPM_TIMEBLOCK(prec_apply);
    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        detail::computeLowerSolveLevelSet<field_type, blocksize_>(m_gpuMatrixReordered.getNonZeroValues().data(),
                                                                  m_gpuMatrixReordered.getRowIndices().data(),
                                                                  m_gpuMatrixReordered.getColumnIndices().data(),
                                                                  m_gpuReorderToNatural.data(),
                                                                  levelStartIdx,
                                                                  numOfRowsInLevel,
                                                                  m_gpuDInv.data(),
                                                                  d.data(),
                                                                  v.data());
        levelStartIdx += numOfRowsInLevel;
    }

    levelStartIdx = m_cpuMatrix.N();
    //  upper triangular solve: (D + U_A) v = Dy
    for (int level = m_levelSets.size() - 1; level >= 0; --level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        levelStartIdx -= numOfRowsInLevel;
        detail::computeUpperSolveLevelSet<field_type, blocksize_>(m_gpuMatrixReordered.getNonZeroValues().data(),
                                                                  m_gpuMatrixReordered.getRowIndices().data(),
                                                                  m_gpuMatrixReordered.getColumnIndices().data(),
                                                                  m_gpuReorderToNatural.data(),
                                                                  levelStartIdx,
                                                                  numOfRowsInLevel,
                                                                  m_gpuDInv.data(),
                                                                  v.data());
    }
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
CuDILU<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::update()
{
    OPM_TIMEBLOCK(prec_update);

    m_gpuMatrix.updateNonzeroValues(m_cpuMatrix, true); // send updated matrix to the gpu

    detail::copyMatDataToReordered<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                           m_gpuMatrix.getRowIndices().data(),
                                                           m_gpuMatrixReordered.getNonZeroValues().data(),
                                                           m_gpuMatrixReordered.getRowIndices().data(),
                                                           m_gpuNaturalToReorder.data(),
                                                           m_gpuMatrixReordered.N());

    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();

        detail::computeDiluDiagonal<field_type, blocksize_>(m_gpuMatrixReordered.getNonZeroValues().data(),
                                                            m_gpuMatrixReordered.getRowIndices().data(),
                                                            m_gpuMatrixReordered.getColumnIndices().data(),
                                                            m_gpuReorderToNatural.data(),
                                                            m_gpuNaturalToReorder.data(),
                                                            levelStartIdx,
                                                            numOfRowsInLevel,
                                                            m_gpuDInv.data());
        levelStartIdx += numOfRowsInLevel;
    }
}

} // namespace Opm::cuistl
#define INSTANTIATE_CUDILU_DUNE(realtype, blockdim)                                                                    \
    template class ::Opm::cuistl::CuDILU<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,            \
                                         ::Opm::cuistl::CuVector<realtype>,                                            \
                                         ::Opm::cuistl::CuVector<realtype>>;                                           \
    template class ::Opm::cuistl::CuDILU<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,             \
                                         ::Opm::cuistl::CuVector<realtype>,                                            \
                                         ::Opm::cuistl::CuVector<realtype>>

INSTANTIATE_CUDILU_DUNE(double, 1);
INSTANTIATE_CUDILU_DUNE(double, 2);
INSTANTIATE_CUDILU_DUNE(double, 3);
INSTANTIATE_CUDILU_DUNE(double, 4);
INSTANTIATE_CUDILU_DUNE(double, 5);
INSTANTIATE_CUDILU_DUNE(double, 6);

INSTANTIATE_CUDILU_DUNE(float, 1);
INSTANTIATE_CUDILU_DUNE(float, 2);
INSTANTIATE_CUDILU_DUNE(float, 3);
INSTANTIATE_CUDILU_DUNE(float, 4);
INSTANTIATE_CUDILU_DUNE(float, 5);
INSTANTIATE_CUDILU_DUNE(float, 6);
