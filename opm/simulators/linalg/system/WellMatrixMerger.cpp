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
#include <config.h>
#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <numeric>

namespace {

template<class Matrix>
bool hasSameMatrixSparsity(const Matrix& matrix,
                           const Opm::MatrixSparsityPattern& pattern)
{
    if (matrix.N() != pattern.rows || matrix.M() != pattern.cols) {
        return false;
    }

    if (pattern.rowOffsets.size() != pattern.rows + 1
        || pattern.rowOffsets.empty()
        || pattern.rowOffsets.front() != 0)
    {
        return false;
    }

    std::size_t entryOffset = 0;
    for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
        if (pattern.rowOffsets[rowIdx] != entryOffset) {
            return false;
        }

        for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
            if (entryOffset >= pattern.columnIndices.size()
                || static_cast<std::size_t>(pattern.columnIndices[entryOffset]) != colIt.index())
            {
                return false;
            }
            ++entryOffset;
        }

        if (pattern.rowOffsets[rowIdx + 1] != entryOffset) {
            return false;
        }
    }

    return entryOffset == pattern.columnIndices.size();
}

template<class Row>
std::size_t countRowEntries(const Row& row)
{
    return static_cast<std::size_t>(std::distance(row.begin(), row.end()));
}

template<class Matrix>
void initializeEmptyMatrix(Matrix& matrix, std::size_t rows, std::size_t cols)
{
    matrix.setSize(rows, cols);
    matrix.setBuildMode(Matrix::random);
    for (std::size_t row = 0; row < rows; ++row) {
        matrix.setrowsize(row, 0);
    }
    matrix.endrowsizes();
    matrix.endindices();
}

template<class MatrixVectorT, class DimensionFn>
std::size_t sumMatrixDimension(const MatrixVectorT& matrices, DimensionFn&& dimension)
{
    return std::accumulate(matrices.begin(), matrices.end(), std::size_t{0},
                           [&dimension](const auto acc, const auto& matrix)
                           { return acc + dimension(matrix); });
}

template<class Row, class Matrix, class ColumnMapper>
void assignRowValues(const Row& row,
                     Matrix& mergedMatrix,
                     std::size_t mergedRow,
                     ColumnMapper&& mapColumn)
{
    for (auto colIt = row.begin(); colIt != row.end(); ++colIt) {
        mergedMatrix[mergedRow][mapColumn(colIt.index())] = *colIt;
    }
}

}

namespace Opm {

template<typename Scalar>
WellMatrixMerger<Scalar>::
WellMatrixMerger(const std::size_t nResDofs,
                 const std::vector<BMatrix>& bMatrices,
                 const std::vector<CMatrix>& cMatrices,
                 const std::vector<DMatrix>& dMatrices,
                 const Opm::SparseTable<int>& wellCells)
    : numResDofs_(nResDofs)
    , bMatrices_(bMatrices)
    , cMatrices_(cMatrices)
    , dMatrices_(dMatrices)
    , wellCells_(wellCells)
{
    assert(inputsAreValid());
}

template<typename Scalar>
bool WellMatrixMerger<Scalar>::
hasSameStructure(const WellMatrixStructure& cachedStructure) const
{
    const auto numWells = bMatrices_.size();
    if (cachedStructure.numResDofs != numResDofs_
        || static_cast<std::size_t>(cachedStructure.wellCells.size()) != numWells
        || cachedStructure.bPatterns.size() != numWells
        || cachedStructure.cPatterns.size() != numWells
        || cachedStructure.dPatterns.size() != numWells)
    {
        return false;
    }

    std::size_t totalWellDofs = 0;
    for (std::size_t well = 0; well < numWells; ++well) {
        totalWellDofs += dMatrices_[well].N();

        const auto& cachedRow = cachedStructure.wellCells[well];
        const auto& currentRow = wellCells_[well];
        if (!std::ranges::equal(cachedRow.begin(), cachedRow.end(), currentRow.begin(), currentRow.end())
            || !hasSameMatrixSparsity(bMatrices_[well], cachedStructure.bPatterns[well])
            || !hasSameMatrixSparsity(cMatrices_[well], cachedStructure.cPatterns[well])
            || !hasSameMatrixSparsity(dMatrices_[well], cachedStructure.dPatterns[well]))
        {
            return false;
        }
    }

    return cachedStructure.totalWellBlocks == totalWellDofs;
}

template<typename Scalar>
WellMatrixStructure
WellMatrixMerger<Scalar>::
buildStructure() const
{
    WellMatrixStructure structure;
    structure.numResDofs = numResDofs_;
    structure.totalWellBlocks = 0;
    structure.wellCells = wellCells_;
    structure.bPatterns.reserve(bMatrices_.size());
    structure.cPatterns.reserve(cMatrices_.size());
    structure.dPatterns.reserve(dMatrices_.size());

    for (std::size_t well = 0; well < bMatrices_.size(); ++well) {
        structure.bPatterns.push_back(captureMatrixSparsity(bMatrices_[well]));
        structure.cPatterns.push_back(captureMatrixSparsity(cMatrices_[well]));
        structure.dPatterns.push_back(captureMatrixSparsity(dMatrices_[well]));
        structure.totalWellBlocks += dMatrices_[well].N();
    }

    return structure;
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
buildMatrices(BMatrix& mergedB,
              CMatrix& mergedC,
              DMatrix& mergedD) const
{
    mergeBMatrix(mergedB);
    mergeCMatrix(mergedC);
    mergeDMatrix(mergedD);
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
updateValues(BMatrix& mergedB,
             CMatrix& mergedC,
             DMatrix& mergedD) const
{
    fillBValues(mergedB);
    fillCValues(mergedC);
    fillDValues(mergedD);
}

template<typename Scalar>
bool
WellMatrixMerger<Scalar>::
inputsAreValid() const
{
    const auto numWells = bMatrices_.size();

    if (cMatrices_.size() != numWells) {
        return false;
    }

    if (dMatrices_.size() != numWells) {
        return false;
    }

    if (static_cast<std::size_t>(wellCells_.size()) != numWells) {
        return false;
    }

    for (std::size_t well = 0; well < numWells; ++well) {
        const auto& B = bMatrices_[well];
        const auto& C = cMatrices_[well];
        const auto& D = dMatrices_[well];
        const auto& cells = wellCells_[well];

        if (cells.size() != B.M()) {
            return false;
        }

        if (cells.size() != C.N()) {
            return false;
        }

        if (B.N() != C.M()) {
            return false;
        }

        if (B.N() != D.N()) {
            return false;
        }

        if (C.M() != D.M()) {
            return false;
        }

        if (D.N() != D.M()) {
            return false;
        }
    }

    return true;
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
fillBValues(BMatrix& mergedMatrix) const
{
    std::size_t rowOffset = 0;
    for (std::size_t well = 0; well < bMatrices_.size(); ++well) {
        const auto& matrix = bMatrices_[well];
        const auto& cells = wellCells_[well];
        for (std::size_t row = 0; row < matrix.N(); ++row) {
            assignRowValues(matrix[row], mergedMatrix, rowOffset + row,
                            [&cells](auto localColumn) {
                                assert(cells[localColumn] >= 0);
                                return static_cast<std::size_t>(cells[localColumn]);
                            });
        }
        rowOffset += matrix.N();
    }
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
fillCValues(CMatrix& mergedMatrix) const
{
    std::size_t colOffset = 0;
    for (std::size_t well = 0; well < cMatrices_.size(); ++well) {
        const auto& matrix = cMatrices_[well];
        const auto& cells = wellCells_[well];
        for (std::size_t row = 0; row < matrix.N(); ++row) {
            assert(cells[row] >= 0);
            const auto cell = static_cast<std::size_t>(cells[row]);
            assignRowValues(matrix[row], mergedMatrix, cell,
                            [colOffset](auto localColumn) { return colOffset + localColumn; });
        }
        colOffset += matrix.M();
    }
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
fillDValues(DMatrix& mergedMatrix) const
{
    std::size_t rowOffset = 0;
    for (const auto& matrix : dMatrices_) {
        for (std::size_t row = 0; row < matrix.N(); ++row) {
            assignRowValues(matrix[row], mergedMatrix, rowOffset + row,
                            [rowOffset](auto localColumn) { return rowOffset + localColumn; });
        }
        rowOffset += matrix.N();
    }
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
mergeBMatrix(BMatrix& mergedMatrix) const
{
    if (bMatrices_.empty()) {
        initializeEmptyMatrix(mergedMatrix, 0, numResDofs_);
        return;
    }

    const auto totalRows = sumMatrixDimension(bMatrices_, [](const auto& matrix) { return matrix.N(); });

    mergedMatrix.setSize(totalRows, numResDofs_);
    mergedMatrix.setBuildMode(BMatrix::random);

    std::size_t rowOffset = 0;
    for (const auto& matrix : bMatrices_) {
        for (std::size_t row = 0; row < matrix.N(); ++row) {
            mergedMatrix.setrowsize(rowOffset + row, countRowEntries(matrix[row]));
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endrowsizes();

    rowOffset = 0;
    for (std::size_t well = 0; well < bMatrices_.size(); ++well) {
        const auto& matrix = bMatrices_[well];
        const auto& cells = wellCells_[well];
        for (std::size_t row = 0; row < matrix.N(); ++row) {
            for (auto colIt = matrix[row].begin(); colIt != matrix[row].end(); ++colIt) {
                assert(cells[colIt.index()] >= 0);
                mergedMatrix.addindex(rowOffset + row,
                                      static_cast<std::size_t>(cells[colIt.index()]));
            }
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endindices();

    fillBValues(mergedMatrix);
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
mergeCMatrix(CMatrix& mergedMatrix) const
{
    if (cMatrices_.empty()) {
        initializeEmptyMatrix(mergedMatrix, numResDofs_, 0);
        return;
    }

    const auto totalCols = sumMatrixDimension(cMatrices_, [](const auto& matrix) { return matrix.M(); });

    mergedMatrix.setSize(numResDofs_, totalCols);
    mergedMatrix.setBuildMode(CMatrix::random);

    std::vector<std::size_t> rowSizes(numResDofs_, 0);
    for (std::size_t well = 0; well < cMatrices_.size(); ++well) {
        const auto& matrix = cMatrices_[well];
        const auto& cells = wellCells_[well];
        for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
            assert(cells[rowIdx] >= 0);
            rowSizes[static_cast<std::size_t>(cells[rowIdx])] += countRowEntries(matrix[rowIdx]);
        }
    }
    for (std::size_t row = 0; row < numResDofs_; ++row) {
        mergedMatrix.setrowsize(row, rowSizes[row]);
    }
    mergedMatrix.endrowsizes();

    std::size_t colOffset = 0;
    for (std::size_t well = 0; well < cMatrices_.size(); ++well) {
        const auto& matrix = cMatrices_[well];
        const auto& cells = wellCells_[well];
        for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
            assert(cells[rowIdx] >= 0);
            const auto cell = static_cast<std::size_t>(cells[rowIdx]);
            for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
                mergedMatrix.addindex(cell, colOffset + colIt.index());
            }
        }
        colOffset += matrix.M();
    }
    mergedMatrix.endindices();

    fillCValues(mergedMatrix);
}

template<typename Scalar>
void
WellMatrixMerger<Scalar>::
mergeDMatrix(DMatrix& mergedMatrix) const
{
    if (dMatrices_.empty()) {
        initializeEmptyMatrix(mergedMatrix, 0, 0);
        return;
    }

    const auto totalSize = sumMatrixDimension(dMatrices_, [](const auto& matrix) { return matrix.N(); });

    mergedMatrix.setSize(totalSize, totalSize);
    mergedMatrix.setBuildMode(DMatrix::random);

    std::size_t rowOffset = 0;
    for (const auto& matrix : dMatrices_) {
        for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
            mergedMatrix.setrowsize(rowOffset + rowIdx, countRowEntries(matrix[rowIdx]));
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endrowsizes();

    rowOffset = 0;
    for (const auto& matrix : dMatrices_) {
        for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
            for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
                mergedMatrix.addindex(rowOffset + rowIdx, rowOffset + colIt.index());
            }
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endindices();

    fillDValues(mergedMatrix);
}

template class WellMatrixMerger<double>;
#if FLOW_INSTANTIATE_FLOAT
template class WellMatrixMerger<float>;
#endif

} // namespace Opm
