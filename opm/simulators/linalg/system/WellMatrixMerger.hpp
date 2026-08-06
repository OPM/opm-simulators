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
#ifndef OPM_WELLMATRIXMERGER_HEADER_INCLUDED
#define OPM_WELLMATRIXMERGER_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/grid/utility/SparseTable.hpp>

#include <cstddef>
#include <vector>

namespace Opm
{

// Exact sparsity signature for one assembled well matrix.
//
// Matrix dimensions alone are not sufficient for MSW wells: the D block
// connectivity depends on the current segment topology, so we keep the full
// per-row column pattern and compare that directly.
struct MatrixSparsityPattern
{
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<std::size_t> rowOffsets;
    std::vector<int> columnIndices;

    bool operator==(const MatrixSparsityPattern& other) const
    {
        return rows == other.rows
            && cols == other.cols
            && rowOffsets == other.rowOffsets
            && columnIndices == other.columnIndices;
    }

    bool operator!=(const MatrixSparsityPattern& other) const
    {
        return !(*this == other);
    }
};

// Structural cache key for the merged well part of the system matrix.
//  totalWellBlocks is the sum of all individual well D-matrix dimensions,
//  i.e., the total number of well degrees of freedom.  It is stored here
//  so that the well-vector size and the structure-rebuild decision stay
//  in the same place.
struct WellMatrixStructure
{
    std::size_t numResDofs = 0;
    std::size_t totalWellBlocks = 0;  // Aggregated well DOFs (sum of D_i.N())
    Opm::SparseTable<int> wellCells;
    std::vector<MatrixSparsityPattern> bPatterns;
    std::vector<MatrixSparsityPattern> cPatterns;
    std::vector<MatrixSparsityPattern> dPatterns;

    bool operator==(const WellMatrixStructure& other) const
    {
        return numResDofs == other.numResDofs
            && totalWellBlocks == other.totalWellBlocks
            && wellCells == other.wellCells
            && bPatterns == other.bPatterns
            && cPatterns == other.cPatterns
            && dPatterns == other.dPatterns;
    }

    bool operator!=(const WellMatrixStructure& other) const
    {
        return !(*this == other);
    }

    // Whether everything sized by the wells still has the same size.  A change
    // that keeps these can be absorbed by rebuilding patterns; a change that
    // does not introduces unknowns the initial build never saw.
    bool hasSameDimensions(const WellMatrixStructure& other) const
    {
        return numResDofs == other.numResDofs
            && totalWellBlocks == other.totalWellBlocks
            && wellCells.size() == other.wellCells.size();
    }
};

template<class Matrix>
MatrixSparsityPattern
captureMatrixSparsity(const Matrix& matrix)
{
    MatrixSparsityPattern pattern;
    pattern.rows = matrix.N();
    pattern.cols = matrix.M();
    pattern.rowOffsets.reserve(matrix.N() + 1);
    if constexpr (requires { matrix.nonzeroes(); }) {
        pattern.columnIndices.reserve(matrix.nonzeroes());
    }
    pattern.rowOffsets.push_back(0);

    for (std::size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
        for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
            pattern.columnIndices.push_back(colIt.index());
        }
        pattern.rowOffsets.push_back(pattern.columnIndices.size());
    }

    return pattern;
}

// WellMatrixMerger assembles the global coupled well part of
//
//     [ A  C ]
//     [ B  D ]
//
// from the per-well blocks B_j, C_j and D_j. It preserves each well's local
// sparsity pattern and only does two structural operations: concatenate the
// well blocks and remap perforation-related rows/columns through the list of
// perforated reservoir cells for each well.

// Give each block a distinctive value pattern so it is easy to see where it
// ended up after merging.
template<typename Scalar>
class WellMatrixMerger
{
public:
    using BMatrix = WRMatrix<Scalar>;
    using CMatrix = RWMatrix<Scalar>;
    using DMatrix = WWMatrix<Scalar>;

    WellMatrixMerger(const std::size_t nResDofs,
                     const std::vector<BMatrix>& bMatrices,
                     const std::vector<CMatrix>& cMatrices,
                     const std::vector<DMatrix>& dMatrices,
                     const Opm::SparseTable<int>& wellCells);

    bool hasSameStructure(const WellMatrixStructure& cachedStructure) const;

    WellMatrixStructure buildStructure() const;

    void buildMatrices(BMatrix& mergedB,
                       CMatrix& mergedC,
                       DMatrix& mergedD) const;

    void updateValues(BMatrix& mergedB,
                      CMatrix& mergedC,
                      DMatrix& mergedD) const;

private:
    bool inputsAreValid() const;

    void fillBValues(BMatrix& mergedMatrix) const;
    void fillCValues(CMatrix& mergedMatrix) const;
    void fillDValues(DMatrix& mergedMatrix) const;

    void mergeBMatrix(BMatrix& mergedMatrix) const;
    void mergeCMatrix(CMatrix& mergedMatrix) const;
    void mergeDMatrix(DMatrix& mergedMatrix) const;

    std::size_t numResDofs_;
    const std::vector<BMatrix>& bMatrices_;
    const std::vector<CMatrix>& cMatrices_;
    const std::vector<DMatrix>& dMatrices_;
    const Opm::SparseTable<int>& wellCells_;
};

} // namespace Opm

#endif // OPM_WELLMATRIXMERGER_HEADER_INCLUDED
