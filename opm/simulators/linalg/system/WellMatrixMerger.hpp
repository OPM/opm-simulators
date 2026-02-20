#pragma once

#include "SystemTypes.hpp"

#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace Opm
{

template <class Row>
int
countRowEntries(const Row& row)
{
    return std::distance(row.begin(), row.end());
}

template<typename Scalar>
void
mergeWRMatrices(const std::vector<WRMatrixT<Scalar>>& matrices,
                WRMatrixT<Scalar>& mergedMatrix,
                const std::vector<std::vector<int>>& wellIndices,
                int numResDof)
{
    if (matrices.empty()) {
        return;
    }

    int numCols = numResDof;
    int totalRows = 0;
    for (const auto& matrix : matrices) {
        totalRows += matrix.N();
    }

    mergedMatrix.setSize(totalRows, numCols);
    mergedMatrix.setBuildMode(WRMatrixT<Scalar>::random);

    int rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            int nnz = countRowEntries(matrix[i]);
            mergedMatrix.setrowsize(rowOffset + i, nnz);
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endrowsizes();

    rowOffset = 0;
    assert(wellIndices.size() == matrices.size());
    for (size_t well = 0; well < matrices.size(); ++well) {
        const auto& matrix = matrices[well];
        const auto& cells = wellIndices[well];
        assert(cells.size() == matrix.M());
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                int cell = cells[colIt.index()];
                mergedMatrix.addindex(rowOffset + i, cell);
            }
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endindices();

    rowOffset = 0;
    for (size_t well = 0; well < matrices.size(); ++well) {
        const auto& matrix = matrices[well];
        const auto& cells = wellIndices[well];
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                int cell = cells[colIt.index()];
                mergedMatrix[rowOffset + i][cell] = *colIt;
            }
        }
        rowOffset += matrix.N();
    }
}

template<typename Scalar>
void
mergeRWMatrices(const std::vector<RWMatrixT<Scalar>>& matrices,
                RWMatrixT<Scalar>& mergedMatrix,
                const std::vector<std::vector<int>>& wellIndices,
                int numResDof)
{
    if (matrices.empty()) {
        return;
    }

    // All matrices should have the same number of rows (reservoir cells)
    int numRows = numResDof;

    // Count total columns (wells)
    int totalCols = 0;
    for (const auto& matrix : matrices) {
        totalCols += matrix.M();
    }

    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(numRows, totalCols);
    mergedMatrix.setBuildMode(RWMatrixT<Scalar>::random);

    // Phase 1: Set row sizes
    std::vector<int> row_sizes(numRows, 0);
    for (const auto& cells : wellIndices) {
        for (const auto& cell : cells) {
            row_sizes[cell] += 1; // maybe larger than need rowsize
        }
    }
    for (int row = 0; row < numRows; ++row) {
        mergedMatrix.setrowsize(row, row_sizes[row]);
    }
    mergedMatrix.endrowsizes();

    // Phase 2: Set column indices
    int colOffset = 0;
    for (size_t well = 0; well < matrices.size(); ++well) {
        const auto& matrix = matrices[well];
        const auto& cells = wellIndices[well];
        for (size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {
            int cell = cells[rowIdx];
            for (auto colIt = matrix[rowIdx].begin(); colIt != matrix[rowIdx].end(); ++colIt) {
                mergedMatrix.addindex(cell, colOffset + colIt.index());
            }
        }
        colOffset += matrix.M();
    }
    mergedMatrix.endindices();

    // Phase 3: Copy values
    colOffset = 0;
    for (size_t well = 0; well < matrices.size(); ++well) {
        const auto& matrix = matrices[well];
        const auto& cells = wellIndices[well];
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            int cell = cells[i];
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix[cell][colOffset + colIt.index()] = *colIt;
            }
        }
        colOffset += matrix.M();
    }
}

// Create diagonal block matrix from D matrices (well-to-well)
template<typename Scalar>
void
createDiagonalBlockMatrix(const std::vector<WWMatrixT<Scalar>>& matrices, WWMatrixT<Scalar>& mergedMatrix)
{
    if (matrices.empty()) {
        return;
    }

    int totalSize = 0;
    for (const auto& matrix : matrices) {
        totalSize += matrix.N();
    }

    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(totalSize, totalSize);
    mergedMatrix.setBuildMode(WWMatrixT<Scalar>::random);

    // Helper struct to find which original matrix a row belongs to
    struct MatrixLocation {
        int matrixIdx;
        int localRowIdx;
        int rowOffset;
    };

    auto findMatrixForRow = [&matrices](int globalRowIdx) -> MatrixLocation {
        int matrixIdx = 0;
        int localRowIdx = globalRowIdx;
        int rowOffset = 0;

        for (const auto& matrix : matrices) {
            if (localRowIdx < static_cast<int>(matrix.N())) {
                return {matrixIdx, localRowIdx, rowOffset};
            }
            localRowIdx -= matrix.N();
            rowOffset += matrix.N();
            matrixIdx++;
        }
        // Should never reach here if globalRowIdx is valid
        return {-1, -1, -1};
    };

    // Phase 1: Set row sizes
    for (int rowIdx = 0; rowIdx < totalSize; ++rowIdx) {
        auto loc = findMatrixForRow(rowIdx);
        if (loc.matrixIdx >= 0 && loc.matrixIdx < static_cast<int>(matrices.size())) {
            int nnz = countRowEntries(matrices[loc.matrixIdx][loc.localRowIdx]);
            mergedMatrix.setrowsize(rowIdx, nnz);
        }
    }
    mergedMatrix.endrowsizes();

    // Phase 2: Set column indices
    for (int rowIdx = 0; rowIdx < totalSize; ++rowIdx) {
        auto loc = findMatrixForRow(rowIdx);
        if (loc.matrixIdx >= 0 && loc.matrixIdx < static_cast<int>(matrices.size())) {
            const auto& origRow = matrices[loc.matrixIdx][loc.localRowIdx];
            for (auto colIt = origRow.begin(); colIt != origRow.end(); ++colIt) {
                mergedMatrix.addindex(rowIdx, loc.rowOffset + colIt.index());
            }
        }
    }
    mergedMatrix.endindices();

    // Phase 3: Copy values
    int rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            for (auto colIt = matrix[i].begin(); colIt != matrix[i].end(); ++colIt) {
                mergedMatrix[rowOffset + i][rowOffset + colIt.index()] = *colIt;
            }
        }
        rowOffset += matrix.N();
    }
}

// Class to handle well matrix merging
template<typename Scalar>
class WellMatrixMerger
{
public:
    explicit WellMatrixMerger(int numResDof)
        : numResDof_(numResDof)
    {
    }

    void addWell(const WRMatrixT<Scalar>& B,
                 const RWMatrixT<Scalar>& C,
                 const WWMatrixT<Scalar>& D,
                 const std::vector<int>& cellIndices,
                 int wellIndex,
                 const std::string& wellName,
                 const WellVectorT<Scalar>& residual)
    {
        B_matrices_.push_back(B);
        C_matrices_.push_back(C);
        D_matrices_.push_back(D);
        wellIndices_.push_back(wellIndex);
        wellCells_.push_back(cellIndices);
        wellNames_.push_back(wellName);
        well_residuals_.push_back(residual);

        // Map cell indices to well indices
        for (const auto& cellIdx : cellIndices) {
            cellToWellMap_[cellIdx] = wellIndex;
        }
    }

    // Finalize the merger by creating the merged matrices and residual
    void finalize()
    {
        mergeWRMatrices<Scalar>(B_matrices_, mergedB_, wellCells_, numResDof_);
        mergeRWMatrices<Scalar>(C_matrices_, mergedC_, wellCells_, numResDof_);
        createDiagonalBlockMatrix<Scalar>(D_matrices_, mergedD_);
        mergeWellResiduals();
    }

    WellVectorT<Scalar>& getMergedWellResidual() { return mergedWellResidual_; }
    const WellVectorT<Scalar>& getMergedWellResidual() const { return mergedWellResidual_; }

    std::vector<int>& getWellDofOffsets() { return wellDofOffsets_; }
    const std::vector<int>& getWellDofOffsets() const { return wellDofOffsets_; }

    WRMatrixT<Scalar>& getMergedB() { return mergedB_; }
    const WRMatrixT<Scalar>& getMergedB() const { return mergedB_; }

    RWMatrixT<Scalar>& getMergedC() { return mergedC_; }
    const RWMatrixT<Scalar>& getMergedC() const { return mergedC_; }

    WWMatrixT<Scalar>& getMergedD() { return mergedD_; }
    const WWMatrixT<Scalar>& getMergedD() const { return mergedD_; }

    int getWellIndexForCell(int cellIndex) const
    {
        auto it = cellToWellMap_.find(cellIndex);
        if (it != cellToWellMap_.end()) {
            return it->second;
        }
        return -1; // Not found
    }

    const std::string& getWellName(int wellIndex) const
    {
        static const std::string notFound = "Unknown";
        for (size_t i = 0; i < wellIndices_.size(); ++i) {
            if (wellIndices_[i] == wellIndex) {
                return wellNames_[i];
            }
        }
        return notFound;
    }

private:
    void mergeWellResiduals()
    {
        std::size_t totalSize = 0;
        wellDofOffsets_.clear();
        wellDofOffsets_.push_back(0);
        for (const auto& res : well_residuals_) {
            totalSize += res.N();
            wellDofOffsets_.push_back(static_cast<int>(totalSize));
        }
        mergedWellResidual_.resize(totalSize);
        std::size_t offset = 0;
        for (const auto& res : well_residuals_) {
            for (std::size_t i = 0; i < res.N(); ++i) {
                mergedWellResidual_[offset + i] = res[i];
            }
            offset += res.N();
        }
    }

    int numResDof_;
    std::vector<WRMatrixT<Scalar>> B_matrices_;
    std::vector<RWMatrixT<Scalar>> C_matrices_;
    std::vector<WWMatrixT<Scalar>> D_matrices_;
    std::vector<WellVectorT<Scalar>> well_residuals_;
    std::vector<int> wellIndices_;
    std::vector<std::vector<int>> wellCells_;
    std::vector<std::string> wellNames_;
    std::map<int, int> cellToWellMap_;

    WRMatrixT<Scalar> mergedB_;
    RWMatrixT<Scalar> mergedC_;
    WWMatrixT<Scalar> mergedD_;
    WellVectorT<Scalar> mergedWellResidual_;
    std::vector<int> wellDofOffsets_;
};

} // namespace Opm
