#pragma once
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
//#include <opm/simulators/linalg/WeightsCalculator.hpp>

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <map>
#include <string>
#include <functional>
#include <opm/simulators/linalg/matrixblock.hh>
// Define constants
const int numResDofs = 3;
const int numWellDofs = 4;
namespace Opm {
// Define matrix and vector types
using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;

// Define system matrix and vector types
using SystemMatrix = Dune::MultiTypeBlockMatrix<
    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
>;
using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

// Define helper constants for accessing MultiTypeBlockMatrix elements
constexpr auto _0 = Dune::Indices::_0;
constexpr auto _1 = Dune::Indices::_1;

// Helper function to create a matrix with a specific pattern and values
template<class Matrix, class Block>
void createMatrix(Matrix& matrix, int numRows, int numCols, const std::vector<int>& colIndices, const Block& blockValue)
{
    matrix.setSize(numRows, numCols);
    matrix.setBuildMode(Matrix::row_wise);
    
    // Define the sparsity pattern
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (const auto& col : colIndices) {
            if (col < numCols) {
                row.insert(col);
            }
        }
    }
    
    // Set the values
    for (int row = 0; row < numRows; ++row) {
        for (const auto& col : colIndices) {
            if (col < numCols) {
                matrix[row][col] = blockValue;
            }
        }
    }
}

// Helper function to print a matrix
template<class Matrix>
void printMatrix(const std::string& name, const Matrix& matrix)
{
    std::cout << "Matrix " << name << " (" << matrix.N() << " x " << matrix.M() << "):" << std::endl;
    
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        std::cout << "  Row " << row.index() << ": ";
        for (auto col = row->begin(); col != row->end(); ++col) {
            std::cout << "Col " << col.index() << " ";
            // Print the first element of each block
            std::cout << "(" << (*col)[0][0] << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Helper function to print a vector
template<class Vector>
void printVector(const std::string& name, const Vector& vector)
{
    std::cout << "Vector " << name << " (size " << vector.size() << "):" << std::endl;
    std::cout << "  [";
    for (size_t i = 0; i < vector.size(); ++i) {
        std::cout << "(";
        for (size_t j = 0; j < vector[i].size(); ++j) {
            std::cout << std::fixed << std::setprecision(4) << vector[i][j];
            if (j < vector[i].size() - 1) std::cout << ", ";
        }
        std::cout << ")";
        if (i < vector.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl << std::endl;
}


// Count non-zeros in a matrix row
template<class Row>
int countRowEntries(const Row& row) {
    return std::distance(row.begin(), row.end());
}

// Merge B matrices (well-to-reservoir) by vertical stacking
inline void mergeWRMatrices(const std::vector<WRMatrix>& matrices, WRMatrix& mergedMatrix, const std::vector<std::vector<int>>& wellIndices, int numResDof)
{
    if (matrices.empty()) {
        return;
    }
    
    // Count total rows and get number of columns
    int numCols = numResDof;
    int totalRows = 0;
    for (const auto& matrix : matrices) {
        totalRows += matrix.N();
    }
    
    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(totalRows, numCols);
    mergedMatrix.setBuildMode(WRMatrix::random);
    
    // Phase 1: Set row sizes
    int rowOffset = 0;
    for (const auto& matrix : matrices) {
        for (int i = 0; i < static_cast<int>(matrix.N()); ++i) {
            int nnz = countRowEntries(matrix[i]);
            mergedMatrix.setrowsize(rowOffset + i, nnz);
        }
        rowOffset += matrix.N();
    }
    mergedMatrix.endrowsizes();
    
    // Phase 2: Set column indices
    rowOffset = 0;
    assert(wellIndices.size() == matrices.size());
    for(size_t well=0; well<matrices.size(); ++well) {
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
    
    // Phase 3: Copy values
    rowOffset = 0;
    for(size_t well=0; well<matrices.size(); ++well) {
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

// Merge C matrices (reservoir-to-well) by horizontal stacking
inline void mergeRWMatrices(const std::vector<RWMatrix>& matrices, RWMatrix& mergedMatrix, const std::vector<std::vector<int>>& wellIndices, int numResDof)
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
    mergedMatrix.setBuildMode(RWMatrix::random);
    
    // Phase 1: Set row sizes
    std::vector<int> row_sizes(numRows, 0);
    for(const auto& cells : wellIndices) {
        for(const auto& cell : cells) {
            row_sizes[cell] += 1; // maybe larger than need rowsize
        }
    }
    for (int row = 0; row < numRows; ++row) {
        mergedMatrix.setrowsize(row, row_sizes[row]);
    }
    mergedMatrix.endrowsizes();
    
    // Phase 2: Set column indices
    int colOffset = 0;
    for(size_t well=0; well < matrices.size(); ++well) {
        const auto& matrix = matrices[well];  
        const auto& cells = wellIndices[well];
        for(size_t rowIdx = 0; rowIdx < matrix.N(); ++rowIdx) {     
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
    for(size_t well=0; well < matrices.size(); ++well) {
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
inline void createDiagonalBlockMatrix(const std::vector<WWMatrix>& matrices, WWMatrix& mergedMatrix)
{
    if (matrices.empty()) {
        return;
    }
    
    // Count total size
    int totalSize = 0;
    for (const auto& matrix : matrices) {
        totalSize += matrix.N();
    }
    
    // Set up the merged matrix with random build mode (required for setrowsize)
    mergedMatrix.setSize(totalSize, totalSize);
    mergedMatrix.setBuildMode(WWMatrix::random);
    
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
class WellMatrixMerger
{
public:
    WellMatrixMerger(int numResDof): numResDof_(numResDof) {}
    
    // Add a well to the merger
    void addWell(const WRMatrix& B, const RWMatrix& C, const WWMatrix& D, 
                 const std::vector<int>& cellIndices, int wellIndex, const std::string& wellName)
    {
        B_matrices_.push_back(B);
        C_matrices_.push_back(C);
        D_matrices_.push_back(D);
        wellIndices_.push_back(wellIndex);
        wellCells_.push_back(cellIndices);
        wellNames_.push_back(wellName);
         
        // Map cell indices to well indices
        for (const auto& cellIdx : cellIndices) {
            cellToWellMap_[cellIdx] = wellIndex;
        }
    }
    
    // Finalize the merger by creating the merged matrices
    void finalize()
    {
        // Merge B matrices (well-to-reservoir)
        mergeWRMatrices(B_matrices_, mergedB_, wellCells_, numResDof_);
        
        // Merge C matrices (reservoir-to-well)
        mergeRWMatrices(C_matrices_, mergedC_, wellCells_, numResDof_);
        
        // Create diagonal block matrix for D matrices (well-to-well)
        createDiagonalBlockMatrix(D_matrices_, mergedD_);
    }
    
    // Get the merged matrices
    const WRMatrix& getMergedB() const { return mergedB_; }
    const RWMatrix& getMergedC() const { return mergedC_; }
    const WWMatrix& getMergedD() const { return mergedD_; }
    
    // Get well index for a cell
    int getWellIndexForCell(int cellIndex) const
    {
        auto it = cellToWellMap_.find(cellIndex);
        if (it != cellToWellMap_.end()) {
            return it->second;
        }
        return -1; // Not found
    }
    
    // Get well name by index
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
    
    // Create system matrix
    SystemMatrix createSystemMatrix(const RRMatrix& A) const
    {
        // Create a system matrix with the correct structure
        SystemMatrix S;
        
        // Set the submatrices
        S[_0][_0] = A;
        S[_0][_1] = mergedC_;
        S[_1][_0] = mergedB_;
        S[_1][_1] = mergedD_;
        
        return S;
    }
    
private:
    int numResDof_;
    std::vector<WRMatrix> B_matrices_;
    std::vector<RWMatrix> C_matrices_;
    std::vector<WWMatrix> D_matrices_;
    std::vector<int> wellIndices_;
    std::vector<std::vector<int>> wellCells_;
    std::vector<std::string> wellNames_;
    std::map<int, int> cellToWellMap_;
    
    WRMatrix mergedB_;
    RWMatrix mergedC_;
    WWMatrix mergedD_;
};

// Helper for creating field matrices with a pattern
template<int ROWS, int COLS>
Dune::FieldMatrix<double, ROWS, COLS> makeBlock(const std::function<double(int, int)>& valueFn) {
    Dune::FieldMatrix<double, ROWS, COLS> block;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            block[i][j] = valueFn(i, j);
        }
    }
    return block;
}

// Helper for creating diagonal blocks with constant off-diagonal values
template<int SIZE>
Dune::FieldMatrix<double, SIZE, SIZE> makeDiagonalBlock(double diagValue, double offDiagValue = 0.0) {
    return makeBlock<SIZE, SIZE>([=](int i, int j) {
        return (i == j) ? diagValue : offDiagValue;
    });
}

// Helper for creating scaled blocks
template<int ROWS, int COLS>
Dune::FieldMatrix<double, ROWS, COLS> makeScaledBlock(double scale) {
    return makeBlock<ROWS, COLS>([=](int i, int j) {
        return scale * (i + 1) * (j + 1);
    });
}

// Create a tridiagonal matrix
inline void createTridiagonalMatrix(RRMatrix& matrix, int numRows, double diagValue, double offDiagValue) {
    matrix.setSize(numRows, numRows);
    matrix.setBuildMode(RRMatrix::random);
    
    // Phase 1: Set row sizes
    for (int i = 0; i < numRows; ++i) {
        int nnz = 1; // Diagonal element
        if (i > 0) nnz++; // Sub-diagonal
        if (i < numRows - 1) nnz++; // Super-diagonal
        matrix.setrowsize(i, nnz);
    }
    matrix.endrowsizes();
    
    // Phase 2: Set column indices
    for (int i = 0; i < numRows; ++i) {
        if (i > 0) {
            matrix.addindex(i, i-1); // Sub-diagonal
        }
        matrix.addindex(i, i); // Diagonal
        if (i < numRows - 1) {
            matrix.addindex(i, i+1); // Super-diagonal
        }
    }
    matrix.endindices();
    
    // Phase 3: Set values
    for (int i = 0; i < numRows; ++i) {
        if (i > 0) {
            matrix[i][i-1] = offDiagValue; // Sub-diagonal
        }
        matrix[i][i] = diagValue; // Diagonal
        if (i < numRows - 1) {
            matrix[i][i+1] = offDiagValue; // Super-diagonal
        }
    }
}
} // namespace Opm