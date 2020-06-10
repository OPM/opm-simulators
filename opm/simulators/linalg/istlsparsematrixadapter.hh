// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::IstlSparseMatrixAdapter
 */
#ifndef EWOMS_ISTL_SPARSE_MATRIX_ADAPTER_HH
#define EWOMS_ISTL_SPARSE_MATRIX_ADAPTER_HH

#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Opm {
namespace Linear {

/*!
 * \ingroup Linear
 * \brief A sparse matrix interface backend for BCRSMatrix from dune-istl.
 */
template <class MatrixBlockType, class AllocatorType=std::allocator<MatrixBlockType> >
class IstlSparseMatrixAdapter
{
public:
    //! \brief Implementation of matrix
    using IstlMatrix = Dune::BCRSMatrix<MatrixBlockType, AllocatorType>;

    //! \brief block type forming the matrix entries
    using MatrixBlock = typename IstlMatrix::block_type;
    static_assert(std::is_same<MatrixBlock, MatrixBlockType>::value,
                  "IstlMatrix::block_type and MatrixBlockType must be identical");

    //! \brief type of scalar
    using Scalar = typename MatrixBlock::field_type;

    /*!
     * \brief Constructor creating an empty matrix.
     */
    IstlSparseMatrixAdapter(const size_t rows, const size_t columns)
        : rows_(rows)
        , columns_(columns)
        , istlMatrix_()
    {}

    /*!
     * \brief Constructor taking simulator and creating an empty matrix .
     */
    template <class Simulator>
    IstlSparseMatrixAdapter(const Simulator& simulator)
        : IstlSparseMatrixAdapter(simulator.model().numTotalDof(), simulator.model().numTotalDof())
    {}

    /*!
     * \brief Allocate matrix structure give a sparsity pattern.
     */
    template <class Set>
    void reserve(const std::vector<Set>& sparsityPattern)
    {
        // allocate raw matrix
        istlMatrix_.reset(new IstlMatrix(rows_, columns_, IstlMatrix::random));

        // make sure sparsityPattern is consistent with number of rows
        assert(rows_ == sparsityPattern.size());

        // allocate space for the rows of the matrix
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx)
            istlMatrix_->setrowsize(dofIdx, sparsityPattern[dofIdx].size());

        istlMatrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx) {
            auto nIt    = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                istlMatrix_->addindex(dofIdx, *nIt);
        }
        istlMatrix_->endindices();
    }

    /*!
     * \brief Return constant reference to matrix implementation.
     */
    IstlMatrix& istlMatrix()
    { return *istlMatrix_; }
    const IstlMatrix& istlMatrix() const
    { return *istlMatrix_; }

    /*!
     * \brief Return number of rows of the matrix.
     */
    size_t rows() const
    { return rows_; }

    /*!
     * \brief Return number of columns of the matrix.
     */
    size_t cols() const
    { return columns_; }

    /*!
     * \brief Set all matrix entries to zero.
     */
    void clear()
    { (*istlMatrix_) = Scalar(0.0); }

    /*!
     * \brief Set given row to zero except for the main-diagonal entry (if it exists).
     *
     * If the sparsity pattern of the matrix features an explicit block on the main
     * diagonal, the diagonal on that block is set to the second agument of the function.
     */
    void clearRow(const size_t row, const Scalar diag = 1.0)
    {
        MatrixBlock diagBlock(Scalar(0));
        for (int i = 0; i < diagBlock.rows; ++i)
            diagBlock[i][i] = diag;

        auto& matRow = (*istlMatrix_)[row];
        auto colIt = matRow.begin();
        const auto& colEndIt = matRow.end();
        for (; colIt != colEndIt; ++colIt) {
            if (colIt.index() == row)
                *colIt = diagBlock;
            else
                *colIt = Scalar(0.0);
        }
    }

    /*!
     * \brief Fill given block with entries stored in the matrix.
     */
    void block(const size_t rowIdx, const size_t colIdx, MatrixBlock& value) const
    { value = (*istlMatrix_)[rowIdx][colIdx]; }

    /*!
     * \brief Set matrix block to given block.
     */
    void setBlock(const size_t rowIdx, const size_t colIdx, const MatrixBlock& value)
    { (*istlMatrix_)[rowIdx][colIdx] = value; }

    /*!
     * \brief Add block to matrix block.
     */
    void addToBlock(const size_t rowIdx, const size_t colIdx, const MatrixBlock& value)
    { (*istlMatrix_)[rowIdx][colIdx] += value; }

    /*!
     * \brief Commit matrix from local caches into matrix native structure.
     *
     * For the ISTL adapter this is unnecessary because there is no caching mechanism.
     */
    void commit()
    { }

    /*!
     * \brief Finish modifying the matrix, i.e., convert the data structure from one
     *        tailored for linearization to one aimed at the linear solver.
     *
     * This may compress the matrix if the build mode is implicit. For the ISTL adapter
     * this is not required.
     */
    void finalize()
    { }

protected:
    size_t rows_;
    size_t columns_;

    std::unique_ptr<IstlMatrix> istlMatrix_;
};

}} // namespace Linear, Opm

#endif
