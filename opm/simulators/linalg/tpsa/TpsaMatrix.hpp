// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef OPM_TPSA_MATRIX_HPP
#define OPM_TPSA_MATRIX_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/tpsa/TpsaTypes.hpp>

#include <array>
#include <cstddef>
#include <numeric>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm::Linear
{

template <class Scalar>
class TpsaMatrix;

/*!
 * \brief Handle on one block of the TPSA matrix
 *
 * Wrapper around MatrixBlock to handle the distribution of incoming 7x7 dense block to the TPSA
 * sub-matrix fields.
 *
 * \tparam Scalar Field type of the matrix entries.
 */
template <class Scalar>
class TpsaBlockRef
{
public:
    using MatrixBlock = Opm::MatrixBlock<Scalar, numTpsaEq, numTpsaEq>;

    //! \brief Construct a handle that does not refer to any block.
    TpsaBlockRef() = default;

    /*!
     * \brief Construct a handle on one block of a TPSA matrix.
     *
     * \param[in] matrix Matrix owning the block
     * \param[in] flatIdx Flat index of the block within the shared sparsity pattern, as
     *                    returned by TpsaMatrix::flatIndex_().
     */
    TpsaBlockRef(const TpsaMatrix<Scalar>& matrix, std::size_t flatIdx)
        : matrix_(&matrix)
        , k_(flatIdx)
    {
    }

    /*!
     * \brief Dereference the handle, so that it can be used where a pointer to a block is
     *        expected.
     *
     * \return Reference to this handle itself.
     */
    TpsaBlockRef& operator*()
    {
        return *this;
    }

    /*!
     * \brief Scatter a dense 7x7 contribution into the sub-matrices.
     *
     * \param[in] b Dense block whose entries are added to the stored entries.
     *
     * \return Reference to this handle.
     */
    TpsaBlockRef& operator+=(const MatrixBlock& b);

    /*!
     * \brief Overwrite the stored entries with a dense 7x7 block.
     *
     * \param[in] b Dense block the stored entries are copied from.
     *
     * \return Reference to this handle.
     */
    TpsaBlockRef& operator=(const MatrixBlock& b);

    /*!
     * \brief Set every stored entry of this block to a scalar.
     *
     * \param[in] value Value assigned to all stored entries of the block.
     *
     * \return Reference to this handle.
     */
    TpsaBlockRef& operator=(Scalar value);

    /*!
     * \brief Gather the stored entries into a dense 7x7 block.
     *
     * \param[out] b Dense block the stored entries are written to. Fully overwritten.
     */
    void gather(MatrixBlock& b) const;

private:
    /*!
     * \brief Run an operation over every stored entry of this block, paired with the
     *        corresponding entry of a dense 7x7 block.
     *
     * Templated on the block type so that the same index map serves both the scattering
     * (const block) and the gathering (mutable block) direction.
     *
     * \tparam Op Callable invoked as op(storedEntry, denseEntry).
     * \tparam Block Dense block type, const for scattering and mutable for gathering.
     *
     * \param[in] op Operation applied to each (stored, dense) entry pair.
     * \param[in,out] b Dense block that supplies or receives the entries, depending on \p op.
     */
    template <class Op, class Block>
    void apply_(Op op, Block& b) const;

    const TpsaMatrix<Scalar>* matrix_{nullptr};
    std::size_t k_{0};
};

/*!
 * \brief TPSA matrix for linear elasticity.
 *
 * This is the equivalent as IstlSparseMatrixAdapter class is for flow linearizations.
*  The TPSA linearizer still provides dense 7x7 blocks, but internally in the class the entries go
*  into 19 sub-matrices. Displacement-displacement sub-matrix have been divided up in 1x1 fields
*  for the Hypre BoomerAMG preconditioner.
*  Overview of sub-matrix fields:
 *
 *     field 0: u_x            (1 dof)
 *     field 1: u_y            (1 dof)
 *     field 2: u_z            (1 dof)
 *     field 3: rot_x/y/z      (3 dofs)
 *     field 4: p_solid        (1 dof)
 *
 * Note that, displacement-displacement off-diagonals are zero in linear elasticity
 *
 * \tparam Scalar Field type of the matrix entries.
 */
template <class Scalar>
class TpsaMatrix
{
    friend class TpsaBlockRef<Scalar>;

public:
    //! \brief What the linear solver operates on.
    using IstlMatrix = TpsaMatrixView<Scalar>;

    //! \brief Dense local block the linearizer accumulates into.
    using MatrixBlock = Opm::MatrixBlock<Scalar, numTpsaEq, numTpsaEq>;

    //! \brief What blockAddress() returns.
    using BlockAddress = TpsaBlockRef<Scalar>;

    //! \brief Field type of the matrix entries.
    using field_type = Scalar;

    /*!
     * \brief Construct a matrix of the given block dimensions.
     *
     * No storage is allocated here; call reserve() with the sparsity pattern first.
     *
     * \param[in] rows Number of block rows, i.e. degrees of freedom.
     * \param[in] columns Number of block columns, i.e. degrees of freedom.
     */
    TpsaMatrix(std::size_t rows, std::size_t columns)
        : rows_(rows)
        , columns_(columns)
    {
    }

    /*!
     * \brief Construct a square matrix sized from a simulator's degrees of freedom.
     *
     * \tparam Simulator Simulator type exposing model().numTotalDof().
     *
     * \param[in] simulator Simulator the number of degrees of freedom is taken from. Only
     *                      read during construction.
     */
    template <class Simulator>
    explicit TpsaMatrix(const Simulator& simulator)
        : TpsaMatrix(simulator.model().numTotalDof(), simulator.model().numTotalDof())
    {
    }

    // The cached value-array base pointers and the sub-matrix view point into
    // this object, so it must not be copied or moved.
    TpsaMatrix(const TpsaMatrix&) = delete;

    TpsaMatrix(TpsaMatrix&&) = delete;

    TpsaMatrix& operator=(const TpsaMatrix&) = delete;

    TpsaMatrix& operator=(TpsaMatrix&&) = delete;

    /*!
     * \brief Allocate all sub-matrices from a common sparsity pattern.
     *
     * Also flattens the pattern, caches the value-array base pointer of every sub-matrix and
     * sets up the solver view. Must be called before any block is accessed.
     *
     * \tparam Set Ordered container of column indices, e.g. std::set<unsigned>.
     *
     * \param[in] sparsityPattern One entry per block row, holding that row's column indices
     *                            in ascending order.
     */
    template <class Set>
    void reserve(const std::vector<Set>& sparsityPattern)
    {
        if (sparsityPattern.size() != rows_) {
            OPM_THROW(std::logic_error,
                      "TPSA: sparsity pattern does not match the number of matrix rows");
        }

        // Flatten the pattern once; the column indices of a std::set are
        // already ascending, which is the order BCRSMatrix stores them in.
        rowStart_.resize(rows_ + 1);
        rowStart_[0] = 0;
        const auto sizes = sparsityPattern
            | std::views::transform([](const auto& a) { return a.size(); });
        std::partial_sum(sizes.begin(), sizes.end(), rowStart_.begin() + 1);
        nnz_ = rowStart_[rows_];

        colIdx_.clear();
        colIdx_.reserve(nnz_);
        for (std::size_t row = 0; row < rows_; ++row) {
            colIdx_.insert(colIdx_.end(),
                           sparsityPattern[row].begin(),
                           sparsityPattern[row].end());
        }

        forEachSubMatrix_([&](auto& subMatrix) {
            reserveSubMatrix_(subMatrix, sparsityPattern);
        });

        // Initialize base_ pointer to sub-matrices and set up TpsaMatrixView
        cacheValueArrays_();
        setMatrixView_();
    }

    /*!
     * \brief Handle on the block at (rowIdx, colIdx).
     *
     * Only called while the sparsity pattern is set up, so the linear scan over
     * the row is not on any hot path.
     *
     * \param[in] rowIdx Block row index.
     * \param[in] colIdx Block column index.
     *
     * \return Handle on the requested block.
     *
     * \throw std::logic_error If the block is not part of the sparsity pattern.
     */
    BlockAddress blockAddress(const std::size_t rowIdx, const std::size_t colIdx) const
    {
        return BlockAddress(*this, flatIndex_(rowIdx, colIdx));
    }

    /*!
     * \brief Set all matrix entries to zero.
     *
     * Does nothing when called before reserve().
     */
    void clear();

    /*!
     * \brief Set the given row to zero, except for the main diagonal.
     *
     * The main diagonal of the block on the diagonal is set to \p diag. Written
     * as dense blocks so that the field split is applied by TpsaBlockRef, i.e. by
     * the same index map assembly goes through, rather than by a second
     * description of which slots are field-diagonal.
     *
     * \param[in] row Block row to clear.
     * \param[in] diag Value put on the main diagonal of the diagonal block.
     */
    void clearRow(const std::size_t row, const Scalar diag = 1.0);

    /*!
     * \brief Zero out the overlap rows and put the identity on their diagonal.
     *
     * \param[in] overlapRows Block rows that are not owned by this process.
     */
    void makeOverlapRowsInvalid(const std::vector<int>& overlapRows);

    /*!
     * \brief Scale the system by one factor per field: A <- D_row * A * D_col.
     *
     * Every sub-matrix couples exactly one field row to one field column, so a
     * field scaling is a single factor per sub-matrix and can be applied to the
     * flat value arrays without walking the sparsity pattern.
     *
     * \param[in] rowFac Row (equation) factor of each field.
     * \param[in] colFac Column (unknown) factor of each field.
     *
     * \note Scales in place, so it must be called once per assembly. Does nothing when called
     *       before reserve().
     */
    void scaleFields(const std::array<Scalar, numTpsaFields>& rowFac,
                     const std::array<Scalar, numTpsaFields>& colFac);

    /*!
     * \brief Fill \p value with the stored entries of the given block.
     *
     * The displacement-displacement off-diagonals are not stored and come back as zero.
     *
     * \param[in] rowIdx Block row index.
     * \param[in] colIdx Block column index.
     * \param[out] value Dense block the entries are written to. Fully overwritten.
     */
    void block(const std::size_t rowIdx, const std::size_t colIdx, MatrixBlock& value) const
    {
        blockAddress(rowIdx, colIdx).gather(value);
    }

    /*!
     * \brief Overwrite the given block with a dense 7x7 block.
     *
     * \param[in] rowIdx Block row index.
     * \param[in] colIdx Block column index.
     * \param[in] value Dense block the stored entries are copied from. Its
     *                  displacement-displacement off-diagonals are ignored.
     */
    void setBlock(const std::size_t rowIdx, const std::size_t colIdx, const MatrixBlock& value)
    {
        blockAddress(rowIdx, colIdx) = value;
    }

    /*!
     * \brief Add a dense 7x7 block to the given block.
     *
     * \param[in] rowIdx Block row index.
     * \param[in] colIdx Block column index.
     * \param[in] value Dense block added to the stored entries. Its
     *                  displacement-displacement off-diagonals are ignored.
     */
    void addToBlock(const std::size_t rowIdx, const std::size_t colIdx, const MatrixBlock& value)
    {
        blockAddress(rowIdx, colIdx) += value;
    }

    //! \brief No local caching, so nothing to commit.
    void commit()
    {
    }

    //! \brief The structure is already solver-ready after reserve().
    void finalize()
    {
    }

    /*!
     * \brief The sub-matrix view the linear solver operates on.
     *
     * Only valid after reserve(); the view points into this object.
     *
     * \return Reference to the 5x5 view over the sub-matrices.
     */
    IstlMatrix& istlMatrix()
    {
        return view_;
    }

    //! \copydoc istlMatrix()
    const IstlMatrix& istlMatrix() const
    {
        return view_;
    }

    /*!
     * \brief Number of block rows.
     *
     * \return Number of block rows of the matrix.
     */
    std::size_t rows() const
    {
        return rows_;
    }

    /*!
     * \brief Number of block columns.
     *
     * \return Number of block columns of the matrix.
     */
    std::size_t cols() const
    {
        return columns_;
    }

    //! \copydoc rows()
    std::size_t N() const
    {
        return rows_;
    }

    //! \copydoc cols()
    std::size_t M() const
    {
        return columns_;
    }

    /*!
     * \brief Number of nonzero blocks in the shared sparsity pattern.
     *
     * \return Number of nonzero blocks, or zero before reserve() has been called.
     */
    std::size_t nonzeroes() const
    {
        return nnz_;
    }

    // Sub-matrix accessors.  dd00/dd11/dd22 and spsp are the scalar blocks Hypre
    // can precondition.

    /*!
     * \brief Access the u_x-u_x sub-matrix.
     *
     * \return Reference to the sub-matrix. Only meaningful after reserve().
     */
    DispDispMatrix00T<Scalar>& dd00()
    {
        return dd00_;
    }

    /*!
     * \brief Access the u_y-u_y sub-matrix.
     *
     * \return Reference to the sub-matrix. Only meaningful after reserve().
     */
    DispDispMatrix11T<Scalar>& dd11()
    {
        return dd11_;
    }

    /*!
     * \brief Access the u_z-u_z sub-matrix.
     *
     * \return Reference to the sub-matrix. Only meaningful after reserve().
     */
    DispDispMatrix22T<Scalar>& dd22()
    {
        return dd22_;
    }

    /*!
     * \brief Access the rotation-rotation sub-matrix.
     *
     * \return Reference to the sub-matrix. Only meaningful after reserve().
     */
    RotRotMatrixT<Scalar>& rr()
    {
        return rr_;
    }

    /*!
     * \brief Access the solid pressure-solid pressure sub-matrix.
     *
     * \return Reference to the sub-matrix. Only meaningful after reserve().
     */
    SPresSPresMatrixT<Scalar>& spsp()
    {
        return spsp_;
    }

    //! \copydoc dd00()
    const DispDispMatrix00T<Scalar>& dd00() const
    {
        return dd00_;
    }

    //! \copydoc dd11()
    const DispDispMatrix11T<Scalar>& dd11() const
    {
        return dd11_;
    }

    //! \copydoc dd22()
    const DispDispMatrix22T<Scalar>& dd22() const
    {
        return dd22_;
    }

    //! \copydoc rr()
    const RotRotMatrixT<Scalar>& rr() const
    {
        return rr_;
    }

    //! \copydoc spsp()
    const SPresSPresMatrixT<Scalar>& spsp() const
    {
        return spsp_;
    }

private:
    /*!
     * \brief Slots in the base_ array, one per sub-matrix.
     *
     * The order must match the tuple returned by subMatrices_().
     */
    enum SubMatrixIdx : std::size_t
    {
        DD00 = 0, DD11, DD22,
        DR0, DR1, DR2,
        DSP0, DSP1, DSP2,
        RD0, RD1, RD2,
        RR, RSP,
        SPD0, SPD1, SPD2,
        SPR, SPSP,
        numSubMatrices
    };

    /*!
     * \brief The (field row, field column) each sub-matrix couples.
     *
     * The order must match the SubMatrixIdx slots above, i.e. the tuple returned
     * by subMatrices_().
     */
    static constexpr std::array<std::pair<std::size_t, std::size_t>, numSubMatrices>
    subMatrixFields_ {{
        {0, 0}, {1, 1}, {2, 2},   // DD00, DD11, DD22
        {0, 3}, {1, 3}, {2, 3},   // DR0,  DR1,  DR2
        {0, 4}, {1, 4}, {2, 4},   // DSP0, DSP1, DSP2
        {3, 0}, {3, 1}, {3, 2},   // RD0,  RD1,  RD2
        {3, 3}, {3, 4},           // RR,   RSP
        {4, 0}, {4, 1}, {4, 2},   // SPD0, SPD1, SPD2
        {4, 3}, {4, 4}            // SPR,  SPSP
    }};

    /*!
     * \brief Scalars per block of a sub-matrix, i.e. the stride of its flat
     *        value array.
     *
     * Read off the sub-matrix' own block type, so there is no parallel table to
     * keep in step with the slot order.
     *
     * \tparam SubMatrix Sub-matrix type whose block_type carries the block dimensions.
     *
     * \param[in] (unnamed) Sub-matrix the stride is read from. Only its type is used.
     *
     * \return Number of scalars per block of that sub-matrix.
     */
    template <class SubMatrix>
    static constexpr std::size_t blockScalars_(const SubMatrix&)
    {
        using Block = typename SubMatrix::block_type;

        return static_cast<std::size_t>(Block::rows) * static_cast<std::size_t>(Block::cols);
    }

    /*!
     * \brief All sub-matrices in slot order.
     *
     * \return Tuple of references to the sub-matrices, ordered as in SubMatrixIdx.
     */
    auto subMatrices_()
    {
        return std::tie(dd00_,
                        dd11_,
                        dd22_,
                        dr0_,
                        dr1_,
                        dr2_,
                        dsp0_,
                        dsp1_,
                        dsp2_,
                        rd0_,
                        rd1_,
                        rd2_,
                        rr_,
                        rsp_,
                        spd0_,
                        spd1_,
                        spd2_,
                        spr_,
                        spsp_);
    }

    /*!
     * \brief Apply an operation to every sub-matrix.
     *
     * \tparam Op Callable invoked as op(subMatrix).
     *
     * \param[in] op Operation applied to each sub-matrix in slot order.
     */
    template <class Op>
    void forEachSubMatrix_(Op op)
    {
        std::apply([&op](auto&... subMatrix) {
                       (op(subMatrix), ...);
                   },
                   subMatrices_());
    }

    /*!
     * \brief Apply an operation to every sub-matrix together with its slot index.
     *
     * \tparam Op Callable invoked as op(subMatrix, subIdx).
     *
     * \param[in] op Operation applied to each sub-matrix and its SubMatrixIdx slot, in slot
     *               order.
     */
    template <class Op>
    void forEachSubMatrixWithIndex_(Op op)
    {
        std::apply([&op, subIdx = std::size_t{0}](auto&... subMatrix) mutable {
                       (op(subMatrix, subIdx++), ...);
                   },
                   subMatrices_());
    }

    /*!
     * \brief Build one sub-matrix from the shared sparsity pattern.
     *
     * All entries are zero afterwards.
     *
     * \tparam SubMatrix BCRSMatrix type of the sub-matrix.
     * \tparam Set Ordered container of column indices.
     *
     * \param[out] subMatrix Sub-matrix that is sized and given its sparsity pattern.
     * \param[in] sparsityPattern One entry per block row, holding that row's column indices
     *                            in ascending order.
     */
    template <class SubMatrix, class Set>
    void reserveSubMatrix_(SubMatrix& subMatrix, const std::vector<Set>& sparsityPattern)
    {
        subMatrix.setBuildMode(SubMatrix::random);
        subMatrix.setSize(rows_, columns_);

        for (std::size_t row = 0; row < rows_; ++row) {
            subMatrix.setrowsize(row, sparsityPattern[row].size());
        }
        subMatrix.endrowsizes();

        for (std::size_t row = 0; row < rows_; ++row) {
            for (const auto& col : sparsityPattern[row]) {
                subMatrix.addindex(row, col);
            }
        }
        // Note: all entries in subMatrix are Scalar(0.0) by default construction in endindices()
        subMatrix.endindices();
    }

    /*!
     * \brief Cache the base pointer of every sub-matrix' value array.
     *
     * Also verifies that the k-th nonzero really is at `base + k*stride` in
     * every sub-matrix, i.e. that the sub-matrices share both the pattern and
     * the storage order.  Everything else in this class rests on that.
     */
    void cacheValueArrays_();

    //! \brief Point the solver view at the sub-matrices owned by this object.
    void setMatrixView_();

    /*!
     * \brief Look up the flat index of a block in the shared sparsity pattern.
     *
     * \param[in] rowIdx Block row index.
     * \param[in] colIdx Block column index.
     *
     * \return Flat index of the block, counting nonzeroes row by row.
     */
    std::size_t flatIndex_(std::size_t rowIdx, std::size_t colIdx) const;

    std::size_t rows_{0};
    std::size_t columns_{0};
    std::size_t nnz_{0};

    // Flattened sparsity pattern, shared by all sub-matrices.
    std::vector<std::size_t> rowStart_{};
    std::vector<unsigned> colIdx_{};

    // Base pointers into the sub-matrices' contiguous value arrays.
    std::array<Scalar*, numSubMatrices> base_{};

    DispDispMatrix00T<Scalar> dd00_{};
    DispDispMatrix11T<Scalar> dd11_{};
    DispDispMatrix22T<Scalar> dd22_{};

    DispRotMatrix0T<Scalar> dr0_{};
    DispRotMatrix1T<Scalar> dr1_{};
    DispRotMatrix2T<Scalar> dr2_{};

    DispSPresMatrix0T<Scalar> dsp0_{};
    DispSPresMatrix1T<Scalar> dsp1_{};
    DispSPresMatrix2T<Scalar> dsp2_{};

    RotDispMatrix0T<Scalar> rd0_{};
    RotDispMatrix1T<Scalar> rd1_{};
    RotDispMatrix2T<Scalar> rd2_{};

    RotRotMatrixT<Scalar> rr_{};
    RotSPresMatrixT<Scalar> rsp_{};

    SPresDispMatrix0T<Scalar> spd0_{};
    SPresDispMatrix1T<Scalar> spd1_{};
    SPresDispMatrix2T<Scalar> spd2_{};

    SPresRotMatrixT<Scalar> spr_{};

    SPresSPresMatrixT<Scalar> spsp_{};

    IstlMatrix view_{};
};

//
// TpsaBlockRef implementation
//
template <class Scalar>
template <class Op, class Block>
void
TpsaBlockRef<Scalar>::apply_(Op op, Block& b) const
{
    using Matrix = TpsaMatrix<Scalar>;
    const auto& base = matrix_->base_;
    const std::size_t k = k_;

    // Displacement-displacement (diagonal components only)
    op(base[Matrix::DD00][k], b[0][0]);
    op(base[Matrix::DD11][k], b[1][1]);
    op(base[Matrix::DD22][k], b[2][2]);

    // Displacement-rotation: three 1x3 blocks
    for (int j = 0; j < numRotDofs; ++j) {
        op(base[Matrix::DR0][3 * k + j], b[0][3 + j]);
        op(base[Matrix::DR1][3 * k + j], b[1][3 + j]);
        op(base[Matrix::DR2][3 * k + j], b[2][3 + j]);
    }

    // Displacement-solid pressure: three 1x1 blocks
    op(base[Matrix::DSP0][k], b[0][6]);
    op(base[Matrix::DSP1][k], b[1][6]);
    op(base[Matrix::DSP2][k], b[2][6]);

    // Rotation-displacement: three 3x1 blocks
    for (int i = 0; i < numRotDofs; ++i) {
        op(base[Matrix::RD0][3 * k + i], b[3 + i][0]);
        op(base[Matrix::RD1][3 * k + i], b[3 + i][1]);
        op(base[Matrix::RD2][3 * k + i], b[3 + i][2]);
    }

    // Rotation-rotation: one 3x3 block
    for (int i = 0; i < numRotDofs; ++i) {
        for (int j = 0; j < numRotDofs; ++j) {
            op(base[Matrix::RR][9 * k + 3 * i + j], b[3 + i][3 + j]);
        }
    }

    // Rotation-solid pressure: one 3x1 block
    for (int i = 0; i < numRotDofs; ++i) {
        op(base[Matrix::RSP][3 * k + i], b[3 + i][6]);
    }

    // Solid pressure-displacement: three 1x1 blocks
    op(base[Matrix::SPD0][k], b[6][0]);
    op(base[Matrix::SPD1][k], b[6][1]);
    op(base[Matrix::SPD2][k], b[6][2]);

    // Solid pressure-rotation: one 1x3 block
    for (int j = 0; j < numRotDofs; ++j) {
        op(base[Matrix::SPR][3 * k + j], b[6][3 + j]);
    }

    // Solid pressure-solid pressure: one 1x1 block
    op(base[Matrix::SPSP][k], b[6][6]);
}

} // namespace Opm::Linear

#endif // OPM_TPSA_MATRIX_HPP
