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
#ifndef OPM_GPUDILU_HPP
#define OPM_GPUDILU_HPP

#include <memory>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/detail/kernel_enums.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_resources.hpp>
#include <vector>
#include <map>
#include <utility>


namespace Opm::gpuistl
{
//! \brief DILU preconditioner on the GPU.
//!
//! \tparam CPUMatrixT Type of the matrix on the CPU
//! \tparam X Type of the update
//! \tparam Y Type of the defect
//! \tparam l Ignored. Just there to have the same number of template arguments
//!    as other preconditioners.
//!
//! \todo Remove the reliance on CPUMatrix. We should be able to use the GPU matrix type
//!    directly.
//!
//! \note We assume X and Y are both GpuVector<real_type>, but we leave them as template
//! arguments in case of future additions.
template <class CPUMatrixT, class X, class Y, int l = 1>
class GpuDILU : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    using GPUMatrix = GpuSparseMatrix<field_type>;
    using FloatMat = GpuSparseMatrix<float>;
    using FloatVec = GpuVector<float>;

    //! \brief The matrix type the preconditioner is for.
    using matrix_type = GPUMatrix;

    
    //! \brief Constructor.
    //!
    //!  Constructor gets all parameters to operate the prec.
    explicit GpuDILU(const GPUMatrix& gpuMatrix, const CPUMatrixT& cpuMatrix, bool splitMatrix, bool tuneKernels, int mixedPrecisionScheme);

    //! \brief Prepare the preconditioner.
    //! \note Does nothing at the time being.
    void pre(X& x, Y& b) override;

    //! \brief Apply the preconditoner.
    void apply(X& v, const Y& d) override;

    //! \brief Post processing
    //! \note Does nothing at the moment
    void post(X& x) override;

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override;

    //! \brief Updates the matrix data.
    void update() final;

    //! \brief perform matrix splitting and reordering
    void reorderAndSplitMatrix(int moveThreadBlockSize);

    //! \brief Compute the diagonal of the DILU, and update the data of the reordered matrix
    void computeDiagonal(int factorizationThreadBlockSize);

    //! \brief function that will experimentally tune the thread block sizes of the important cuda kernels
    void tuneThreadBlockSizes();


    //! \returns false
    static constexpr bool shouldCallPre()
    {
        return false;
    }

    //! \returns false
    static constexpr bool shouldCallPost()
    {
        return false;
    }

    virtual bool hasPerfectUpdate() const override {
        return true;
    }


private:
    //! \brief Apply the preconditoner.
    void apply(X& v, const Y& d, int lowerSolveThreadBlockSize, int upperSolveThreadBlockSize);
    //! \brief Updates the matrix data.
    void update(int moveThreadBlockSize, int factorizationThreadBlockSize);
    //! \brief size_t describing the dimensions of the square block elements
    static constexpr const size_t blocksize_ = CPUMatrixT::block_type::cols;
    //! \brief SparseTable storing each row by level
    Opm::SparseTable<size_t> m_levelSets;
    //! \brief converts from index in reordered structure to index natural ordered structure
    std::vector<int> m_reorderedToNatural;
    //! \brief converts from index in natural ordered structure to index reordered strucutre
    std::vector<int> m_naturalToReordered;
    //! \brief The A matrix stored on the gpu, and its reordred version
    const GPUMatrix& m_gpuMatrix;
    //! \brief Stores the matrix in its entirety reordered. Optional in case splitting is used
    std::unique_ptr<GPUMatrix> m_gpuMatrixReordered;
    //! \brief If matrix splitting is enabled, then we store the lower and upper part separately
    std::unique_ptr<GPUMatrix> m_gpuMatrixReorderedLower;
    std::unique_ptr<GPUMatrix> m_gpuMatrixReorderedUpper;
    //! \brief If matrix splitting is enabled, we also store the diagonal separately
    std::unique_ptr<GpuVector<field_type>> m_gpuMatrixReorderedDiag;
    //! \brief If mixed precision is enabled, store a float matrix
    std::unique_ptr<FloatMat> m_gpuMatrixReorderedLowerFloat;
    std::unique_ptr<FloatMat> m_gpuMatrixReorderedUpperFloat;
    std::unique_ptr<FloatVec> m_gpuMatrixReorderedDiagFloat;
    std::unique_ptr<FloatVec> m_gpuDInvFloat;
    //! row conversion from natural to reordered matrix indices stored on the GPU
    GpuVector<int> m_gpuNaturalToReorder;
    //! row conversion from reordered to natural matrix indices stored on the GPU
    GpuVector<int> m_gpuReorderToNatural;
    //! \brief Stores the inverted diagonal that we use in DILU
    GpuVector<field_type> m_gpuDInv;
    //! \brief Bool storing whether or not we should store matrices in a split format
    bool m_splitMatrix;
    //! \brief Bool storing whether or not we will tune the threadblock sizes. Only used for AMD cards
    bool m_tuneThreadBlockSizes;
    //! \brief Enum describing how we should store the factorized matrix
    const MatrixStorageMPScheme m_mixedPrecisionScheme;
    //! \brief variables storing the threadblocksizes to use if using the tuned sizes and AMD cards
    //! The default value of -1 indicates that we have not calibrated and selected a value yet
    int m_upperSolveThreadBlockSize = -1;
    int m_lowerSolveThreadBlockSize = -1;
    int m_moveThreadBlockSize = -1;
    int m_DILUFactorizationThreadBlockSize = -1;

    // Graphs for Apply
    std::map<std::pair<field_type*, const field_type*>, GPUGraph> m_apply_graphs;
    std::map<std::pair<field_type*, const field_type*>, GPUGraphExec> m_executableGraphs;

    // Stream for the DILU operations on the GPU
    GPUStream m_stream{};
    // Events for synchronization with main stream
    GPUEvent m_before{};
    GPUEvent m_after{};
};
} // end namespace Opm::gpuistl

#endif
