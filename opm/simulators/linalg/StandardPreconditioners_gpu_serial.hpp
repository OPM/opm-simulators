/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_STANDARDPRECONDITIONERS_GPU_SERIAL_HEADER
#define OPM_STANDARDPRECONDITIONERS_GPU_SERIAL_HEADER


#include <dune/istl/bcrsmatrix.hh>
#include <type_traits>


namespace Opm
{


template <class Operator>
struct StandardPreconditioners<Operator,
                               Dune::Amg::SequentialInformation,
                               typename std::enable_if_t<Opm::is_gpu_operator_v<Operator>>> 
{

    using O = Operator;
    using C = Dune::Amg::SequentialInformation;
    using F = PreconditionerFactory<O, C>;
    using M = typename F::Matrix;
    using V = typename F::Vector;
    using P = PropertyTree;
    using PrecPtr = typename F::PrecPtr;

    using field_type = typename V::field_type;

    static constexpr int maxblocksize = 6;
    static void add()
    {
        F::addCreator("ilu0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using GpuILU0 =
                typename gpuistl::GpuSeqILU0<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            return std::make_shared<GpuILU0>(op.getmat(), w);
        });

        F::addCreator("jac", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using GPUJac = typename gpuistl::GpuJac<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            return std::make_shared<GPUJac>(op.getmat(), w);
        });

        // The next two (OPMILU0 and DILU) are the GPU preconditioners that need CPU matrices.
        // Since we are not storing the block size compile time for the GPU matrices 
        // **and** we need the CPU matrix for the creation, we need to
        // dispatch the creation of the preconditioner based on the block size. 
        //
        // Note that this dispatching is not needed in the future, since we will have a constructor taking GPU matrices directly.
        F::addCreator("opmilu0", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) {
            return dispatchOnBlockSize<maxblocksize>(op, [&](const auto& cpuMatrix) {
                const bool split_matrix = prm.get<bool>("split_matrix", true);
                const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
                const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
                using CPUMatrixType = std::remove_const_t<std::remove_reference_t<decltype(cpuMatrix)>>;
                using GPUILU0 =
                    typename gpuistl::OpmGpuILU0<CPUMatrixType, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;

                return std::make_shared<GPUILU0>(op.getmat(), cpuMatrix, split_matrix, tune_gpu_kernels, mixed_precision_scheme);
            });
        });

        F::addCreator("dilu", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) {
            return dispatchOnBlockSize<maxblocksize>(op, [&](const auto& cpuMatrix) {
                const bool split_matrix = prm.get<bool>("split_matrix", true);
                const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
                const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
                const bool reorder = prm.get<bool>("reorder", true);
                using CPUMatrixType = std::remove_const_t<std::remove_reference_t<decltype(cpuMatrix)>>;
                using GPUDILU =
                    typename gpuistl::GpuDILU<CPUMatrixType, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
                return std::make_shared<GPUDILU>(op.getmat(), cpuMatrix, split_matrix, tune_gpu_kernels, mixed_precision_scheme, reorder);
            });
        });

        // Only add AMG preconditioners to the factory if the operator
        // is an actual matrix operator.
        if constexpr (std::is_same_v<O, Dune::MatrixAdapter<M, V, V>>) {
#if HAVE_AMGX
            // Only add AMGX for scalar matrices
            F::addCreator("amgx", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                // Only create AMGX preconditioner for scalar matrices
                if (op.getmat().blockSize() == 1) {
                    auto prm_copy = prm;
                    prm_copy.put("setup_frequency", Opm::Parameters::Get<Opm::Parameters::CprReuseInterval>());
                    return std::make_shared<Amgx::AmgxPreconditioner<M, V, V>>(op.getmat(), prm_copy);
                } else {
                    OPM_THROW(std::logic_error, "AMGX preconditioner only works with scalar matrices (block size 1)");
                }
            });
#endif

#if HAVE_HYPRE && HYPRE_USING_CUDA || HYPRE_USING_HIP
            // Only register Hypre preconditioner for double precision
            if constexpr (std::is_same_v<HYPRE_Real, typename V::field_type>) {
                F::addCreator("hypre", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                    // Only create Hypre preconditioner for scalar matrices
                    if (op.getmat().blockSize() == 1) {
                        return std::make_shared<Hypre::HyprePreconditioner<M, V, V, Dune::Amg::SequentialInformation>>(op.getmat(), prm, Dune::Amg::SequentialInformation());
                    } else {
                        OPM_THROW(std::logic_error, "Hypre preconditioner only works with scalar matrices (block size 1).");
                    }
                });
            }
#endif

        }
    }


private:

    /**
     * This function creates a CPU matrix from the operator holding a GPU matrix.
     * 
     * This is a workaround for now since some of the GPU preconditioners need a 
     * CPU matrix for the intial setup (graph coloring). The CPU matrix is only
     * used in the constructor, **not** in the update function or the apply function.
     */
    template<class BlockType>
    static Dune::BCRSMatrix<BlockType> makeCPUMatrix(const O& op) {
        // TODO: Make this more efficient. Maybe we can simply copy the memory areas directly?
        //       Do note that this function is anyway going away when we have a GPU
        //       constructor for the preconditioners, so it is not a priority.
        const auto& gpuMatrix = op.getmat();

        const auto nonZeros = gpuMatrix.getNonZeroValues().asStdVector();
        const auto rowIndices = gpuMatrix.getRowIndices().asStdVector();
        const auto columnIndices = gpuMatrix.getColumnIndices().asStdVector();

        const auto numberOfNonZeroes = gpuMatrix.nonzeroes();
        const auto N = gpuMatrix.N();

        Dune::BCRSMatrix<BlockType> matrix(N, N, numberOfNonZeroes, Dune::BCRSMatrix<BlockType>::row_wise);
        for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {

            for (auto j = rowIndices[row.index()]; j != rowIndices[row.index() + 1]; ++j) {
                const auto columnIndex = columnIndices[j];
                row.insert(columnIndex);
            }
        }

        for (std::size_t i = 0; i < N; ++i) {
            for (auto j = rowIndices[i]; j != rowIndices[i + 1]; ++j) {
                const auto columnIndex = columnIndices[j];
                // Now it gets a bit tricky, first we need to fetch the block matrix
                BlockType blockMatrix;
                constexpr static auto rows = BlockType::rows;

                for (std::size_t k = 0; k < rows; ++k) {
                    for (std::size_t l = 0; l < rows; ++l) {
                        blockMatrix[k][l] = nonZeros[j * rows * rows + k * rows + l];
                    }
                }
                matrix[i][columnIndex] = blockMatrix;
            }
        }

        return matrix;
    }

    /**
     * This function dispatches the creation of the preconditioner based on the block size.
     * 
     * Note that this is needed since the GPU operators/matrices do not hold the block size compile time.
     * 
     * Also note that this function is not expected to be used in the future, since we will 
     * have a GPU constructor for the preconditioners (DILU and OPMILU0), hence removing the need
     * for a CPU matrix and for this function.
     */
    template<int blocksizeCompileTime, class CreateType>
    static PrecPtr dispatchOnBlockSize(const O& op, CreateType create) {
        if (op.getmat().blockSize() == blocksizeCompileTime) {
            const auto cpuMatrix = makeCPUMatrix<Dune::FieldMatrix<field_type, blocksizeCompileTime, blocksizeCompileTime>>(op);
            return create(cpuMatrix);
        } 
        if constexpr (blocksizeCompileTime > 1) {
            return dispatchOnBlockSize<blocksizeCompileTime - 1>(op, create);
        }
        else {
            throw std::runtime_error(fmt::format("Unsupported block size: {}. Max blocksize supported is {}.", op.getmat().blockSize(), maxblocksize));
        }
    }
    
};


} // namespace Opm

#endif // OPM_STANDARDPRECONDITIONERS_GPU_SERIAL_HEADER
