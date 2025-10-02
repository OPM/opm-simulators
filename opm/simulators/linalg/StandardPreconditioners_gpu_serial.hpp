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

#include <opm/simulators/linalg/gpuistl/detail/gpu_preconditioner_utils.hpp>

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
        F::addCreator("opmilu0", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) -> PrecPtr {
            return op.getmat().dispatchOnBlocksize([&](auto blockSizeVal) -> PrecPtr {
                constexpr int blockSize = decltype(blockSizeVal)::value;
                const auto cpuMatrix = gpuistl::detail::makeCPUMatrix<O, Dune::FieldMatrix<field_type, blockSize, blockSize>>(op);
                const bool split_matrix = prm.get<bool>("split_matrix", true);
                const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
                const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
                using CPUMatrixType = std::remove_const_t<std::remove_reference_t<decltype(cpuMatrix)>>;
                using GPUILU0 =
                    typename gpuistl::OpmGpuILU0<CPUMatrixType, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;

                return std::make_shared<GPUILU0>(op.getmat(), cpuMatrix, split_matrix, tune_gpu_kernels, mixed_precision_scheme);
            });
        });

        F::addCreator("dilu", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) -> PrecPtr {
            return op.getmat().dispatchOnBlocksize([&](auto blockSizeVal) -> PrecPtr {
                constexpr int blockSize = decltype(blockSizeVal)::value;
                const auto cpuMatrix = gpuistl::detail::makeCPUMatrix<O, Dune::FieldMatrix<field_type, blockSize, blockSize>>(op);
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

        F::addCreator("cpr", [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
            if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
            }
            using Scalar = typename V::field_type;
            using GpuVector = gpuistl::GpuVector<Scalar>;
            using LevelTransferPolicy = Opm::gpuistl::GpuPressureTransferPolicy<O, Dune::Amg::SequentialInformation, Scalar, false>;
            return std::make_shared<Dune::OwningTwoLevelPreconditioner<O, GpuVector, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });

        F::addCreator("cprt", [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
            if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
            }
            using Scalar = typename V::field_type;
            using GpuVector = gpuistl::GpuVector<Scalar>;
            using LevelTransferPolicy = Opm::gpuistl::GpuPressureTransferPolicy<O, Dune::Amg::SequentialInformation, Scalar, true>;
            return std::make_shared<Dune::OwningTwoLevelPreconditioner<O, GpuVector, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });
    }


};


} // namespace Opm

#endif // OPM_STANDARDPRECONDITIONERS_GPU_SERIAL_HEADER
