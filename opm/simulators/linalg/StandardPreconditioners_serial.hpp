/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_STANDARDPRECONDITIONERS_SERIAL_HPP
#define OPM_STANDARDPRECONDITIONERS_SERIAL_HPP

#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/PreconditionerCPUMatrixToGPUMatrix.hpp>
#endif

namespace Opm {


template <class Operator>
struct StandardPreconditioners<Operator, Dune::Amg::SequentialInformation, typename std::enable_if_t<!Opm::is_gpu_operator_v<Operator>>> 
{
    static void add()
    {
        using namespace Dune;
        using O = Operator;
        using C = Dune::Amg::SequentialInformation;
        using F = PreconditionerFactory<O, C>;
        using M = typename F::Matrix;
        using V = typename F::Vector;
        using P = PropertyTree;
        F::addCreator("ILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), 0, w, MILU_VARIANT::ILU);
        });
        F::addCreator("DuneILU", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            const bool resort = prm.get<bool>("resort", false);
            return getRebuildOnUpdateWrapper<Dune::SeqILU<M, V, V>>(op.getmat(), n, w, resort);
        });
        F::addCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            return std::make_shared<ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), n, w, MILU_VARIANT::ILU);
        });
        F::addCreator("ILUn", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), n, w, MILU_VARIANT::ILU);
        });
        F::addCreator("DILU", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            DUNE_UNUSED_PARAMETER(prm);
            return std::make_shared<MultithreadDILU<M, V, V>>(op.getmat());
        });
        F::addCreator("Jac", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return getDummyUpdateWrapper<SeqJac<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("GS", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return getDummyUpdateWrapper<SeqGS<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("SOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return getDummyUpdateWrapper<SeqSOR<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("SSOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return getDummyUpdateWrapper<SeqSSOR<M, V, V>>(op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is an actual matrix operator.
        if constexpr (std::is_same_v<O, Dune::MatrixAdapter<M, V, V>>) {
            F::addCreator("amg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
                    using Smoother = SeqILU<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "GS") {
                    using Smoother = SeqGS<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "DILU") {
                    using Smoother = MultithreadDILU<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "ILUn") {
                    using Smoother = SeqILU<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " + smoother + ".");
                }
            });
            F::addCreator("kamg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
                    using Smoother = SeqILU<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "GS") {
                    using Smoother = SeqGS<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "ILUn") {
                    using Smoother = SeqILU<M, V, V>;
                    return AMGHelper<O, C, M, V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " + smoother + ".");
                }
            });
            F::addCreator("famg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                if  constexpr (std::is_same_v<typename V::field_type, float>) {
                    OPM_THROW(std::logic_error, "famg requires UMFPack which is not available for floats");
                    return nullptr;
                } else {
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    Dune::Amg::Parameters parms;
                    parms.setNoPreSmoothSteps(1);
                    parms.setNoPostSmoothSteps(1);
                    return getRebuildOnUpdateWrapper<Dune::Amg::FastAMG<O, V>>(op, crit, parms);
                }
            });

#if HAVE_AMGX
            // Only add AMGX for scalar matrices
            if constexpr (M::block_type::rows == 1 && M::block_type::cols == 1) {
                F::addCreator("amgx", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                    auto prm_copy = prm;
                    prm_copy.put("setup_frequency", Opm::Parameters::Get<Opm::Parameters::CprReuseInterval>());
                    return std::make_shared<Amgx::AmgxPreconditioner<M, V, V>>(op.getmat(), prm_copy);
                });
            }
#endif

#if HAVE_HYPRE
            // Only add Hypre for scalar matrices
            if constexpr (M::block_type::rows == 1 && M::block_type::cols == 1 &&
                          std::is_same_v<HYPRE_Real, typename V::field_type>) {
                F::addCreator("hypre", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                    return std::make_shared<Hypre::HyprePreconditioner<M, V, V>>(op.getmat(), prm);
                });
            }
#endif
        }

        // Add CPRW only for the WellModelMatrixAdapter, as the method requires that the operator
        // has the addWellPressureEquations() method (and a few more) it can not be combined with
        // a well-less operator such as Dune::MatrixAdapter.  For OPM Flow this corresponds to
        // requiring --matrix-add-well-contributions=false (which is the default).
        if constexpr (std::is_same_v<O, WellModelMatrixAdapter<M, V, V>>) {
            F::addCreator(
                "cprw",
                [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                    if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                        OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                    }
                    using Scalar = typename V::field_type;
                    using LevelTransferPolicy
                        = PressureBhpTransferPolicy<O, Dune::Amg::SequentialInformation, Scalar, false>;
                    return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(
                        op, prm, weightsCalculator, pressureIndex);
                });
        }

        F::addCreator(
            "cpr",
            [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                }
                using Scalar = typename V::field_type;
                using LevelTransferPolicy = PressureTransferPolicy<O, Dune::Amg::SequentialInformation, Scalar, false>;
                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(
                    op, prm, weightsCalculator, pressureIndex);
            });
        F::addCreator(
            "cprt",
            [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                }
                using Scalar = typename V::field_type;
                using LevelTransferPolicy = PressureTransferPolicy<O, Dune::Amg::SequentialInformation, Scalar, true>;
                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(
                    op, prm, weightsCalculator, pressureIndex);
            });

#if HAVE_CUDA
        // Here we create the *wrapped* GPU preconditioners
        // meaning they will act as CPU preconditioners on the outside,
        // but copy data back and forth to the GPU as needed.

        // TODO: Make this use the GPU preconditioner factory once that is up and running.
        F::addCreator("GPUILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using field_type = typename V::field_type;
            using GpuILU0 = typename gpuistl::
                GpuSeqILU0<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            return std::make_shared<gpuistl::PreconditionerAdapter<V, V, GpuILU0>>(
                std::make_shared<GpuILU0>(op.getmat(), w));
        });

        F::addCreator("GPUILU0Float", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using block_type = typename V::block_type;
            using VTo = Dune::BlockVector<Dune::FieldVector<float, block_type::dimension>>;
            using matrix_type_to =
                typename Dune::BCRSMatrix<Dune::FieldMatrix<float, block_type::dimension, block_type::dimension>>;
            using GpuILU0 = typename gpuistl::
                GpuSeqILU0<matrix_type_to, gpuistl::GpuVector<float>, gpuistl::GpuVector<float>>;
            using Adapter = typename gpuistl::PreconditionerAdapter<VTo, VTo, GpuILU0>;
            using Converter = typename gpuistl::PreconditionerConvertFieldTypeAdapter<Adapter, M, V, V>;
            auto converted = std::make_shared<Converter>(op.getmat());
            auto adapted = std::make_shared<Adapter>(std::make_shared<GpuILU0>(converted->getConvertedMatrix(), w));
            converted->setUnderlyingPreconditioner(adapted);
            return converted;
        });

        F::addCreator("GPUJAC", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using field_type = typename V::field_type;
            using GPUJac =
                typename gpuistl::GpuJac<gpuistl::GpuSparseMatrix<field_type>,
                    gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, GPUJac, M>;
            return std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(
                std::make_shared<MatrixOwner>(op.getmat(), w));
        });

        F::addCreator("OPMGPUILU0", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) {
            const bool split_matrix = prm.get<bool>("split_matrix", true);
            const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
            const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);

            using field_type = typename V::field_type;
            using GPUILU0 = typename gpuistl::OpmGpuILU0<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, GPUILU0, M>;
            return std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(
                // Note: op.getmat() is passed twice, because the ILU0 needs both the CPU and GPU matrix. 
                // The first argument will be converted to a GPU matrix, and the second one is used as a CPU matrix.
                std::make_shared<MatrixOwner>(op.getmat(), op.getmat(), split_matrix, tune_gpu_kernels, mixed_precision_scheme));
        });


        F::addCreator("GPUDILU", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) {
            const bool split_matrix = prm.get<bool>("split_matrix", true);
            const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
            const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
            using field_type = typename V::field_type;
            using GPUDILU = typename gpuistl::GpuDILU<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, GPUDILU, M>;
            return std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(
                // Note: op.getmat() is passed twice, because the DILU needs both the CPU and GPU matrix. 
                // The first argument will be converted to a GPU matrix, and the second one is used as a CPU matrix.
                std::make_shared<MatrixOwner>(op.getmat(), op.getmat(), split_matrix, tune_gpu_kernels, mixed_precision_scheme));
        });

        F::addCreator("GPUDILUFloat", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t) {
            const bool split_matrix = prm.get<bool>("split_matrix", true);
            const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
            const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);

            using block_type = typename V::block_type;
            using VTo = Dune::BlockVector<Dune::FieldVector<float, block_type::dimension>>;
            using matrix_type_to = typename Dune::BCRSMatrix<Dune::FieldMatrix<float, block_type::dimension, block_type::dimension>>;
            using GpuDILU = typename gpuistl::GpuDILU<matrix_type_to, gpuistl::GpuVector<float>, gpuistl::GpuVector<float>>;
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<float>, 
                gpuistl::GpuVector<float>, GpuDILU, matrix_type_to>;
            using Adapter = typename gpuistl::PreconditionerAdapter<VTo, VTo, MatrixOwner>;
            using Converter = typename gpuistl::PreconditionerConvertFieldTypeAdapter<Adapter, M, V, V>;

           
            auto converted = std::make_shared<Converter>(op.getmat());
            // Note: converted->getConvertedMatrix() is passed twice, because the DILU needs both the CPU and GPU matrix. 
            // The first argument will be converted to a GPU matrix, and the second one is used as a CPU matrix.
            auto adapted = std::make_shared<Adapter>(std::make_shared<MatrixOwner>(
                converted->getConvertedMatrix(), converted->getConvertedMatrix(),
                split_matrix, tune_gpu_kernels, mixed_precision_scheme));
            converted->setUnderlyingPreconditioner(adapted);
            return converted;
        });
#endif // HAVE_CUDA
    }
};


} // namespace Opm


#endif // OPM_STANDARDPRECONDITIONERS_SERIAL_HPP
