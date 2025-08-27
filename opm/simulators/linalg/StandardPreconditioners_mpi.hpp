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

#ifndef OPM_STANDARDPRECONDITIONERS_MPI_HEADER
#define OPM_STANDARDPRECONDITIONERS_MPI_HEADER

#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/PreconditionerCPUMatrixToGPUMatrix.hpp>
#endif


namespace Opm {


template <class Smoother>
struct AMGSmootherArgsHelper
{
    static auto args(const PropertyTree& prm)
    {
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("iterations", 1);
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("relaxation", 1.0);
        return smootherArgs;
    }
};

template <class M, class V, class C>
struct AMGSmootherArgsHelper<ParallelOverlappingILU0<M, V, V, C>>
{
    static auto args(const PropertyTree& prm)
    {
        using Smoother = ParallelOverlappingILU0<M, V, V, C>;
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("iterations", 1);
        const int iluwitdh = prm.get<int>("iluwidth", 0);
        smootherArgs.setN(iluwitdh);
        const MILU_VARIANT milu = convertString2Milu(prm.get<std::string>("milutype", std::string("ilu")));
        smootherArgs.setMilu(milu);
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("relaxation", 1.0);
        return smootherArgs;
    }
};

// trailing return type with decltype used for detecting existence of setUseFixedOrder member function by overloading the setUseFixedOrder function
template <typename C>
auto setUseFixedOrder(C& criterion, bool booleanValue) -> decltype(criterion.setUseFixedOrder(booleanValue))
{
    return criterion.setUseFixedOrder(booleanValue); // Set flag to ensure that the matrices in the AMG hierarchy are constructed with deterministic indices.
}
template <typename C>
void setUseFixedOrder(C&, ...)
{
    // do nothing, since the function setUseFixedOrder does not exist yet
}

template <class Operator, class Comm, class Matrix, class Vector>
typename AMGHelper<Operator, Comm, Matrix, Vector>::Criterion
AMGHelper<Operator, Comm, Matrix, Vector>::criterion(const PropertyTree& prm)
{
    Criterion criterion(15, prm.get<int>("coarsenTarget", 1200));
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm.get<double>("alpha", 0.33));
    criterion.setBeta(prm.get<double>("beta", 1e-5));
    criterion.setMaxLevel(prm.get<int>("maxlevel", 15));
    criterion.setSkipIsolated(prm.get<bool>("skip_isolated", false));
    criterion.setNoPreSmoothSteps(prm.get<int>("pre_smooth", 1));
    criterion.setNoPostSmoothSteps(prm.get<int>("post_smooth", 1));
    criterion.setDebugLevel(prm.get<int>("verbosity", 0));
    // As the default we request to accumulate data to 1 process always as our matrix
    // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
    // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
    criterion.setAccumulate(static_cast<Dune::Amg::AccumulationMode>(prm.get<int>("accumulate", 1)));
    criterion.setProlongationDampingFactor(prm.get<double>("prolongationdamping", 1.6));
    criterion.setMaxDistance(prm.get<int>("maxdistance", 2));
    criterion.setMaxConnectivity(prm.get<int>("maxconnectivity", 15));
    criterion.setMaxAggregateSize(prm.get<int>("maxaggsize", 6));
    criterion.setMinAggregateSize(prm.get<int>("minaggsize", 4));
    setUseFixedOrder(criterion, true); // If possible, set flag to ensure that the matrices in the AMG hierarchy are constructed with deterministic indices.
    return criterion;
}

template <class Operator, class Comm, class Matrix, class Vector>
template <class Smoother>
typename AMGHelper<Operator, Comm, Matrix, Vector>::PrecPtr
AMGHelper<Operator, Comm, Matrix, Vector>::makeAmgPreconditioner(const Operator& op,
                                                                 const PropertyTree& prm,
                                                                 bool useKamg)
{
    auto crit = criterion(prm);
    auto sargs = AMGSmootherArgsHelper<Smoother>::args(prm);
    if (useKamg) {
        using Type = Dune::DummyUpdatePreconditioner<Dune::Amg::KAMG<Operator, Vector, Smoother>>;
        return std::make_shared<Type>(
            op, crit, sargs, prm.get<std::size_t>("max_krylov", 1), prm.get<double>("min_reduction", 1e-1));
    } else {
        using Type = Dune::Amg::AMGCPR<Operator, Vector, Smoother>;
        return std::make_shared<Type>(op, crit, sargs);
    }
}

template <class Operator, class Comm, typename = void> // Note: Last argument is to allow partial specialization for GPU
struct StandardPreconditioners 
{
    static void add()
    {
        using namespace Dune;
        using O = Operator;
        using C = Comm;
        using F = PreconditionerFactory<O, C>;
        using M = typename F::Matrix;
        using V = typename F::Vector;
        using P = PropertyTree;
        F::addCreator("ilu0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            return createParILU(op, prm, comm, 0);
        });
        F::addCreator("paroverilu0",
                      [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
                          return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
                      });
        F::addCreator("ilun", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
        });
        F::addCreator("duneilu", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            const bool resort = prm.get<bool>("resort", false);
            return wrapBlockPreconditioner<RebuildOnUpdatePreconditioner<Dune::SeqILU<M, V, V>>>(
                comm, op.getmat(), n, w, resort);
        });
        F::addCreator("dilu", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            DUNE_UNUSED_PARAMETER(prm);
            return wrapBlockPreconditioner<MultithreadDILU<M, V, V>>(comm, op.getmat());
        });
        F::addCreator("jac", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqJac<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("gs", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqGS<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("sor", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("ssor", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is the overlapping schwarz operator or GhostLastMatrixAdapter. This could be extended
        // later, but at this point no other operators are compatible
        // with the AMG hierarchy construction.
        if constexpr (std::is_same_v<O, Dune::OverlappingSchwarzOperator<M, V, V, C>> ||
                      std::is_same_v<O, Opm::GhostLastMatrixAdapter<M, V, V, C>>) {
            F::addCreator("amg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
                using PrecPtr = std::shared_ptr<Dune::PreconditionerWithUpdate<V, V>>;
                std::string smoother = prm.get<std::string>("smoother", "paroverilu0");
                // Make the smoother type lowercase for internal canonical representation
                std::transform(smoother.begin(), smoother.end(), smoother.begin(), ::tolower);
                // TODO: merge this with ILUn, and possibly simplify the factory to only work with ILU?
                if (smoother == "ilu0" || smoother == "paroverilu0") {
                    using Smoother = ParallelOverlappingILU0<M, V, V, C>;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    auto sargs = AMGSmootherArgsHelper<Smoother>::args(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "dilu") {
                    using SeqSmoother = Dune::MultithreadDILU<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "jac") {
                    using SeqSmoother = SeqJac<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "gs") {
                    using SeqSmoother = SeqGS<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "sor") {
                    using SeqSmoother = SeqSOR<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "ssor") {
                    using SeqSmoother = SeqSSOR<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else if (smoother == "ilun") {
                    using SeqSmoother = SeqILU<M, V, V>;
                    using Smoother = Dune::BlockPreconditioner<V, V, C, SeqSmoother>;
                    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
                    SmootherArgs sargs;
                    auto crit = AMGHelper<O, C, M, V>::criterion(prm);
                    PrecPtr prec = std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                    return prec;
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " + smoother + ".");
                }
            });
#if HAVE_HYPRE
            if constexpr (M::block_type::rows == 1 && M::block_type::cols == 1
                          && std::is_same_v<HYPRE_Real, typename V::field_type>) {
                F::addCreator(
                    "hypre", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
                        return std::make_shared<Hypre::HyprePreconditioner<M, V, V, C>>(op.getmat(), prm, comm);
                    });
            }
#endif
        }


        F::addCreator("cpr",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()> weightsCalculator,
                         std::size_t pressureIndex,
                         const C& comm) {
                          assert(weightsCalculator);
                          if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                              OPM_THROW(std::logic_error,
                                        "Pressure index out of bounds. It needs to specified for CPR");
                          }
                          using Scalar = typename V::field_type;
                          using LevelTransferPolicy = PressureTransferPolicy<O, Comm, Scalar, false>;
                          return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(
                              op, prm, weightsCalculator, pressureIndex, comm);
                      });
        F::addCreator("cprt",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()> weightsCalculator,
                         std::size_t pressureIndex,
                         const C& comm) {
                          assert(weightsCalculator);
                          if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                              OPM_THROW(std::logic_error,
                                        "Pressure index out of bounds. It needs to specified for CPR");
                          }
                          using Scalar = typename V::field_type;
                          using LevelTransferPolicy = PressureTransferPolicy<O, Comm, Scalar, true>;
                          return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(
                              op, prm, weightsCalculator, pressureIndex, comm);
                      });

        // Add CPRW only for the WellModelGhostLastMatrixAdapter, as the method requires that the
        // operator has the addWellPressureEquations() method (and a few more) it can not be combined
        // with a well-less operator such as GhostLastMatrixAdapter or OverlappingSchwarzOperator.
        // For OPM Flow this corresponds to requiring --matrix-add-well-contributions=false
        // (which is the default).
        if constexpr (std::is_same_v<O, WellModelGhostLastMatrixAdapter<M, V, V, true>>) {
            F::addCreator("cprw",
                          [](const O& op,
                             const P& prm,
                             const std::function<V()> weightsCalculator,
                             std::size_t pressureIndex,
                             const C& comm) {
                              assert(weightsCalculator);
                              if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                                  OPM_THROW(std::logic_error,
                                            "Pressure index out of bounds. It needs to specified for CPR");
                              }
                              using Scalar = typename V::field_type;
                              using LevelTransferPolicy = PressureBhpTransferPolicy<O, Comm, Scalar, false>;
                              return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(
                                  op, prm, weightsCalculator, pressureIndex, comm);
                          });
        }

#if HAVE_CUDA
        // Here we create the *wrapped* GPU preconditioners
        // meaning they will act as CPU preconditioners on the outside,
        // but copy data back and forth to the GPU as needed.

        // TODO: Make this use the GPU preconditioner factory once that is up and running.
        F::addCreator("gpuilu0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const double w = prm.get<double>("relaxation", 1.0);
            using field_type = typename V::field_type;
            using GpuILU0 = typename gpuistl::
                GpuSeqILU0<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            auto gpuILU0 = std::make_shared<GpuILU0>(op.getmat(), w);

            auto adapted = std::make_shared<gpuistl::PreconditionerAdapter<V, V, GpuILU0>>(gpuILU0);
            auto wrapped = std::make_shared<gpuistl::GpuBlockPreconditioner<V, V, Comm>>(adapted, comm);
            return wrapped;
        });

        F::addCreator("gpujac", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const double w = prm.get<double>("relaxation", 1.0);
            using field_type = typename V::field_type;
            using GpuJac =
                typename gpuistl::GpuJac<gpuistl::GpuSparseMatrix<field_type>, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, GpuJac, M>;
           
            auto gpuJac = std::make_shared<MatrixOwner>(op.getmat(), w);

            auto adapted = std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(gpuJac);
            auto wrapped = std::make_shared<gpuistl::GpuBlockPreconditioner<V, V, Comm>>(adapted, comm);
            return wrapped;
        });

        F::addCreator("gpudilu", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const bool split_matrix = prm.get<bool>("split_matrix", true);
            const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
            const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
            const bool reorder = prm.get<bool>("reorder", true);
            using field_type = typename V::field_type;
            using GpuDILU = typename gpuistl::GpuDILU<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;
            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, GpuDILU, M>;
        
            // Note: op.getmat() is passed twice, because the GpuDILU needs both the CPU and GPU matrix.
            // The first argument will be converted to a GPU matrix, and the second one is used as a CPU matrix.
            auto gpuDILU = std::make_shared<MatrixOwner>(op.getmat(), op.getmat(), split_matrix, tune_gpu_kernels, mixed_precision_scheme, reorder);

            auto adapted = std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(gpuDILU);
            auto wrapped = std::make_shared<gpuistl::GpuBlockPreconditioner<V, V, Comm>>(adapted, comm);
            return wrapped;
        });

        F::addCreator("opmgpuilu0", [](const O& op, [[maybe_unused]] const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const bool split_matrix = prm.get<bool>("split_matrix", true);
            const bool tune_gpu_kernels = prm.get<bool>("tune_gpu_kernels", true);
            const int mixed_precision_scheme = prm.get<int>("mixed_precision_scheme", 0);
            using field_type = typename V::field_type;
            using OpmGpuILU0 = typename gpuistl::OpmGpuILU0<M, gpuistl::GpuVector<field_type>, gpuistl::GpuVector<field_type>>;

            using MatrixOwner = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<gpuistl::GpuVector<field_type>, 
                gpuistl::GpuVector<field_type>, OpmGpuILU0, M>;
    
            // Note: op.getmat() is passed twice, because the OPMGPUILU0 needs both the CPU and GPU matrix.
            // The first argument will be converted to a GPU matrix, and the second one is used as a CPU matrix.
            auto gpuilu0 = std::make_shared<MatrixOwner>(op.getmat(), op.getmat(), split_matrix, tune_gpu_kernels, mixed_precision_scheme);

            auto adapted = std::make_shared<gpuistl::PreconditionerAdapter<V, V, MatrixOwner>>(gpuilu0);
            auto wrapped = std::make_shared<gpuistl::GpuBlockPreconditioner<V, V, Comm>>(adapted, comm);
            return wrapped;
        });
#endif // HAVE_CUDA
    }


    static typename PreconditionerFactory<Operator, Comm>::PrecPtr
    createParILU(const Operator& op, const PropertyTree& prm, const Comm& comm, const int ilulevel)
    {
        using F = PreconditionerFactory<Operator, Comm>;
        using M = typename F::Matrix;
        using V = typename F::Vector;

        const double w = prm.get<double>("relaxation", 1.0);
        const bool redblack = prm.get<bool>("redblack", false);
        const bool reorder_spheres = prm.get<bool>("reorder_spheres", false);
        // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
        if (ilulevel == 0) {
            const std::size_t num_interior = interiorIfGhostLast(comm);
            assert(num_interior <= op.getmat().N());
            return std::make_shared<ParallelOverlappingILU0<M, V, V, Comm>>(
                op.getmat(), comm, w, MILU_VARIANT::ILU, num_interior, redblack, reorder_spheres);
        } else {
            return std::make_shared<ParallelOverlappingILU0<M, V, V, Comm>>(
                op.getmat(), comm, ilulevel, w, MILU_VARIANT::ILU, redblack, reorder_spheres);
        }
    }

    /// Helper method to determine if the local partitioning has the
    /// K interior cells from [0, K-1] and ghost cells from [K, N-1].
    /// Returns K if true, otherwise returns N. This is motivated by
    /// usage in the ParallelOverlappingILU0 preconditioner.
    static std::size_t interiorIfGhostLast(const Comm& comm)
    {
        std::size_t interior_count = 0;
        std::size_t highest_interior_index = 0;
        const auto& is = comm.indexSet();
        for (const auto& ind : is) {
            if (Comm::OwnerSet::contains(ind.local().attribute())) {
                ++interior_count;
                highest_interior_index = std::max(highest_interior_index, ind.local().local());
            }
        }
        if (highest_interior_index + 1 == interior_count) {
            return interior_count;
        } else {
            return is.size();
        }
    }
};


} // namespace Opm

#endif // OPM_STANDARDPRECONDITIONERS_MPI_HEADER
