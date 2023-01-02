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

#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/amgcpr.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/fastamg.hh>

namespace Opm {

template<class Smoother>
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

template<class M, class V, class C>
struct AMGSmootherArgsHelper<Opm::ParallelOverlappingILU0<M,V,V,C>>
{
    static auto args(const PropertyTree& prm)
    {
        using Smoother = Opm::ParallelOverlappingILU0<M, V, V, C>;
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

template <class Operator, class Comm, class Matrix, class Vector>
typename AMGHelper<Operator, Comm, Matrix, Vector>::Criterion
AMGHelper<Operator,Comm,Matrix,Vector>::criterion(const PropertyTree& prm)
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
    return criterion;
}

template <class Operator, class Comm, class Matrix, class Vector>
template <class Smoother>
typename AMGHelper<Operator, Comm, Matrix, Vector>::PrecPtr
AMGHelper<Operator,Comm,Matrix,Vector>::
makeAmgPreconditioner(const Operator& op,
                      const PropertyTree& prm,
                      bool useKamg)
{
    auto crit = criterion(prm);
    auto sargs = AMGSmootherArgsHelper<Smoother>::args(prm);
    if (useKamg) {
        using Type = Dune::DummyUpdatePreconditioner<Dune::Amg::KAMG<Operator, Vector, Smoother>>;
        return std::make_shared<Type>(op, crit, sargs,
                                      prm.get<size_t>("max_krylov", 1),
                                      prm.get<double>("min_reduction", 1e-1));
    } else {
        using Type = Dune::Amg::AMGCPR<Operator, Vector, Smoother>;
        return std::make_shared<Type>(op, crit, sargs);
    }
}

template<class Operator, class Comm>
struct StandardPreconditioners
{
    static void add()
    {
        using namespace Dune;
        using O = Operator;
        using C = Comm;
        using F = PreconditionerFactory<O,C>;
        using M = typename F::Matrix;
        using V = typename F::Vector;
        using P = PropertyTree;
        F::addCreator("ILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          return createParILU(op, prm, comm, 0);
        });
        F::addCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
        });
        F::addCreator("ILUn", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
        });
        F::addCreator("Jac", [](const O& op, const P& prm, const std::function<V()>&,
                     std::size_t, const C& comm) {
          const int n = prm.get<int>("repeats", 1);
          const double w = prm.get<double>("relaxation", 1.0);
          return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqJac<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("GS", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          const int n = prm.get<int>("repeats", 1);
          const double w = prm.get<double>("relaxation", 1.0);
          return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqGS<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("SOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          const int n = prm.get<int>("repeats", 1);
          const double w = prm.get<double>("relaxation", 1.0);
          return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        F::addCreator("SSOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
          const int n = prm.get<int>("repeats", 1);
          const double w = prm.get<double>("relaxation", 1.0);
          return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is the overlapping schwarz operator. This could be extended
        // later, but at this point no other operators are compatible
        // with the AMG hierarchy construction.
        if constexpr (std::is_same_v<O, Dune::OverlappingSchwarzOperator<M, V, V, C>>) {
          F::addCreator("amg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t, const C& comm) {
            const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
            if (smoother == "ILU0" || smoother == "ParOverILU0") {
              using Smoother = Opm::ParallelOverlappingILU0<M, V, V, C>;
              auto crit = AMGHelper<O,C,M,V>::criterion(prm);
              auto sargs = AMGSmootherArgsHelper<Smoother>::args(prm);
              return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
            } else {
              OPM_THROW(std::invalid_argument, "Properties: No smoother with name " + smoother + ".");
            }
          });
        }

        F::addCreator("cpr", [](const O& op, const P& prm, const std::function<V()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
          assert(weightsCalculator);
          if (pressureIndex == std::numeric_limits<std::size_t>::max())
          {
            OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
          }
          using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Comm, false>;
          return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(op, prm, weightsCalculator, pressureIndex, comm);
        });
        F::addCreator("cprt", [](const O& op, const P& prm, const std::function<V()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
          assert(weightsCalculator);
          if (pressureIndex == std::numeric_limits<std::size_t>::max())
          {
            OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
          }
          using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Comm, true>;
          return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(op, prm, weightsCalculator, pressureIndex, comm);
        });

        if constexpr (std::is_same_v<O, WellModelGhostLastMatrixAdapter<M, V, V, true>>) {
          F::addCreator("cprw",
                       [](const O& op, const P& prm, const std::function<V()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
            assert(weightsCalculator);
            if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
              OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
            }
            using LevelTransferPolicy = Opm::PressureBhpTransferPolicy<O, Comm, false>;
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(op, prm, weightsCalculator, pressureIndex, comm);
          });
        }
    }


    static typename PreconditionerFactory<Operator,Comm>::PrecPtr
    createParILU(const Operator& op, const PropertyTree& prm, const Comm& comm, const int ilulevel)
    {
        using F = PreconditionerFactory<Operator,Comm>;
        using M = typename F::Matrix;
        using V = typename F::Vector;

        const double w = prm.get<double>("relaxation", 1.0);
        const bool redblack = prm.get<bool>("redblack", false);
        const bool reorder_spheres = prm.get<bool>("reorder_spheres", false);
        // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
        if (ilulevel == 0) {
            const size_t num_interior = interiorIfGhostLast(comm);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, Comm>>(
                op.getmat(), comm, w, Opm::MILU_VARIANT::ILU, num_interior, redblack, reorder_spheres);
        } else {
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, Comm>>(
                op.getmat(), comm, ilulevel, w, Opm::MILU_VARIANT::ILU, redblack, reorder_spheres);
        }
    }

    /// Helper method to determine if the local partitioning has the
    /// K interior cells from [0, K-1] and ghost cells from [K, N-1].
    /// Returns K if true, otherwise returns N. This is motivated by
    /// usage in the ParallelOverlappingILU0 preconditioner.
    static size_t interiorIfGhostLast(const Comm& comm)
    {
        size_t interior_count = 0;
        size_t highest_interior_index = 0;
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

template<class Operator>
struct StandardPreconditioners<Operator,Dune::Amg::SequentialInformation>
{
    static void add()
    {
        using namespace Dune;
        using O = Operator;
        using C = Dune::Amg::SequentialInformation;
        using F = PreconditionerFactory<O,C>;
        using M = typename F::Matrix;
        using V = typename F::Vector;
        using P = PropertyTree;
        F::addCreator("ILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), 0, w, Opm::MILU_VARIANT::ILU);
        });
        F::addCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        F::addCreator("ILUn", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        F::addCreator("Jac", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqJac<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("GS", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqGS<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("SOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSOR<M, V, V>>(op.getmat(), n, w);
        });
        F::addCreator("SSOR", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSSOR<M, V, V>>(op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is an actual matrix operator.
        if constexpr (std::is_same_v<O, Dune::MatrixAdapter<M, V, V>>) {
            F::addCreator("amg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
    #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
    #else
                    using Smoother = SeqILU0<M, V, V>;
    #endif
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "ILUn") {
    #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
    #else
                    using Smoother = SeqILUn<M, V, V>;
    #endif
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm);
                } else {
                    OPM_THROW(std::invalid_argument,
                              "Properties: No smoother with name " + smoother + ".");
                }
            });
            F::addCreator("kamg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
    #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
    #else
                        using Smoother = SeqILU0<M, V, V>;
    #endif
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                    // } else if (smoother == "GS") {
                    //     using Smoother = SeqGS<M, V, V>;
                    //     return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "ILUn") {
    #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
    #else
                    using Smoother = SeqILUn<M, V, V>;
    #endif
                    return AMGHelper<O,C,M,V>::template makeAmgPreconditioner<Smoother>(op, prm, true);
                } else {
                    OPM_THROW(std::invalid_argument,
                              "Properties: No smoother with name " + smoother + ".");
                }
            });
            F::addCreator("famg", [](const O& op, const P& prm, const std::function<V()>&, std::size_t) {
                auto crit = AMGHelper<O,C,M,V>::criterion(prm);
                Dune::Amg::Parameters parms;
                parms.setNoPreSmoothSteps(1);
                parms.setNoPostSmoothSteps(1);
                return wrapPreconditioner<Dune::Amg::FastAMG<O, V>>(op, crit, parms);
            });
        }
        if constexpr (std::is_same_v<O, WellModelMatrixAdapter<M, V, V, false>>) {
            F::addCreator("cprw", [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                }
                using LevelTransferPolicy = Opm::PressureBhpTransferPolicy<O, Dune::Amg::SequentialInformation, false>;
                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
            });
            }

        F::addCreator("cpr", [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                                if (pressureIndex == std::numeric_limits<std::size_t>::max())
                                {
                                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                                }
                                using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Dune::Amg::SequentialInformation, false>;
                                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });
        F::addCreator("cprt", [](const O& op, const P& prm, const std::function<V()>& weightsCalculator, std::size_t pressureIndex) {
                                if (pressureIndex == std::numeric_limits<std::size_t>::max())
                                {
                                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                                }
                                using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Dune::Amg::SequentialInformation, true>;
                                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });
    }
};

template <class Operator, class Comm>
PreconditionerFactory<Operator,Comm>::PreconditionerFactory()
{
}


template <class Operator, class Comm>
PreconditionerFactory<Operator,Comm>&
PreconditionerFactory<Operator,Comm>::instance()
{
    static PreconditionerFactory singleton;
    return singleton;
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator,Comm>::PrecPtr
PreconditionerFactory<Operator,Comm>::
doCreate(const Operator& op, const PropertyTree& prm,
         const std::function<Vector()> weightsCalculator,
         std::size_t pressureIndex)
{
    if (!defAdded_) {
      StandardPreconditioners<Operator,Comm>::add();
      defAdded_ = true;
    }
    const std::string& type = prm.get<std::string>("type", "ParOverILU0");
    auto it = creators_.find(type);
    if (it == creators_.end()) {
        std::ostringstream msg;
        msg << "Preconditioner type " << type << " is not registered in the factory. Available types are: ";
        for (const auto& prec : creators_) {
            msg << prec.first << ' ';
        }
        msg << std::endl;
        OPM_THROW(std::invalid_argument, msg.str());
    }
    return it->second(op, prm, weightsCalculator, pressureIndex);
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator,Comm>::PrecPtr
PreconditionerFactory<Operator,Comm>::
doCreate(const Operator& op, const PropertyTree& prm,
         const std::function<Vector()> weightsCalculator,
         std::size_t pressureIndex, const Comm& comm)
{
    if (!defAdded_) {
        StandardPreconditioners<Operator,Comm>::add();
        defAdded_ = true;
    }
    const std::string& type = prm.get<std::string>("type", "ParOverILU0");
    auto it = parallel_creators_.find(type);
    if (it == parallel_creators_.end()) {
        std::ostringstream msg;
        msg << "Parallel preconditioner type " << type
            << " is not registered in the factory. Available types are: ";
        for (const auto& prec : parallel_creators_) {
            msg << prec.first << ' ';
        }
        msg << std::endl;
        OPM_THROW(std::invalid_argument, msg.str());
    }
    return it->second(op, prm, weightsCalculator, pressureIndex, comm);
}

template <class Operator, class Comm>
void PreconditionerFactory<Operator,Comm>::
doAddCreator(const std::string& type, Creator c)
{
    creators_[type] = c;
}

template <class Operator, class Comm>
void PreconditionerFactory<Operator,Comm>::
doAddCreator(const std::string& type, ParCreator c)
{
    parallel_creators_[type] = c;
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator,Comm>::PrecPtr
PreconditionerFactory<Operator,Comm>::
create(const Operator& op, const PropertyTree& prm,
       const std::function<Vector()>& weightsCalculator,
       std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, weightsCalculator, pressureIndex);
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator,Comm>::PrecPtr
PreconditionerFactory<Operator,Comm>::
create(const Operator& op, const PropertyTree& prm,
       const std::function<Vector()>& weightsCalculator, const Comm& comm,
       std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, weightsCalculator, pressureIndex, comm);
}


template <class Operator, class Comm>
typename PreconditionerFactory<Operator,Comm>::PrecPtr
PreconditionerFactory<Operator,Comm>::
create(const Operator& op, const PropertyTree& prm, const Comm& comm,
       std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, std::function<Vector()>(), pressureIndex, comm);
}

template <class Operator, class Comm>
void PreconditionerFactory<Operator,Comm>::
addCreator(const std::string& type, Creator creator)
{
    instance().doAddCreator(type, creator);
}

template <class Operator, class Comm>
void PreconditionerFactory<Operator,Comm>::
addCreator(const std::string& type, ParCreator creator)
{
    instance().doAddCreator(type, creator);
}

using CommSeq = Dune::Amg::SequentialInformation;

template<int Dim>
using OpFSeq = Dune::MatrixAdapter<Dune::BCRSMatrix<Dune::FieldMatrix<double,Dim,Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<double,Dim>>>;
template<int Dim>
using OpBSeq = Dune::MatrixAdapter<Dune::BCRSMatrix<Opm::MatrixBlock<double,Dim,Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<double,Dim>>>;

template<int Dim, bool overlap>
using OpW = WellModelMatrixAdapter<Dune::BCRSMatrix<MatrixBlock<double,Dim,Dim>>,
                                      Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                      Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                      overlap>;

template<int Dim, bool overlap>
using OpWG = WellModelGhostLastMatrixAdapter<Dune::BCRSMatrix<MatrixBlock<double,Dim,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                             overlap>;

#if HAVE_MPI
using CommPar = Dune::OwnerOverlapCopyCommunication<int,int>;

template<int Dim>
using OpFPar = Dune::OverlappingSchwarzOperator<Dune::BCRSMatrix<Dune::FieldMatrix<double,Dim,Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                                CommPar>;

template<int Dim>
using OpBPar = Dune::OverlappingSchwarzOperator<Dune::BCRSMatrix<MatrixBlock<double,Dim,Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<double,Dim>>,
                                                CommPar>;

#define INSTANCE_PF_PAR(Dim) \
    template class PreconditionerFactory<OpBSeq<Dim>,CommPar>; \
    template class PreconditionerFactory<OpFPar<Dim>,CommPar>; \
    template class PreconditionerFactory<OpBPar<Dim>,CommPar>; \
    template class PreconditionerFactory<OpW<Dim,false>,CommPar>; \
    template class PreconditionerFactory<OpWG<Dim,true>,CommPar>; \
    template class PreconditionerFactory<OpBPar<Dim>,CommSeq>;
#endif

#define INSTANCE_PF_SEQ(Dim) \
    template class PreconditionerFactory<OpFSeq<Dim>,CommSeq>; \
    template class PreconditionerFactory<OpBSeq<Dim>,CommSeq>; \
    template class PreconditionerFactory<OpW<Dim,false>,CommSeq>; \
    template class PreconditionerFactory<OpWG<Dim,true>,CommSeq>;

#if HAVE_MPI
#define INSTANCE_PF(Dim) \
    INSTANCE_PF_PAR(Dim) \
    INSTANCE_PF_SEQ(Dim)
#else
#define INSTANCE_PF(Dim) \
    INSTANCE_PF_SEQ(Dim)
#endif
}
