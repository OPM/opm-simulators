
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

#ifndef OPM_PRECONDITIONERFACTORY_HEADER
#define OPM_PRECONDITIONERFACTORY_HEADER

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <dune/istl/paamg/aggregates.hh>
#include <dune/istl/paamg/matrixhierarchy.hh>

#include <cstddef>
#include <map>
#include <memory>
#include <limits>
#include <string>

#if HAVE_CUDA
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#endif

namespace Opm
{

class PropertyTree;

template <class Operator, class Comm, class Matrix, class Vector>
struct AMGHelper
{
    using PrecPtr = std::shared_ptr<Dune::PreconditionerWithUpdate<Vector, Vector>>;
    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Matrix, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    static Criterion criterion(const PropertyTree& prm);

    template <class Smoother>
    static PrecPtr makeAmgPreconditioner(const Operator& op,
                                         const PropertyTree& prm,
                                         bool useKamg = false);
};

/// This is an object factory for creating preconditioners.  The
/// user need only interact with the factory through the static
/// methods addStandardPreconditioners() and create(). In addition
/// a user can call the addCreator() static method to add further
/// preconditioners.
template <class Operator, class Comm>
class PreconditionerFactory
{
public:
    /// Linear algebra types.
    using Matrix = typename Operator::matrix_type;
    using Vector = typename Operator::domain_type; // Assuming symmetry: that domain and range types are the same.

    /// The type of pointer returned by create().
    using PrecPtr = std::shared_ptr<Dune::PreconditionerWithUpdate<Vector, Vector>>;

    /// The type of creator functions passed to addCreator().
    using Creator = std::function<PrecPtr(const Operator&, const PropertyTree&,
                                          const std::function<Vector()>&, std::size_t)>;
    using ParCreator = std::function<PrecPtr(const Operator&, const PropertyTree&,
                                             const std::function<Vector()>&, std::size_t, const Comm&)>;

    /// Create a new serial preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param weightsCalculator Calculator for weights used in CPR.
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const PropertyTree& prm,
                          const std::function<Vector()>& weightsCalculator = {},
                          std::size_t pressureIndex = std::numeric_limits<std::size_t>::max());

    /// Create a new parallel preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param comm  communication object (typically OwnerOverlapCopyCommunication).
    /// \param weightsCalculator Calculator for weights used in CPR.
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const PropertyTree& prm,
                          const std::function<Vector()>& weightsCalculator, const Comm& comm,
                          std::size_t pressureIndex = std::numeric_limits<std::size_t>::max());

    /// Create a new parallel preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param comm  communication object (typically OwnerOverlapCopyCommunication).
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const PropertyTree& prm, const Comm& comm,
                          std::size_t pressureIndex = std::numeric_limits<std::size_t>::max());

    /// Add a creator for a serial preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, Creator creator);

    /// Add a creator for a parallel preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, ParCreator creator);

    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Matrix, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;


private:

    /// Helper struct to explicitly overload amgSmootherArgs() version for
    /// ParallelOverlappingILU0, since in-class specialization is not allowed.
    template <typename X> struct Id { using Type = X; };

    template <typename Smoother>
    static auto amgSmootherArgs(const PropertyTree& prm)
    {
        return amgSmootherArgs(prm, Id<Smoother>());
    }

    template <typename Smoother>
    static auto amgSmootherArgs(const PropertyTree& prm,
                                Id<Smoother>)
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

    static auto amgSmootherArgs(const PropertyTree& prm,
                                Id<Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Comm>>)
    {
        using Smoother = Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Comm>;
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

    template <class Smoother>
    static PrecPtr makeAmgPreconditioner(const Operator& op, const PropertyTree& prm, bool useKamg = false)
    {
        auto crit = amgCriterion(prm);
        auto sargs = amgSmootherArgs<Smoother>(prm);
	if(useKamg){
	    return std::make_shared<
		Dune::DummyUpdatePreconditioner<
		    Dune::Amg::KAMG< Operator, Vector, Smoother>
		    >
		>(op, crit, sargs,
		  prm.get<size_t>("max_krylov", 1),
		  prm.get<double>("min_reduction", 1e-1)  );
	}else{
            return std::make_shared<Dune::Amg::AMGCPR<Operator, Vector, Smoother>>(op, crit, sargs);
        }
    }

    /// Helper method to determine if the local partitioning has the
    /// K interior cells from [0, K-1] and ghost cells from [K, N-1].
    /// Returns K if true, otherwise returns N. This is motivated by
    /// usage in the ParallelOverlappingILU0 preconditiner.
    template <class CommArg>
    static size_t interiorIfGhostLast(const CommArg& comm)
    {
        size_t interior_count = 0;
        size_t highest_interior_index = 0;
        const auto& is = comm.indexSet();
        for (const auto& ind : is) {
            if (CommArg::OwnerSet::contains(ind.local().attribute())) {
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

    static PrecPtr
    createParILU(const Operator& op, const PropertyTree& prm, const Comm& comm, const int ilulevel)
    {
        const double w = prm.get<double>("relaxation", 1.0);
        const bool redblack = prm.get<bool>("redblack", false);
        const bool reorder_spheres = prm.get<bool>("reorder_spheres", false);
        // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
        if (ilulevel == 0) {
            const size_t num_interior = interiorIfGhostLast(comm);
            return std::make_shared<Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Comm>>(
                op.getmat(), comm, w, Opm::MILU_VARIANT::ILU, num_interior, redblack, reorder_spheres);
        } else {
            return std::make_shared<Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Comm>>(
                op.getmat(), comm, ilulevel, w, Opm::MILU_VARIANT::ILU, redblack, reorder_spheres);
        }
    }

    // Add a useful default set of preconditioners to the factory.
    // This is the default template, used for parallel preconditioners.
    // (Serial specialization below).
    template <class CommArg>
    void addStandardPreconditioners(const CommArg*)
    {
        using namespace Dune;
        using O = Operator;
        using M = Matrix;
        using V = Vector;
        using P = PropertyTree;
        using C = Comm;
        doAddCreator("ILU0", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            return createParILU(op, prm, comm, 0);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
        });
        doAddCreator("ILUn", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            return createParILU(op, prm, comm, prm.get<int>("ilulevel", 0));
        });
        doAddCreator("Jac", [](const O& op, const P& prm, const std::function<Vector()>&,
                               std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqJac<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqGS<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is the overlapping schwarz operator. This could be extended
        // later, but at this point no other operators are compatible
        // with the AMG hierarchy construction.
        if constexpr (std::is_same_v<O, Dune::OverlappingSchwarzOperator<M, V, V, C>>) {
            doAddCreator("amg", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t, const C& comm) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
                    using Smoother = Opm::ParallelOverlappingILU0<M, V, V, C>;
                    auto crit = amgCriterion(prm);
                    auto sargs = amgSmootherArgs<Smoother>(prm);
                    return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother << ".");
                }
            });
        }

        doAddCreator("cpr", [](const O& op, const P& prm, const std::function<Vector()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
            assert(weightsCalculator);
            if (pressureIndex == std::numeric_limits<std::size_t>::max())
            {
                OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
            }
            using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Comm, false>;
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(op, prm, weightsCalculator, pressureIndex, comm);
        });
        doAddCreator("cprt", [](const O& op, const P& prm, const std::function<Vector()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
            assert(weightsCalculator);
            if (pressureIndex == std::numeric_limits<std::size_t>::max())
            {
                OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
            }
            using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Comm, true>;
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(op, prm, weightsCalculator, pressureIndex, comm);
        });

        if constexpr (std::is_same_v<O, WellModelGhostLastMatrixAdapter<M, V, V, true>>) {
            doAddCreator("cprw",
                         [](const O& op, const P& prm, const std::function<Vector()> weightsCalculator, std::size_t pressureIndex, const C& comm) {
                             assert(weightsCalculator);
                             if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                                 OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                             }
                             using LevelTransferPolicy = Opm::PressureBhpTransferPolicy<O, Comm, false>;
                             return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy, Comm>>(
                                 op, prm, weightsCalculator, pressureIndex, comm);
                         });
        }
    }

    // Add a useful default set of preconditioners to the factory.
    // This is the specialization for the serial case.
    void addStandardPreconditioners(const Dune::Amg::SequentialInformation*)
    {
        using namespace Dune;
        using O = Operator;
        using M = Matrix;
        using V = Vector;
        using P = PropertyTree;
        doAddCreator("ILU0", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), 0, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            const int n = prm.get<int>("ilulevel", 0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ILUn", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const int n = prm.get<int>("ilulevel", 0);
            const double w = prm.get<double>("relaxation", 1.0);
            return std::make_shared<Opm::ParallelOverlappingILU0<M, V, V>>(
                op.getmat(), n, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("Jac", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqJac<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqGS<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSOR<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const int n = prm.get<int>("repeats", 1);
            const double w = prm.get<double>("relaxation", 1.0);
            return wrapPreconditioner<SeqSSOR<M, V, V>>(op.getmat(), n, w);
        });

        // Only add AMG preconditioners to the factory if the operator
        // is an actual matrix operator.
        if constexpr (std::is_same_v<O, Dune::MatrixAdapter<M, V, V>>) {
            doAddCreator("amg", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
#else
                    using Smoother = SeqILU0<M, V, V>;
#endif
                    return makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm);
                } else if (smoother == "ILUn") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
#else
                            using Smoother = SeqILUn<M, V, V>;
#endif
                    return makeAmgPreconditioner<Smoother>(op, prm);
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother << ".");
                }
            });
            doAddCreator("kamg", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
                const std::string smoother = prm.get<std::string>("smoother", "ParOverILU0");
                if (smoother == "ILU0" || smoother == "ParOverILU0") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
#else
                        using Smoother = SeqILU0<M, V, V>;
#endif
                    return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "Jac") {
                    using Smoother = SeqJac<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SOR") {
                    using Smoother = SeqSOR<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm, true);
                    // } else if (smoother == "GS") {
                    //     using Smoother = SeqGS<M, V, V>;
                    //     return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "SSOR") {
                    using Smoother = SeqSSOR<M, V, V>;
                    return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else if (smoother == "ILUn") {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
                    using Smoother = SeqILU<M, V, V>;
#else
                        using Smoother = SeqILUn<M, V, V>;
#endif
                    return makeAmgPreconditioner<Smoother>(op, prm, true);
                } else {
                    OPM_THROW(std::invalid_argument, "Properties: No smoother with name " << smoother << ".");
                }
            });
            doAddCreator("famg", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
                auto crit = amgCriterion(prm);
                Dune::Amg::Parameters parms;
                parms.setNoPreSmoothSteps(1);
                parms.setNoPostSmoothSteps(1);
                return wrapPreconditioner<Dune::Amg::FastAMG<O, V>>(op, crit, parms);
            });
        }
        if constexpr (std::is_same_v<O, WellModelMatrixAdapter<M, V, V, false>>) {
            doAddCreator("cprw", [](const O& op, const P& prm, const std::function<Vector()>& weightsCalculator, std::size_t pressureIndex) {
                if (pressureIndex == std::numeric_limits<std::size_t>::max()) {
                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                }
                using LevelTransferPolicy = Opm::PressureBhpTransferPolicy<O, Dune::Amg::SequentialInformation, false>;
                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
            });
            }

        doAddCreator("cpr", [](const O& op, const P& prm, const std::function<Vector()>& weightsCalculator, std::size_t pressureIndex) {
                                if (pressureIndex == std::numeric_limits<std::size_t>::max())
                                {
                                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                                }
                                using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Dune::Amg::SequentialInformation, false>;
                                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });
        doAddCreator("cprt", [](const O& op, const P& prm, const std::function<Vector()>& weightsCalculator, std::size_t pressureIndex) {
                                if (pressureIndex == std::numeric_limits<std::size_t>::max())
                                {
                                    OPM_THROW(std::logic_error, "Pressure index out of bounds. It needs to specified for CPR");
                                }
                                using LevelTransferPolicy = Opm::PressureTransferPolicy<O, Dune::Amg::SequentialInformation, true>;
                                return std::make_shared<OwningTwoLevelPreconditioner<O, V, LevelTransferPolicy>>(op, prm, weightsCalculator, pressureIndex);
        });

        #if HAVE_CUDA
            doAddCreator("CUILU0", [](const O& op, const P& prm, const std::function<Vector()>&, std::size_t) {
            const double w = prm.get<double>("relaxation", 1.0);
            using field_type = typename V::field_type;
            using CuILU0 = typename Opm::cuistl::CuSeqILU0<M, Opm::cuistl::CuVector<field_type>, Opm::cuistl::CuVector<field_type>>;
            //return std::make_shared<Opm::cuistl::PreconditionerAdapter<CuILU0, M, V, V>>(
            //    std::make_shared<CuILU0>(op.getmat(), w));
            wrapPreconditioner<Opm::cuistl::PreconditionerAdapter<CuILU0, M, V, V>>(std::make_shared<CuILU0>(op.getmat(), w));
        });
        #endif
    }


    // The method that implements the singleton pattern,
    // using the Meyers singleton technique.
    static PreconditionerFactory& instance();

    // Private constructor, to keep users from creating a PreconditionerFactory.
    PreconditionerFactory();

    // Actually creates the product object.
    PrecPtr doCreate(const Operator& op, const PropertyTree& prm,
                     const std::function<Vector()> weightsCalculator,
                     std::size_t pressureIndex);

    PrecPtr doCreate(const Operator& op, const PropertyTree& prm,
                     const std::function<Vector()> weightsCalculator,
                     std::size_t pressureIndex, const Comm& comm);

    // Actually adds the creator.
    void doAddCreator(const std::string& type, Creator c);

    // Actually adds the creator.
    void doAddCreator(const std::string& type, ParCreator c);

    // This map contains the whole factory, i.e. all the Creators.
    std::map<std::string, Creator> creators_;
    std::map<std::string, ParCreator> parallel_creators_;
    bool defAdded_= false; //!< True if defaults creators have been added
};

} // namespace Dune

#endif // OPM_PRECONDITIONERFACTORY_HEADER
