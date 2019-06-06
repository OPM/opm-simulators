
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

#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/amgcpr.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/preconditioners.hh>

#include <boost/property_tree/ptree.hpp>

#include <map>
#include <memory>

namespace Dune
{

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
    using Creator = std::function<PrecPtr(const Operator&, const boost::property_tree::ptree&)>;
    using ParCreator = std::function<PrecPtr(const Operator&, const boost::property_tree::ptree&, const Comm&)>;

    /// Create a new serial preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const boost::property_tree::ptree& prm)
    {
        return instance().doCreate(op, prm);
    }

    /// Create a new parallel preconditioner and return a pointer to it.
    /// \param op    operator to be preconditioned.
    /// \param prm   parameters for the preconditioner, in particular its type.
    /// \param comm  communication object (typically OwnerOverlapCopyCommunication).
    /// \return      (smart) pointer to the created preconditioner.
    static PrecPtr create(const Operator& op, const boost::property_tree::ptree& prm, const Comm& comm)
    {
        return instance().doCreate(op, prm, comm);
    }

    /// Add a creator for a serial preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, Creator creator)
    {
        instance().doAddCreator(type, creator);
    }

    /// Add a creator for a parallel preconditioner to the PreconditionerFactory.
    /// After the call, the user may obtain a preconditioner by
    /// calling create() with the given type string as a parameter
    /// contained in the property_tree.
    /// \param type     the type string we want the PreconditionerFactory to
    ///                 associate with the preconditioner.
    /// \param creator  a function or lambda creating a preconditioner.
    static void addCreator(const std::string& type, ParCreator creator)
    {
        instance().doAddCreator(type, creator);
    }

private:
    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<Matrix, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    // Helpers for creation of AMG preconditioner.
    static Criterion amgCriterion(const boost::property_tree::ptree& prm)
    {
        Criterion criterion(15, prm.get<int>("coarsenTarget"));
        criterion.setDefaultValuesIsotropic(2);
        criterion.setAlpha(prm.get<double>("alpha"));
        criterion.setBeta(prm.get<double>("beta"));
        criterion.setMaxLevel(prm.get<int>("maxlevel"));
        criterion.setSkipIsolated(false);
        criterion.setDebugLevel(prm.get<int>("verbosity"));
        return criterion;
    }

    template <typename Smoother>
    static auto amgSmootherArgs(const boost::property_tree::ptree& prm)
    {
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("iterations");
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("relaxation");
        return smootherArgs;
    }

    template <class Smoother>
    static PrecPtr makeAmgPreconditioner(const Operator& op, const boost::property_tree::ptree& prm)
    {
        auto crit = amgCriterion(prm);
        auto sargs = amgSmootherArgs<Smoother>(prm);
        return std::make_shared<Dune::Amg::AMGCPR<Operator, Vector, Smoother>>(op, crit, sargs);
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
        using P = boost::property_tree::ptree;
        using C = Comm;
        doAddCreator("ILU0", [](const O& op, const P& prm, const C& comm) {
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqILU0<M, V, V>>>(comm, op.getmat(), w);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm, const C& comm) {
            const double w = prm.get<double>("relaxation");
            // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
            return wrapPreconditioner<Opm::ParallelOverlappingILU0<M, V, V, C>>(
                op.getmat(), comm, 0, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ILUn", [](const O& op, const P& prm, const C& comm) {
            const int n = prm.get<int>("ilulevel");
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqILUn<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("Jac", [](const O& op, const P& prm, const C& comm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqJac<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm, const C& comm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqGS<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm, const C& comm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm, const C& comm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapBlockPreconditioner<DummyUpdatePreconditioner<SeqSSOR<M, V, V>>>(comm, op.getmat(), n, w);
        });
        doAddCreator("amg", [](const O& op, const P& prm, const C& comm) {
            const std::string smoother = prm.get<std::string>("smoother");
            if (smoother == "ILU0") {
                using Smoother = Opm::ParallelOverlappingILU0<M, V, V, C>;
                auto crit = amgCriterion(prm);
                auto sargs = amgSmootherArgs<Smoother>(prm);
                return std::make_shared<Dune::Amg::AMGCPR<O, V, Smoother, C>>(op, crit, sargs, comm);
            } else {
                std::string msg("No such smoother: ");
                msg += smoother;
                throw std::runtime_error(msg);
            }
        });
        doAddCreator("cpr", [](const O& op, const P& prm, const C& comm) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, false, Comm>>(op, prm, comm);
        });
        doAddCreator("cprt", [](const O& op, const P& prm, const C& comm) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, true, Comm>>(op, prm, comm);
        });
    }

    // Add a useful default set of preconditioners to the factory.
    // This is the specialization for the serial case.
    void addStandardPreconditioners(const Dune::Amg::SequentialInformation*)
    {
        using namespace Dune;
        using O = Operator;
        using M = Matrix;
        using V = Vector;
        using P = boost::property_tree::ptree;
        doAddCreator("ILU0", [](const O& op, const P& prm) {
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqILU0<M, V, V>>(op.getmat(), w);
        });
        doAddCreator("ParOverILU0", [](const O& op, const P& prm) {
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<Opm::ParallelOverlappingILU0<M, V, V>>(op.getmat(), 0, w, Opm::MILU_VARIANT::ILU);
        });
        doAddCreator("ILUn", [](const O& op, const P& prm) {
            const int n = prm.get<int>("ilulevel");
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqILUn<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("Jac", [](const O& op, const P& prm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqJac<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("GS", [](const O& op, const P& prm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqGS<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SOR", [](const O& op, const P& prm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqSOR<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("SSOR", [](const O& op, const P& prm) {
            const int n = prm.get<int>("repeats");
            const double w = prm.get<double>("relaxation");
            return wrapPreconditioner<SeqSSOR<M, V, V>>(op.getmat(), n, w);
        });
        doAddCreator("amg", [](const O& op, const P& prm) {
            const std::string smoother = prm.get<std::string>("smoother");
            if (smoother == "ILU0") {
                using Smoother = SeqILU0<M, V, V>;
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
                using Smoother = SeqILUn<M, V, V>;
                return makeAmgPreconditioner<Smoother>(op, prm);
            } else {
                std::string msg("No such smoother: ");
                msg += smoother;
                throw std::runtime_error(msg);
            }
        });
        doAddCreator("famg", [](const O& op, const P& prm) {
            auto crit = amgCriterion(prm);
            Dune::Amg::Parameters parms;
            parms.setNoPreSmoothSteps(1);
            parms.setNoPostSmoothSteps(1);
            return wrapPreconditioner<Dune::Amg::FastAMG<O, V>>(op, crit, parms);
        });
        doAddCreator("cpr", [](const O& op, const P& prm) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, false>>(op, prm);
        });
        doAddCreator("cprt", [](const O& op, const P& prm) {
            return std::make_shared<OwningTwoLevelPreconditioner<O, V, true>>(op, prm);
        });
    }


    // The method that implements the singleton pattern,
    // using the Meyers singleton technique.
    static PreconditionerFactory& instance()
    {
        static PreconditionerFactory singleton;
        return singleton;
    }

    // Private constructor, to keep users from creating a PreconditionerFactory.
    PreconditionerFactory()
    {
        Comm* dummy = nullptr;
        addStandardPreconditioners(dummy);
    }

    // Actually creates the product object.
    PrecPtr doCreate(const Operator& op, const boost::property_tree::ptree& prm)
    {
        const std::string& type = prm.get<std::string>("type");
        auto it = creators_.find(type);
        if (it == creators_.end()) {
            std::ostringstream msg;
            msg << "Preconditioner type " << type << " is not registered in the factory. Available types are: ";
            for (const auto& prec : creators_) {
                msg << prec.first << ' ';
            }
            msg << std::endl;
            throw std::runtime_error(msg.str());
        }
        return it->second(op, prm);
    }

    PrecPtr doCreate(const Operator& op, const boost::property_tree::ptree& prm, const Comm& comm)
    {
        const std::string& type = prm.get<std::string>("type");
        auto it = parallel_creators_.find(type);
        if (it == parallel_creators_.end()) {
            std::ostringstream msg;
            msg << "Parallel preconditioner type " << type
                << " is not registered in the factory. Available types are: ";
            for (const auto& prec : parallel_creators_) {
                msg << prec.first << ' ';
            }
            msg << std::endl;
            throw std::runtime_error(msg.str());
        }
        return it->second(op, prm, comm);
    }

    // Actually adds the creator.
    void doAddCreator(const std::string& type, Creator c)
    {
        creators_[type] = c;
    }

    // Actually adds the creator.
    void doAddCreator(const std::string& type, ParCreator c)
    {
        parallel_creators_[type] = c;
    }

    // This map contains the whole factory, i.e. all the Creators.
    std::map<std::string, Creator> creators_;
    std::map<std::string, ParCreator> parallel_creators_;
};

} // namespace Dune

#endif // OPM_PRECONDITIONERFACTORY_HEADER
