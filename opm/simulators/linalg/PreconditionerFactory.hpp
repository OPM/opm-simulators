
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
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <dune/istl/paamg/aggregates.hh>
#include <dune/istl/paamg/matrixhierarchy.hh>

#include <cstddef>
#include <map>
#include <memory>
#include <limits>
#include <string>

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
