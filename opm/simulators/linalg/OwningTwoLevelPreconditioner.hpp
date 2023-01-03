/*
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

#ifndef OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
#define OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <fstream>
#include <type_traits>


namespace Opm
{
// Circular dependency between PreconditionerFactory [which can make an OwningTwoLevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses PreconditionerFactory to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
template <class Operator, class Comm = Dune::Amg::SequentialInformation>
class PreconditionerFactory;
}

namespace Dune
{


// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <class Operator>
class FlexibleSolver;

template <typename T, typename A, int i>
std::ostream& operator<<(std::ostream& out,
                         const BlockVector< FieldVector< T, i >, A >&  vector)
{
   Dune::writeMatrixMarket(vector, out);
   return out;
}

    template <typename T, typename A, int i>
    std::istream& operator>>(std::istream& input,
                             BlockVector< FieldVector< T, i >, A >&  vector)
{
   Dune::readMatrixMarket(vector, input);
   return input;
}

/// A version of the two-level preconditioner that is:
/// - Self-contained, because it owns its policy components.
/// - Flexible, because it uses the runtime-flexible solver
///   and preconditioner factory.
template <class OperatorType,
          class VectorType,
          class LevelTransferPolicy,
          class Communication = Dune::Amg::SequentialInformation>
class OwningTwoLevelPreconditioner : public Dune::PreconditionerWithUpdate<VectorType, VectorType>
{
public:
    using MatrixType = typename OperatorType::matrix_type;
    using PrecFactory = Opm::PreconditionerFactory<OperatorType, Communication>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;

    OwningTwoLevelPreconditioner(const OperatorType& linearoperator, const Opm::PropertyTree& prm,
                                 const std::function<VectorType()> weightsCalculator,
                                 std::size_t pressureIndex)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator,
                                            prm.get_child_optional("finesmoother") ?
                                            prm.get_child("finesmoother") : Opm::PropertyTree(),
                                            std::function<VectorType()>(), pressureIndex))
        , comm_(nullptr)
        , weightsCalculator_(weightsCalculator)
        , weights_(weightsCalculator())
        , levelTransferPolicy_(dummy_comm_, weights_, prm, pressureIndex)
        , coarseSolverPolicy_(prm.get_child_optional("coarsesolver") ? prm.get_child("coarsesolver") : Opm::PropertyTree())
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           prm.get<int>("pre_smooth", 0),
                           prm.get<int>("post_smooth", 1))
        , prm_(prm)
    {
        if (prm.get<int>("verbosity", 0) > 10) {
            std::string filename = prm.get<std::string>("weights_filename", "impes_weights.txt");
            std::ofstream outfile(filename);
            if (!outfile) {
                OPM_THROW(std::ofstream::failure,
                          "Could not write weights to file " + filename + ".");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    OwningTwoLevelPreconditioner(const OperatorType& linearoperator, const Opm::PropertyTree& prm,
                                 const std::function<VectorType()> weightsCalculator,
                                 std::size_t pressureIndex, const Communication& comm)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator,
                                            prm.get_child_optional("finesmoother") ?
                                            prm.get_child("finesmoother"): Opm::PropertyTree(),
                                            std::function<VectorType()>(),
                                            comm, pressureIndex))
        , comm_(&comm)
        , weightsCalculator_(weightsCalculator)
        , weights_(weightsCalculator())
        , levelTransferPolicy_(*comm_, weights_, prm, pressureIndex)
        , coarseSolverPolicy_(prm.get_child_optional("coarsesolver") ? prm.get_child("coarsesolver") : Opm::PropertyTree())
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           prm.get<int>("pre_smooth", 0),
                           prm.get<int>("post_smooth", 1))
        , prm_(prm)
    {
        if (prm.get<int>("verbosity", 0) > 10 && comm.communicator().rank() == 0) {
            auto filename = prm.get<std::string>("weights_filename", "impes_weights.txt");
            std::ofstream outfile(filename);
            if (!outfile) {
                OPM_THROW(std::ofstream::failure,
                          "Could not write weights to file " + filename + ".");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }

    virtual void apply(VectorType& v, const VectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }

    virtual void post(VectorType& x) override
    {
        twolevel_method_.post(x);
    }

    virtual void update() override
    {
        weights_ = weightsCalculator_();
        updateImpl(comm_);
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return linear_operator_.category();
    }

private:
    using CoarseOperator = typename LevelTransferPolicy::CoarseOperator;
    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperator,
                                                               FlexibleSolver<CoarseOperator>,
                                                               LevelTransferPolicy>;
    using TwoLevelMethod
        = Dune::Amg::TwoLevelMethodCpr<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    // Handling parallel vs serial instantiation of preconditioner factory.
    template <class Comm>
    void updateImpl(const Comm*)
    {
        // Parallel case.
        auto child = prm_.get_child_optional("finesmoother");
        finesmoother_ = PrecFactory::create(linear_operator_, child ? *child : Opm::PropertyTree(), *comm_);
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    void updateImpl(const Dune::Amg::SequentialInformation*)
    {
        // Serial case.
        auto child = prm_.get_child_optional("finesmoother");
        finesmoother_ = PrecFactory::create(linear_operator_, child ? *child : Opm::PropertyTree());
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    const OperatorType& linear_operator_;
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    const Communication* comm_;
    std::function<VectorType()> weightsCalculator_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
    Opm::PropertyTree prm_;
    Communication dummy_comm_;
};

} // namespace Dune




#endif // OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
