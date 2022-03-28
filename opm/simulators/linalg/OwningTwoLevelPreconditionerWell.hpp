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

#pragma once
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <type_traits>


namespace Opm
{
// Circular dependency between PreconditionerFactory [which can make an OwningTwoLevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses PreconditionerFactory to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
//template <class Operator, class Comm = Dune::Amg::SequentialInformation>
//class PreconditionerFactory;
}

namespace Dune
{


// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <class Operator>
class FlexibleSolver;

// template <typename T, typename A, int i>
// std::ostream& operator<<(std::ostream& out,
//                          const BlockVector< FieldVector< T, i >, A >&  vector)
// {
//    Dune::writeMatrixMarket(vector, out);
//    return out;
// }

//     template <typename T, typename A, int i>
//     std::istream& operator>>(std::istream& input,
//                              BlockVector< FieldVector< T, i >, A >&  vector)
// {
//    Dune::readMatrixMarket(vector, input);
//    return input;
// }

/// A version of the two-level preconditioner that is:
/// - Self-contained, because it owns its policy components.
/// - Flexible, because it uses the runtime-flexible solver
///   and preconditioner factory.
template <class OperatorType,
          class VectorType,
          bool transpose = false,
          class Communication = Dune::Amg::SequentialInformation>
class OwningTwoLevelPreconditionerWell : public Dune::PreconditionerWithUpdate<VectorType, VectorType>
{
public:
    using pt = boost::property_tree::ptree;
    using MatrixType = typename OperatorType::matrix_type;
    using PrecFactory = Opm::PreconditionerFactory<OperatorType, Communication>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;

    OwningTwoLevelPreconditionerWell(const OperatorType& linearoperator, const pt& prm,
                                 const std::function<VectorType()> weightsCalculator)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator,
                                            prm.get_child_optional("finesmoother")?
                                            prm.get_child("finesmoother"): pt()))
        , comm_(nullptr)
        , weightsCalculator_(weightsCalculator)
        , weights_(weightsCalculator())
        , levelTransferPolicy_(dummy_comm_, weights_, prm)//.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child_optional("coarsesolver")? prm.get_child("coarsesolver") : pt())
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           prm.get<int>("pre_smooth", transpose? 1 : 0),
                           prm.get<int>("post_smooth", transpose? 0 : 1))
        , prm_(prm)
    {
        if (prm.get<int>("verbosity", 0) > 10) {
            std::string filename = prm.get<std::string>("weights_filename", "impes_weights.txt");
            std::ofstream outfile(filename);
            if (!outfile) {
                OPM_THROW(std::ofstream::failure, "Could not write weights to file " << filename << ".");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    OwningTwoLevelPreconditionerWell(const OperatorType& linearoperator,
                                     const pt& prm,
                                     const std::function<VectorType()> weightsCalculator, const Communication& comm)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator,
                                            prm.get_child_optional("finesmoother")?
                                            prm.get_child("finesmoother"): pt(), comm))
        , comm_(&comm)
        , weightsCalculator_(weightsCalculator)
        , weights_(weightsCalculator())
        , levelTransferPolicy_(*comm_, weights_, prm)//.get<int>("pressure_var_index", 1))
        , coarseSolverPolicy_(prm.get_child_optional("coarsesolver")? prm.get_child("coarsesolver") : pt())
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           prm.get<int>("pre_smooth", transpose? 1 : 0),
                           prm.get<int>("post_smooth", transpose? 0 : 1))
        , prm_(prm)
    {
        if (prm.get<int>("verbosity", 0) > 10 && comm.communicator().rank() == 0) {
            auto filename = prm.get<std::string>("weights_filename", "impes_weights.txt");
            std::ofstream outfile(filename);
            if (!outfile) {
                OPM_THROW(std::ofstream::failure, "Could not write weights to file " << filename << ".");
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
    using PressureMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using PressureVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using SeqCoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using ParCoarseOperatorType
        = Dune::OverlappingSchwarzOperator<PressureMatrixType, PressureVectorType, PressureVectorType, Communication>;
    using CoarseOperatorType = std::conditional_t<std::is_same<Communication, Dune::Amg::SequentialInformation>::value,
                                                  SeqCoarseOperatorType,
                                                  ParCoarseOperatorType>;
    using LevelTransferPolicy = Opm::PressureBhpTransferPolicy<OperatorType, CoarseOperatorType, Communication, transpose>;
    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperatorType,
                                                               FlexibleSolver<CoarseOperatorType>,
                                                               LevelTransferPolicy>;
    using TwoLevelMethod
        = Dune::Amg::TwoLevelMethodCpr<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    // Handling parallel vs serial instantiation of preconditioner factory.
    template <class Comm>
    void updateImpl(const Comm*)
    {
        // Parallel case.
        // using ParOperatorType = Dune::OverlappingSchwarzOperator<MatrixType, VectorType, VectorType, Comm>;
        // auto op_prec = std::make_shared<ParOperatorType>(linear_operator_.getmat(), *comm_);
        auto child = prm_.get_child_optional("finesmoother");
        finesmoother_ = PrecFactory::create(linear_operator_, child ? *child : pt(), *comm_);
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
        // linearoperator_for_precond_ = op_prec;
    }

    void updateImpl(const Dune::Amg::SequentialInformation*)
    {
        // Serial case.
        // using SeqOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
        // auto op_prec = std::make_shared<SeqOperatorType>(linear_operator_.getmat());
        auto child = prm_.get_child_optional("finesmoother");
        finesmoother_ = PrecFactory::create(linear_operator_, child ? *child : pt());
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
        // linearoperator_for_precond_ = op_prec;
    }

    const OperatorType& linear_operator_;
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    const Communication* comm_;
    std::function<VectorType()> weightsCalculator_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
    boost::property_tree::ptree prm_;
    Communication dummy_comm_;
    //std::shared_ptr<AbstractOperatorType> linearoperator_for_precond_;
};

} // namespace Dune





