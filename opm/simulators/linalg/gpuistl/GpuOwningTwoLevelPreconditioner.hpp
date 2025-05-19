/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_GPUOWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
#define OPM_GPUOWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/gpuistl/GpuTwoLevelMethodCpr.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <memory>

namespace Opm::gpuistl
{

template <class OperatorType,
          class VectorType,
          class LevelTransferPolicy,
          class Communication = Dune::Amg::SequentialInformation>
class GpuOwningTwoLevelPreconditioner : public Dune::PreconditionerWithUpdate<VectorType, VectorType>
{
public:
    using MatrixType = typename OperatorType::matrix_type;
    using PrecFactory = Opm::PreconditionerFactory<OperatorType, Communication>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;
    using ScalarType = typename VectorType::field_type;

    GpuOwningTwoLevelPreconditioner(const OperatorType& linearoperator, const Opm::PropertyTree& prm,
                                 const std::function<VectorType()> weightsCalculator,
                                 std::size_t pressureIndex)
        : linear_operator_(linearoperator)
        , weightsCalculator_(weightsCalculator)
        , weights_(weightsCalculator ? weightsCalculator() : VectorType(linearoperator.getmat().N()))
        , finesmoother_(PrecFactory::create(linearoperator,
                                            prm.get_child_optional("finesmoother") ?
                                            prm.get_child("finesmoother") : Opm::PropertyTree(),
                                            std::function<VectorType()>(), pressureIndex))
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
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        //twolevel_method_.pre(x, b);
    }

    virtual void apply(VectorType& v, const VectorType& d) override
    {
        //twolevel_method_.apply(v, d);
    }

    virtual void post(VectorType& x) override
    {
        //twolevel_method_.post(x);
    }

    virtual void update() override
    {
        weights_ = weightsCalculator_();
        updateImpl();
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return linear_operator_.category();
    }

    virtual bool hasPerfectUpdate() const override
    {
        //return twolevel_method_.hasPerfectUpdate();
        return false;
    }

    // Static methods to indicate that this preconditioner needs pre() and post() calls
    // TODO: Potentially set those to false once we have a proper implementation of the preconditioner
    static constexpr bool shouldCallPre() { return true; }
    static constexpr bool shouldCallPost() { return true; }

private:
    using CoarseOperator = typename LevelTransferPolicy::CoarseOperator;

    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperator,
                                                               Dune::FlexibleSolver<CoarseOperator>,
                                                               LevelTransferPolicy>;

    using TwoLevelMethod = Dune::Amg::GpuTwoLevelMethodCpr<
        OperatorType,
        CoarseSolverPolicy,
        Dune::PreconditionerWithUpdate<VectorType, VectorType>,
        LevelTransferPolicy>;

    void updateImpl()
    {
        //twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    const OperatorType& linear_operator_;
    std::function<VectorType()> weightsCalculator_;
    VectorType weights_;
    std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>> finesmoother_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
    Opm::PropertyTree prm_;
    Communication dummy_comm_;

};

} // namespace Opm::gpuistl




#endif // OPM_OWNINGTWOLEVELPRECONDITIONERGPU_HEADER_INCLUDED
