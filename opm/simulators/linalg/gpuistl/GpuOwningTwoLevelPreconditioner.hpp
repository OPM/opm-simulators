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

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>
#include <opm/simulators/linalg/gpuistl/GpuTwoLevelMethodCpr.hpp>
#include <opm/simulators/linalg/gpuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/gpuistl/PreconditionerHolder.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <fstream>
#include <memory>


namespace Opm
{
// Circular dependency between PreconditionerFactory [which can make an OwningTwoLevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses PreconditionerFactory to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
template <class Operator, class Comm>
class PreconditionerFactory;
}

namespace Dune
{

// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <class Operator>
class FlexibleSolver;

template <class OperatorType,
          class GpuVectorType,
          class LevelTransferPolicy,
          class WeightsVectorType,
          class Communication = Dune::Amg::SequentialInformation>
class GpuOwningTwoLevelPreconditioner : public Dune::PreconditionerWithUpdate<GpuVectorType, GpuVectorType>
{
public:
    using MatrixType = typename OperatorType::matrix_type;
    using PrecFactory = Opm::PreconditionerFactory<OperatorType, Communication>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, GpuVectorType, GpuVectorType>;
    using ScalarType = typename GpuVectorType::field_type;
    using CpuVectorType = typename OperatorType::domain_type;

    GpuOwningTwoLevelPreconditioner(const OperatorType& linearoperator, const Opm::PropertyTree& prm,
                                 const std::function<WeightsVectorType()> weightsCalculator,
                                 std::size_t pressureIndex)
        : linear_operator_(linearoperator)
        , weightsCalculator_(weightsCalculator)
        , cpu_weights_(weightsCalculator_())
        , finesmoother_(PrecFactory::template createGpu<GpuVectorType>(linearoperator,
                                           prm.get_child_optional("finesmoother") ?
                                           prm.get_child("finesmoother") : Opm::PropertyTree(),
                                           std::function<CpuVectorType()>(), pressureIndex))
        , levelTransferPolicy_(dummy_comm_, cpu_weights_, prm, pressureIndex)
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
            Dune::writeMatrixMarket(cpu_weights_, outfile);
        }
    }

    virtual void pre(GpuVectorType& x, GpuVectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }

    virtual void apply(GpuVectorType& v, const GpuVectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }

    virtual void post(GpuVectorType& x) override
    {
        twolevel_method_.post(x);
    }

    virtual void update() override
    {
        cpu_weights_ = weightsCalculator_();
        updateImpl();
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return linear_operator_.category();
    }

    virtual bool hasPerfectUpdate() const override {
        return twolevel_method_.hasPerfectUpdate();
    }

    // Static methods to indicate that this preconditioner needs pre() and post() calls
    // TODO: Potentially set those to false once we have a proper implementation of the preconditioner
    static constexpr bool shouldCallPre() { return true; }
    static constexpr bool shouldCallPost() { return true; }

private:
    using CoarseOperator = typename LevelTransferPolicy::CoarseOperator;

    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperator,
                                                               FlexibleSolver<CoarseOperator>,
                                                               LevelTransferPolicy>;

    using TwoLevelMethod
        = Dune::Amg::GpuTwoLevelMethodCpr<OperatorType, CoarseSolverPolicy, Dune::PreconditionerWithUpdate<GpuVectorType, GpuVectorType>>;

    void updateImpl()
    {
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    const OperatorType& linear_operator_;
    std::function<WeightsVectorType()> weightsCalculator_;
    WeightsVectorType cpu_weights_;
    std::shared_ptr<PreconditionerWithUpdate<GpuVectorType, GpuVectorType>> finesmoother_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
    Opm::PropertyTree prm_;
    Communication dummy_comm_;

};

} // namespace Dune




#endif // OPM_OWNINGTWOLEVELPRECONDITIONERGPU_HEADER_INCLUDED
