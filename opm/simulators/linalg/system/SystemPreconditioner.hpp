#pragma once

#include "MultiComm.hpp"
#include "SystemTypes.hpp"

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm
{

class SystemPreconditioner : public Dune::Preconditioner<SystemVector, SystemVector>
{
public:
    using ResOperator = Dune::MatrixAdapter<RRMatrix, ResVector, ResVector>;
    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOperator>;
    using WellOperator = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemPreconditioner(const SystemMatrix& S,
                         const std::function<ResVector()>& weightCalculator,
                         int pressureIndex,
                         const Opm::PropertyTree& prm);

    void apply(SystemVector& v, const SystemVector& d) override;
    void pre(SystemVector&, SystemVector&) override
    {
    }
    void post(SystemVector&) override
    {
    }
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    void update();

private:
    const SystemMatrix& S_;
    const double well_tol_;
    const double res_tol_;

    std::unique_ptr<ResOperator> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<ResFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    WellVector wSol_;
    ResVector resSol_;
    ResVector dresSol_;
    WellVector dwSol_;
    ResVector tmp_resRes_;
    WellVector tmp_wRes_;
    ResVector resRes_;
    WellVector wRes_;
};

class SystemPreconditionerParallel : public Dune::Preconditioner<SystemVector, SystemVector>
{
public:
    using WellComm = Dune::JacComm;
    using ResComm = Dune::OwnerOverlapCopyCommunication<int, int>;
    using SystemComm = Dune::MultiCommunicator<const ResComm&, const WellComm&>;

    using ResOperator = Dune::OverlappingSchwarzOperator<RRMatrix, ResVector, ResVector, ResComm>;
    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOperator>;
    using WellOperator = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemPreconditionerParallel(const SystemMatrix& S,
                                 const std::function<ResVector()>& weightCalculator,
                                 int pressureIndex,
                                 const Opm::PropertyTree& prm,
                                 const SystemComm& syscomm);

    void apply(SystemVector& v, const SystemVector& d) override;
    void pre(SystemVector&, SystemVector&) override
    {
    }
    void post(SystemVector&) override
    {
    }
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

    void update();

private:
    const SystemMatrix& S_;
    const SystemComm& syscomm_;
    const double well_tol_;
    const double res_tol_;

    std::unique_ptr<ResOperator> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<ResFlexibleSolverType> resSolver_;
    std::unique_ptr<ResFlexibleSolverType> resSmoother_;
    std::unique_ptr<WellFlexibleSolverType> wellSolver_;

    WellVector wSol_;
    ResVector resSol_;
    ResVector dresSol_;
    WellVector dwSol_;
    ResVector tmp_resRes_;
    WellVector tmp_wRes_;
    ResVector resRes_;
    WellVector wRes_;
};

} // namespace Opm
