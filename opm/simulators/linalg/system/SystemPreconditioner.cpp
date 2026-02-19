#include <config.h>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>

namespace Opm
{

SystemPreconditioner::SystemPreconditioner(const SystemMatrix& S,
                                           const std::function<ResVector()>& weightsCalculator,
                                           int pressureIndex,
                                           const Opm::PropertyTree& prm)
    : S_(S)
    , well_tol_(prm.get<double>("well_solver.tol"))
    , res_tol_(prm.get<double>("reservoir_solver.tol"))
{
    auto rop = std::make_unique<ResOperator>(S[_0][_0]);
    auto wop = std::make_unique<WellOperator>(S[_1][_1]);
    auto resprm = prm.get_child("reservoir_solver");
    auto resprmsmoother = prm.get_child("reservoir_smoother");
    auto wellprm = prm.get_child("well_solver");

    auto rsol
        = std::make_unique<ResFlexibleSolverType>(*rop, resprm, weightsCalculator, pressureIndex);
    auto rsmoother = std::make_unique<ResFlexibleSolverType>(
        *rop, resprmsmoother, weightsCalculator, pressureIndex);

    std::function<WellVector()> weightsCalculatorWell;
    auto wsol = std::make_unique<WellFlexibleSolverType>(
        *wop, wellprm, weightsCalculatorWell, pressureIndex);

    rop_ = std::move(rop);
    wop_ = std::move(wop);
    resSolver_ = std::move(rsol);
    resSmoother_ = std::move(rsmoother);
    wellSolver_ = std::move(wsol);

    const auto numRes = S[_0][_0].N();
    const auto numWell = S[_1][_1].N();
    wSol_.resize(numWell);
    resSol_.resize(numRes);
    dresSol_.resize(numRes);
    dwSol_.resize(numWell);
    tmp_resRes_.resize(numRes);
    tmp_wRes_.resize(numWell);
    resRes_.resize(numRes);
    wRes_.resize(numWell);
}

void
SystemPreconditioner::update()
{
    resSolver_->preconditioner().update();
    resSmoother_->preconditioner().update();
    wellSolver_->preconditioner().update();
}

void
SystemPreconditioner::apply(SystemVector& v, const SystemVector& d)
{
    const auto& A = S_[_0][_0];
    const auto& C = S_[_1][_0];
    const auto& B = S_[_0][_1];
    const auto& D = S_[_1][_1];

    resRes_ = d[_0];
    wRes_ = d[_1];
    resSol_ = 0.0;
    wSol_ = 0.0;

    // Stage 1: Reservoir CPR solve
    {
        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        Dune::InverseOperatorResult res_result;
        resSolver_->apply(dresSol_, tmp_resRes_, res_tol_, res_result);
        resSol_ += dresSol_;
        A.mmv(dresSol_, resRes_);
        C.mmv(dresSol_, wRes_);
    }

    // Stage 2: Well solve + reservoir system smoothing
    {
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        Dune::InverseOperatorResult well_result;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_tol_, well_result);
        wSol_ += dwSol_;
        B.mmv(dwSol_, resRes_);
        D.mmv(dwSol_, wRes_);

        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        Dune::InverseOperatorResult res_result;
        resSmoother_->apply(dresSol_, tmp_resRes_, res_tol_, res_result);
        resSol_ += dresSol_;
        // Only update wRes_ (needed for final well solve); resRes_ is not used again
        C.mmv(dresSol_, wRes_);
    }

    // Stage 3: Final well solve
    {
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        Dune::InverseOperatorResult well_result;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_tol_, well_result);
        wSol_ += dwSol_;
    }

    v[_0] = resSol_;
    v[_1] = wSol_;
}


SystemPreconditionerParallel::SystemPreconditionerParallel(
    const SystemMatrix& S,
    const std::function<ResVector()>& weightsCalculator,
    int pressureIndex,
    const Opm::PropertyTree& prm,
    const SystemComm& syscomm)
    : S_(S)
    , syscomm_(syscomm)
    , well_tol_(prm.get<double>("well_solver.tol"))
    , res_tol_(prm.get<double>("reservoir_solver.tol"))
{
    auto rop = std::make_unique<ResOperator>(S[_0][_0], syscomm_[_0]);
    auto wop = std::make_unique<WellOperator>(S[_1][_1]);
    auto resprm = prm.get_child("reservoir_solver");
    auto wellprm = prm.get_child("well_solver");
    auto resprmsmoother = prm.get_child("reservoir_smoother");

    auto rsol = std::make_unique<ResFlexibleSolverType>(
        *rop, syscomm_[_0], resprm, weightsCalculator, pressureIndex);
    auto rsmoother = std::make_unique<ResFlexibleSolverType>(
        *rop, syscomm_[_0], resprmsmoother, weightsCalculator, pressureIndex);

    std::function<WellVector()> weightsCalculatorWell;
    auto wsol = std::make_unique<WellFlexibleSolverType>(
        *wop, wellprm, weightsCalculatorWell, pressureIndex);

    rop_ = std::move(rop);
    wop_ = std::move(wop);
    resSolver_ = std::move(rsol);
    resSmoother_ = std::move(rsmoother);
    wellSolver_ = std::move(wsol);

    const auto numRes = S[_0][_0].N();
    const auto numWell = S[_1][_1].N();
    wSol_.resize(numWell);
    resSol_.resize(numRes);
    dresSol_.resize(numRes);
    dwSol_.resize(numWell);
    tmp_resRes_.resize(numRes);
    tmp_wRes_.resize(numWell);
    resRes_.resize(numRes);
    wRes_.resize(numWell);
}

void
SystemPreconditionerParallel::update()
{
    resSolver_->preconditioner().update();
    resSmoother_->preconditioner().update();
    wellSolver_->preconditioner().update();
}

void
SystemPreconditionerParallel::apply(SystemVector& v, const SystemVector& d)
{
    const auto& A = S_[_0][_0];
    const auto& C = S_[_1][_0];
    const auto& B = S_[_0][_1];
    const auto& D = S_[_1][_1];

    resRes_ = d[_0];
    wRes_ = d[_1];
    resSol_ = 0.0;
    wSol_ = 0.0;

    // Stage 1: Reservoir CPR solve
    {
        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        syscomm_[_0].copyOwnerToAll(tmp_resRes_, tmp_resRes_);
        Dune::InverseOperatorResult res_result;
        resSolver_->apply(dresSol_, tmp_resRes_, res_tol_, res_result);
        resSol_ += dresSol_;
        A.mmv(dresSol_, resRes_);
        C.mmv(dresSol_, wRes_);
    }

    // Stage 2: Well solve + reservoir system smoothing
    {
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        Dune::InverseOperatorResult well_result;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_tol_, well_result);
        wSol_ += dwSol_;
        B.mmv(dwSol_, resRes_);
        D.mmv(dwSol_, wRes_);

        dresSol_ = 0.0;
        tmp_resRes_ = resRes_;
        syscomm_[_0].copyOwnerToAll(tmp_resRes_, tmp_resRes_);
        Dune::InverseOperatorResult res_result;
        resSmoother_->apply(dresSol_, tmp_resRes_, res_tol_, res_result);
        resSol_ += dresSol_;
        // Only update wRes_ (needed for final well solve); resRes_ is not used again
        C.mmv(dresSol_, wRes_);
    }

    // Stage 3: Final well solve
    {
        dwSol_ = 0.0;
        tmp_wRes_ = wRes_;
        Dune::InverseOperatorResult well_result;
        wellSolver_->apply(dwSol_, tmp_wRes_, well_tol_, well_result);
        wSol_ += dwSol_;
    }

    syscomm_[_0].copyOwnerToAll(resSol_, resSol_);
    v[_0] = resSol_;
    v[_1] = wSol_;
}

} // namespace Opm

template class Dune::FlexibleSolver<
    Dune::MatrixAdapter<Opm::WWMatrix, Opm::WellVector, Opm::WellVector>>;
