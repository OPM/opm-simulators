#pragma once

#include <dune/common/shared_ptr.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>

#include <memory>

namespace Dune
{

/*!
    \brief Adapter exposing a preconditioner through the inverse-operator interface.

    Each apply() call performs one preconditioner application. Unlike
    LoopSolver it does not update the defect with the linear operator or run
    iterative convergence checks.
 */
template<class X>
class Preconditioner2InverseOperator : public InverseOperator<X, X>
{
public:
    using typename InverseOperator<X, X>::domain_type;
    using typename InverseOperator<X, X>::range_type;
    using typename InverseOperator<X, X>::field_type;
    using typename InverseOperator<X, X>::real_type;
    using typename InverseOperator<X, X>::scalar_real_type;

    using PreconditionerType = Preconditioner<X, X>;

    explicit Preconditioner2InverseOperator(PreconditionerType& prec)
        : prec_(stackobject_to_shared_ptr(prec))
    {
    }

    explicit Preconditioner2InverseOperator(std::shared_ptr<PreconditionerType> prec)
        : prec_(std::move(prec))
    {
    }

    void apply(X& x, X& b, InverseOperatorResult& res) override
    {
        prec_->pre(x, b);
        // Preconditioners compute a correction/update, so start from zero to
        // interpret this adapter as x = M^{-1} b.
        x = 0;
        prec_->apply(x, b);
        prec_->post(x);

        // No residual norm is evaluated here; report one successful
        // application rather than an iterative solve.
        res.clear();
        res.iterations = 1;
        res.reduction = 1.0;
        res.converged = true;
    }

    void apply(X& x, X& b, [[maybe_unused]] double reduction, InverseOperatorResult& res) override
    {
        apply(x, b, res);
    }

    SolverCategory::Category category() const override
    {
        return SolverCategory::category(*prec_);
    }

private:
    // Non-owning alias when constructed from a reference, shared ownership
    // when constructed from a shared_ptr.
    std::shared_ptr<PreconditionerType> prec_;
};

} // namespace Dune