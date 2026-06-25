/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_PRECONDITIONER2INVERSEOPERATOR_HEADER_INCLUDED
#define OPM_PRECONDITIONER2INVERSEOPERATOR_HEADER_INCLUDED

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

#endif // OPM_PRECONDITIONER2INVERSEOPERATOR_HEADER_INCLUDED