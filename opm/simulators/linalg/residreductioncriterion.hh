// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::ResidReductionCriterion
 */
#ifndef EWOMS_RESID_REDUCTION_CRITERION_HH
#define EWOMS_RESID_REDUCTION_CRITERION_HH

#include "convergencecriterion.hh"

#include <opm/material/common/Unused.hpp>

#include <dune/istl/scalarproducts.hh>

namespace Opm {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \brief Provides a convergence criterion which looks at the
 *        reduction of the two-norm of the residual for the linear
 *        solvers.
 *
 * For the ResidReductionCriterion, the error of the solution is defined
 * as
 * \f[ e^k = \frac{\left| A x_k - b \right|}{\left| A x_0 - b \right|}\;, \f]
 */
template <class Vector>
class ResidReductionCriterion : public ConvergenceCriterion<Vector>
{
    using Scalar = typename Vector::field_type;

public:
    ResidReductionCriterion(Dune::ScalarProduct<Vector>& scalarProduct,
                            Scalar tolerance = 1e-6)
        : scalarProduct_(scalarProduct), tolerance_(tolerance)
    {}

    /*!
     * \brief Set the maximum allowed weighted maximum of the reduction of the
     * linear residual.
     */
    void setTolerance(Scalar tol)
    { tolerance_ = tol; }

    /*!
     * \brief Return the maximum allowed weighted maximum of the reduction of the linear residual.
     */
    Scalar tolerance() const
    { return tolerance_; }

    /*!
     * \copydoc ConvergenceCriterion::setInitial(const Vector& , const Vector& )
     */
    void setInitial(const Vector& curSol OPM_UNUSED, const Vector& curResid)
    {
        static constexpr Scalar eps = std::numeric_limits<Scalar>::min()*1e10;

        // make sure that we don't allow an initial error of 0 to avoid
        // divisions by zero
        curDefect_ = scalarProduct_.norm(curResid);
        lastDefect_ = curDefect_;
        initialDefect_ = std::max(curDefect_, eps);
    }

    /*!
     * \copydoc ConvergenceCriterion::update(const Vector& , const Vector& )
     */
    void update(const Vector& curSol OPM_UNUSED,
                const Vector& changeIndicator OPM_UNUSED,
                const Vector& curResid)
    {
        lastDefect_ = curDefect_;
        curDefect_ = scalarProduct_.norm(curResid);
    }

    /*!
     * \copydoc ConvergenceCriterion::converged()
     */
    bool converged() const
    { return accuracy() <= tolerance(); }

    /*!
     * \copydoc ConvergenceCriterion::accuracy()
     */
    Scalar accuracy() const
    { return curDefect_/initialDefect_; }

    /*!
     * \copydoc ConvergenceCriterion::printInitial()
     */
    void printInitial(std::ostream& os=std::cout) const
    {
        os << std::setw(20) << "iteration ";
        os << std::setw(20) << "residual ";
        os << std::setw(20) << "accuracy ";
        os << std::setw(20) << "rate ";
        os << std::endl;
    }

    /*!
     * \copydoc ConvergenceCriterion::print()
     */
    void print(Scalar iter, std::ostream& os=std::cout) const
    {
        static constexpr Scalar eps = std::numeric_limits<Scalar>::min()*1e10;

        os << std::setw(20) << iter << " ";
        os << std::setw(20) << curDefect_ << " ";
        os << std::setw(20) << accuracy() << " ";
        os << std::setw(20) << (lastDefect_/std::max(eps, curDefect_)) << " ";
        os << std::endl;
    }

private:
    Dune::ScalarProduct<Vector>& scalarProduct_;

    Scalar tolerance_;
    Scalar initialDefect_;
    Scalar curDefect_;
    Scalar lastDefect_;
};

//! \} end documentation

}} // end namespace Linear, Opm

#endif
