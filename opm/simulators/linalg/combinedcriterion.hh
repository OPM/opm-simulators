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
 * \copydoc Opm::Linear::CombinedCriterion
 */
#ifndef EWOMS_COMBINED_CRITERION_HH
#define EWOMS_COMBINED_CRITERION_HH

#include "convergencecriterion.hh"

#include <iostream>

namespace Opm {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \brief Convergence criterion which looks at the absolute value of the residual and
 *        fails if the linear solver stagnates.
 *
 * For the CombinedCriterion, the error of the solution is defined as \f[ e^k = \max_i\{
 * \left| r^k_i \right| \}\;, \f]
 *
 * where \f$r^k = \mathbf{A} x^k - b \f$ is the residual for the k-th iterative solution
 * vector \f$x^k\f$.
 *
 * In addition, to the reduction of the maximum residual, the linear solver is aborted
 * early if the residual goes below or above absolute limits.
 */
template <class Vector, class CollectiveCommunication>
class CombinedCriterion : public ConvergenceCriterion<Vector>
{
    typedef typename Vector::field_type Scalar;
    typedef typename Vector::block_type BlockType;

public:
    CombinedCriterion(const CollectiveCommunication& comm)
        : comm_(comm)
    {}

    CombinedCriterion(const CollectiveCommunication& comm,
                                       Scalar residualReductionTolerance,
                                       Scalar absResidualTolerance = 0.0,
                                       Scalar maxResidual = 0.0)
        : comm_(comm),
          residualReductionTolerance_(residualReductionTolerance),
          absResidualTolerance_(absResidualTolerance),
          maxResidual_(maxResidual)
    { }

    /*!
     * \brief Sets the residual reduction tolerance.
     */
    void setResidualReductionTolerance(Scalar tol)
    { residualReductionTolerance_ = tol; }

    /*!
     * \brief Returns the tolerance of the residual reduction of the solution.
     */
    Scalar residualReductionTolerance() const
    { return residualReductionTolerance_; }

    /*!
     * \brief Returns the reduction of the maximum of the residual compared to the
     *        initial solution.
     */
    Scalar residualReduction() const
    { return residualError_/std::max<Scalar>(1e-20, initialResidualError_); }

    /*!
     * \brief Sets the maximum absolute tolerated residual.
     */
    void setAbsResidualTolerance(Scalar tol)
    { absResidualTolerance_ = tol; }

    /*!
     * \brief Returns the tolerated maximum of the the infinity norm of the absolute
     *        residual.
     */
    Scalar absResidualTolerance() const
    { return absResidualTolerance_; }

    /*!
     * \brief Returns the infinity norm of the absolute residual.
     */
    Scalar absResidual() const
    { return residualError_; }

    /*!
     * \copydoc ConvergenceCriterion::setInitial(const Vector& , const Vector& )
     */
    void setInitial(const Vector& curSol, const Vector& curResid) override
    {
        updateErrors_(curSol, curSol, curResid);
        stagnates_ = false;

        // to avoid divisions by zero, make sure that we don't use an initial error of 0
        residualError_ = std::max<Scalar>(residualError_,
                                          std::numeric_limits<Scalar>::min()*1e10);
        initialResidualError_ = residualError_;
        lastResidualError_ = residualError_;
    }

    /*!
     * \copydoc ConvergenceCriterion::update(const Vector&, const Vector&, const Vector&)
     */
    void update(const Vector& curSol, const Vector& changeIndicator, const Vector& curResid) override
    { updateErrors_(curSol, changeIndicator, curResid);  }

    /*!
     * \copydoc ConvergenceCriterion::converged()
     */
    bool converged() const override
    {
        // we're converged if the solution is better than the tolerance
        // fix-point and residual tolerance.
        return
            residualReduction() <= residualReductionTolerance() ||
            absResidual() <= absResidualTolerance();
    }

    /*!
     * \copydoc ConvergenceCriterion::failed()
     */
    bool failed() const override
    { return !converged() && (stagnates_ || residualError_ > maxResidual_); }

    /*!
     * \copydoc ConvergenceCriterion::accuracy()
     *
     * For the accuracy we only take the residual into account,
     */
    Scalar accuracy() const override
    { return residualError_/initialResidualError_; }

    /*!
     * \copydoc ConvergenceCriterion::printInitial()
     */
    void printInitial(std::ostream& os = std::cout) const override
    {
        os << std::setw(20) << "iteration ";
        os << std::setw(20) << "residual ";
        os << std::setw(20) << "reduction ";
        os << std::setw(20) << "rate ";
        os << std::endl;
    }

    /*!
     * \copydoc ConvergenceCriterion::print()
     */
    void print(Scalar iter, std::ostream& os = std::cout) const override
    {
        const Scalar eps = std::numeric_limits<Scalar>::min()*1e10;

        os << std::setw(20) << iter << " ";
        os << std::setw(20) << absResidual() << " ";
        os << std::setw(20) << accuracy() << " ";
        os << std::setw(20) << lastResidualError_/std::max<Scalar>(residualError_, eps) << " ";
        os << std::endl << std::flush;
    }

private:
    // update the weighted absolute residual
    void updateErrors_(const Vector& curSol OPM_UNUSED, const Vector& changeIndicator,  const Vector& curResid)
    {
        lastResidualError_ = residualError_;
        residualError_ = 0.0;
        stagnates_ = true;
        for (size_t i = 0; i < curResid.size(); ++i) {
            for (unsigned j = 0; j < BlockType::dimension; ++j) {
                residualError_ =
                    std::max<Scalar>(residualError_,
                                     std::abs(curResid[i][j]));

                if (stagnates_ && changeIndicator[i][j] != 0.0)
                    // only stagnation means that we've failed!
                    stagnates_ = false;
            }
        }

        residualError_ = comm_.max(residualError_);

        // the linear solver only stagnates if all processes stagnate
        stagnates_ = comm_.min(stagnates_);
    }

    const CollectiveCommunication& comm_;

    // the infinity norm of the residual of the last iteration
    Scalar lastResidualError_;

    // the infinity norm of the residual of the current iteration
    Scalar residualError_;

    // the infinity norm of the residual of the initial solution
    Scalar initialResidualError_;

    // the minimum reduction of the residual norm where the solution is to be considered
    // converged
    Scalar residualReductionTolerance_;

    // the maximum residual norm for the residual for the solution to be considered to be
    // converged
    Scalar absResidualTolerance_;

    // The maximum error which is tolerated before we fail.
    Scalar maxResidual_;

    // does the linear solver seem to stagnate, i.e. were the last two solutions
    // identical?
    bool stagnates_;
};

//! \} end documentation

}} // end namespace Linear, Opm

#endif
