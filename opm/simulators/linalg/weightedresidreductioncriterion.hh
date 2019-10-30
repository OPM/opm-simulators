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
 * \copydoc Opm::Linear::WeightedResidualReductionCriterion
 */
#ifndef EWOMS_WEIGHTED_RESIDUAL_REDUCTION_CRITERION_HH
#define EWOMS_WEIGHTED_RESIDUAL_REDUCTION_CRITERION_HH

#include "convergencecriterion.hh"

#include <iostream>

namespace Opm {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \brief Convergence criterion which looks at the weighted absolute
 *        value of the residual
 *
 * For the WeightedResidualReductionCriterion, the error of the
 * solution is defined as
 * \f[ e^k = \max_i\{ \left| w_i r^k_i \right| \}\;, \f]
 *
 * where \f$r^k = \mathbf{A} x^k - b \f$ is the residual for the
 * k-th iterative solution vector \f$x^k\f$ and \f$w_i\f$ is the
 * weight of the \f$i\f$-th linear equation.
 *
 * In addition, to the reduction of the maximum defect, the linear
 * solver is also considered to be converged, if the defect goes below
 * a given absolute limit.
 */
template <class Vector, class CollectiveCommunication>
class WeightedResidualReductionCriterion : public ConvergenceCriterion<Vector>
{
    typedef typename Vector::field_type Scalar;
    typedef typename Vector::block_type BlockType;

public:
    WeightedResidualReductionCriterion(const CollectiveCommunication& comm)
        : comm_(comm)
    {}

    WeightedResidualReductionCriterion(const CollectiveCommunication& comm,
                                       const Vector& residWeights,
                                       Scalar residualReductionTolerance,
                                       Scalar fixPointTolerance,
                                       Scalar absResidualTolerance = 0.0,
                                       Scalar maxError = 0.0)
        : comm_(comm),
          residWeightVec_(residWeights),
          fixPointTolerance_(fixPointTolerance),
          residualReductionTolerance_(residualReductionTolerance),
          absResidualTolerance_(absResidualTolerance),
          maxError_(maxError)
    { }

    /*!
     * \brief Sets the relative weight of each row of the residual.
     *
     * For the WeightedResidualReductionCriterion, the error of the solution is
     * defined as
     * \f[ e^k = \max_i\{ \left| w_i r^k_i \right| \}\;, \f]
     *
     * where \f$r^k = \mathbf{A} x^k - b \f$ is the residual for the
     * k-th iterative solution vector \f$x^k\f$ and \f$w_i\f$ is the
     * weight of the \f$i\f$-th linear equation.
     *
     * This method is not part of the generic ConvergenceCriteria interface.
     *
     * \param residWeightVec A Dune::BlockVector<Dune::FieldVector<Scalar, n> >
     *                  with the relative weights of the linear equations
     */
    void setResidualWeight(const Vector& residWeightVec)
    { residWeightVec_ = residWeightVec; }

    /*!
     * \brief Return the relative weight of a row of the residual.
     *
     * \param outerIdx The index of the outer vector (i.e. Dune::BlockVector)
     * \param innerIdx The index of the inner vector (i.e. Dune::FieldVector)
     */
    Scalar residualWeight(size_t outerIdx, unsigned innerIdx) const
    {
        return (residWeightVec_.size() == 0)
                   ? 1.0
                   : residWeightVec_[outerIdx][innerIdx];
    }

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
     * \brief Sets the maximum absolute tolerated residual.
     */
    void setResidualTolerance(Scalar tol)
    { absResidualTolerance_ = tol; }

    /*!
     * \brief Returns the maximum absolute tolerated residual.
     */
    Scalar absResidualTolerance() const
    { return absResidualTolerance_; }

    /*!
     * \brief Returns the reduction of the weighted maximum of the
     *        residual compared to the initial solution.
     */
    Scalar residualAccuracy() const
    { return residualError_/std::max<Scalar>(1e-20, initialResidualError_); }

    /*!
     * \brief Sets the fix-point tolerance.
     */
    void setFixPointTolerance(Scalar tol)
    { fixPointTolerance_ = tol; }

    /*!
     * \brief Returns the maximum tolerated difference between two
     *        iterations to be met before a solution is considered to be
     *        converged.
     */
    Scalar fixPointTolerance() const
    { return fixPointTolerance_; }

    /*!
     * \brief Returns the weighted maximum of the difference
     *        between the last two iterative solutions.
     */
    Scalar fixPointAccuracy() const
    { return fixPointError_; }

    /*!
     * \copydoc ConvergenceCriterion::setInitial(const Vector& , const Vector& )
     */
    void setInitial(const Vector& curSol, const Vector& curResid)
    {
        lastResidualError_ = 1e100;

        lastSolVec_ = curSol;
        updateErrors_(curSol, curResid);
        // the fix-point error is not applicable for the initial solution!
        fixPointError_ = 1e100;

        // make sure that we don't allow an initial error of 0 to avoid
        // divisions by zero
        residualError_ = std::max<Scalar>(residualError_, 1e-20);
        initialResidualError_ = residualError_;
    }

    /*!
     * \copydoc ConvergenceCriterion::update(const Vector& , const Vector& )
     */
    void update(const Vector& curSol,
                const Vector& changeIndicator OPM_UNUSED,
                const Vector& curResid)
    {
        lastResidualError_ = residualError_;
        updateErrors_(curSol, curResid);
    }

    /*!
     * \copydoc ConvergenceCriterion::converged()
     */
    bool converged() const
    {
        // we're converged if the solution is better than the tolerance
        // fix-point and residual tolerance.
        return
            residualAccuracy() <= residualReductionTolerance() ||
            residualError_ <= absResidualTolerance_;
    }

    /*!
     * \copydoc ConvergenceCriterion::failed()
     */
    bool failed() const
    {
        return
            (!converged() && fixPointError_ <= fixPointTolerance_)
            || residualError_ > maxError_;
    }

    /*!
     * \copydoc ConvergenceCriterion::accuracy()
     *
     * For the accuracy we only take the residual into account,
     */
    Scalar accuracy() const
    { return residualError_; }

    /*!
     * \copydoc ConvergenceCriterion::printInitial()
     */
    void printInitial(std::ostream& os = std::cout) const
    {
        os << std::setw(20) << " Iter ";
        os << std::setw(20) << " Delta ";
        os << std::setw(20) << " Residual ";
        os << std::setw(20) << " ResidRed ";
        os << std::setw(20) << " Rate ";
        os << std::endl;
    }

    /*!
     * \copydoc ConvergenceCriterion::print()
     */
    void print(Scalar iter, std::ostream& os = std::cout) const
    {
        static constexpr Scalar eps = std::numeric_limits<Scalar>::min()*1e10;

        os << std::setw(20) << iter << " ";
        os << std::setw(20) << fixPointAccuracy() << " ";
        os << std::setw(20) << residualError_ << " ";
        os << std::setw(20) << 1.0/residualAccuracy() << " ";
        os << std::setw(20) << lastResidualError_ / std::max<Scalar>(residualError_, eps) << " ";
        os << std::endl << std::flush;
    }

private:
    // update the weighted absolute residual
    void updateErrors_(const Vector& curSol, const Vector& curResid)
    {
        residualError_ = 0.0;
        fixPointError_ = 0.0;
        for (size_t i = 0; i < curResid.size(); ++i) {
            for (unsigned j = 0; j < BlockType::dimension; ++j) {
                residualError_ =
                    std::max<Scalar>(residualError_,
                                     residualWeight(i, j)*std::abs(curResid[i][j]));
                fixPointError_ =
                    std::max<Scalar>(fixPointError_,
                                     std::abs(curSol[i][j] - lastSolVec_[i][j])
                                     /std::max<Scalar>(1.0, curSol[i][j]));
            }
        }
        lastSolVec_ = curSol;

        residualError_ = comm_.max(residualError_);
        fixPointError_ = comm_.max(fixPointError_);
    }

    const CollectiveCommunication& comm_;

    // the weights of the components of the residual
    Vector residWeightVec_;

    // solution vector of the last iteration
    Vector lastSolVec_;

    // the maximum of the weighted difference between the last two
    // iterations
    Scalar fixPointError_;

    // the maximum allowed relative tolerance for difference of the
    // solution of two iterations
    Scalar fixPointTolerance_;

    // the maximum of the absolute weighted residual of the last
    // iteration
    Scalar residualError_;

    // the maximum of the absolute weighted difference of the last
    // iteration
    Scalar lastResidualError_;

    // the maximum of the absolute weighted residual of the initial
    // solution
    Scalar initialResidualError_;

    // the maximum allowed relative tolerance of the residual for the
    // solution to be considered converged
    Scalar residualReductionTolerance_;

    // the maximum allowed absolute tolerance of the residual for the
    // solution to be considered converged
    Scalar absResidualTolerance_;

    // The maximum error which is tolerated before we fail.
    Scalar maxError_;
};

//! \} end documentation

}} // end namespace Linear, Opm

#endif
