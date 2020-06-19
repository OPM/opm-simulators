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
 * \copydoc Opm::FixPointCriterion
 */
#ifndef EWOMS_ISTL_FIXPOINT_CRITERION_HH
#define EWOMS_ISTL_FIXPOINT_CRITERION_HH

#include "convergencecriterion.hh"

#include <opm/material/common/Unused.hpp>

namespace Opm {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \brief Provides a convergence criterion for the linear solvers
 *        which looks at the weighted maximum of the difference
 *        between two iterations.
 *
 * For the FixPointCriterion, the error of the solution is defined
 * as
 * \f[ e^k = \max_i\{ \left| w_i \delta^k_i \right| \}\;, \f]
 *
 * where \f$\delta = x^k - x^{k + 1} \f$ is the difference between
 * two consequtive iterative solution vectors \f$x^k\f$ and \f$x^{k + 1}\f$
 * and \f$w_i\f$ is the weight of the \f$i\f$-th degree of freedom.
 *
 * This criterion requires that the block type of the
 * vector is a Dune::FieldVector
 */
template <class Vector, class CollectiveCommunication>
class FixPointCriterion : public ConvergenceCriterion<Vector>
{
    using Scalar = typename Vector::field_type;
    using BlockType = typename Vector::block_type;

public:
    FixPointCriterion(const CollectiveCommunication& comm) : comm_(comm)
    {}

    FixPointCriterion(const CollectiveCommunication& comm,
                      const Vector& weightVec, Scalar reduction)
        : comm_(comm), weightVec_(weightVec), tolerance_(reduction)
    {}

    /*!
     * \brief Sets the relative weight of a primary variable
     *
     * For the FixPointCriterion, the error of the solution is defined
     * as
     * \f[ e^k = \max_i\{ \left| w_i \delta^k_i \right| \}\;, \f]
     *
     * where \f$\delta = x^k - x^{k + 1} \f$ is the difference between
     * two consequtive iterative solution vectors \f$x^k\f$ and \f$x^{k + 1}\f$
     * and \f$w_i\f$ is the weight of the \f$i\f$-th degree of freedom.
     *
     * This method is specific to the FixPointCriterion.
     *
     * \param weightVec A Dune::BlockVector<Dune::FieldVector<Scalar, n> >
     *                  with the relative weights of the degrees of freedom
     */
    void setWeight(const Vector& weightVec)
    { weightVec_ = weightVec; }

    /*!
     * \brief Return the relative weight of a primary variable
     *
     * For the FixPointCriterion, the error of the solution is defined
     * as
     * \f[ e^k = \max_i\{ \left| w_i \delta^k_i \right| \}\;, \f]
     *
     * where \f$\delta = x^k - x^{k + 1} \f$ is the difference between
     * two consequtive iterative solution vectors \f$x^k\f$ and \f$x^{k + 1}\f$
     * and \f$w_i\f$ is the weight of the \f$i\f$-th degree of freedom.
     *
     * This method is specific to the FixPointCriterion.
     *
     * \param outerIdx The index of the outer vector (i.e. Dune::BlockVector)
     * \param innerIdx The index of the inner vector (i.e. Dune::FieldVector)
     */
    Scalar weight(int outerIdx, int innerIdx) const
    { return (weightVec_.size() == 0) ? 1.0 : weightVec_[outerIdx][innerIdx]; }

    /*!
     * \brief Set the maximum allowed weighted maximum difference between two
     * iterations
     */
    /*!
     * \brief Set the maximum allowed maximum difference between two
     *        iterationsfor the solution considered to be converged.
     */
    void setTolerance(Scalar tol)
    { tolerance_ = tol; }

    /*!
     * \brief Return the maximum allowed weighted difference between two
     *        iterations for the solution considered to be converged.
     */
    Scalar tolerance() const
    { return tolerance_; }

    /*!
     * \copydoc ConvergenceCriterion::setInitial(const Vector&, const Vector&)
     */
    void setInitial(const Vector& curSol, const Vector& curResid OPM_UNUSED)
    {
        lastSol_ = curSol;
        delta_ = 1000 * tolerance_;
    }

    /*!
     * \copydoc ConvergenceCriterion::update(const Vector&, const Vector&, const Vector&)
     */
    void update(const Vector& curSol,
                const Vector& changeIndicator OPM_UNUSED,
                const Vector& curResid OPM_UNUSED)
    {
        assert(curSol.size() == lastSol_.size());

        delta_ = 0.0;
        for (size_t i = 0; i < curSol.size(); ++i) {
            for (size_t j = 0; j < BlockType::dimension; ++j) {
                delta_ =
                    std::max(delta_, weight(i, j)*std::abs(curSol[i][j] - lastSol_[i][j]));
            }
        }

        delta_ = comm_.max(delta_);
        lastSol_ = curSol;
    }

    /*!
     * \copydoc ConvergenceCriterion::converged()
     */
    bool converged() const
    { return accuracy() < tolerance(); }

    /*!
     * \copydoc ConvergenceCriterion::accuracy()
     */
    Scalar accuracy() const
    { return delta_; }

private:
    const CollectiveCommunication& comm_;

    Vector lastSol_;   // solution of the last iteration
    Vector weightVec_; // solution of the last iteration
    Scalar delta_;     // the maximum of the absolute weighted difference of the
                       // last two iterations
    Scalar tolerance_; // the maximum allowed delta for the solution to be
                       // considered converged
};

//! \} end documentation

}} // end namespace Linear, Opm

#endif
