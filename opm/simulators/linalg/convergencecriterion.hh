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
 * \copydoc Opm::ConvergenceCriterion
 */
#ifndef EWOMS_ISTL_CONVERGENCE_CRITERION_HH
#define EWOMS_ISTL_CONVERGENCE_CRITERION_HH

#include <opm/material/common/Unused.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>

#include <cmath>
#include <iostream>
#include <iomanip>

namespace Opm {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \file
 * \brief Define some base class for the convergence criteria of the linear
 * solvers of DUNE-ISTL.
 */

/*!
 * \brief Base class for all convergence criteria which only defines an virtual
 * API.
 */
template <class Vector>
class ConvergenceCriterion
{
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<typename Vector::field_type>::real_type real_type;

    typedef real_type Scalar;

public:
    /*!
     * \brief Destructor.
     *
     * In the ConvergenceCriterion it does not do anything, but it is
     * required to be declared virtual.
     */
    virtual ~ConvergenceCriterion()
    {}

    /*!
     * \brief Set the initial solution of the linear system of equations.
     *
     * This version of the method does NOT take the two-norm of the
     * residual as argument. If the two-norm of the defect is available
     * for the linear solver, the version of the update() method with it
     * should be called.
     *
     * \param curSol The current iterative solution of the linear system
     *               of equations
     * \param curResid The residual vector of the current iterative
     *                 solution of the linear system of equations
     */
    virtual void setInitial(const Vector& curSol, const Vector& curResid) = 0;

    /*!
     * \brief Update the internal members of the convergence criterion
     *        with the current solution.
     *
     * This version of the method does NOT take the two-norm of the
     * residual as argument. If the two-norm of the defect is available
     * for the linear solver, the version of the update() method with it
     * should be called.
     *
     * \param curSol The current iterative solution of the linear system
     *               of equations
     * \param changeIndicator A vector where all non-zero values indicate that the
     *                        solution has changed since the last iteration.
     * \param curResid The residual vector of the current iterative
     *                 solution of the linear system of equations
     */
    virtual void update(const Vector& curSol, const Vector& changeIndicator, const Vector& curResid) = 0;

    /*!
     * \brief Returns true if and only if the convergence criterion is
     *        met.
     */
    virtual bool converged() const = 0;

    /*!
     * \brief Returns true if the convergence criterion cannot be met anymore because the
     *        solver has broken down.
     */
    virtual bool failed() const
    { return false; }

    /*!
     * \brief Returns the accuracy of the solution at the last update.
     *
     * A value of zero means that the solution was exact.
     */
    virtual Scalar accuracy() const = 0;

    /*!
     * \brief Prints the initial information about the convergence behaviour.
     *
     * This method is called after setInitial() if the solver thinks
     * it's a good idea to be verbose. In practice, "printing the
     * initial information" means printing column headers and the
     * initial state.
     *
     * \param os The output stream to which the message gets written.
     */
    virtual void printInitial(std::ostream& os OPM_UNUSED= std::cout) const
    {}

    /*!
     * \brief Prints the information about the convergence behaviour for
     *        the current iteration.
     *
     * \param iter The iteration number. The semantics of this parameter
     *             are chosen by the linear solver.
     * \param os The output stream to which the message gets written.
     */
    virtual void print(Scalar iter OPM_UNUSED, std::ostream& os OPM_UNUSED = std::cout) const
    {}
};

//! \} end documentation

}} // end namespace Linear, Opm

#endif
