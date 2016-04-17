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
 *
 * \brief A traits class which provides basic mathematical functions for arbitrary scalar
 *        floating point values.
 *
 * The reason why this is done in such a complicated way is to enable other approaches,
 * in particular automatic differentiation based ones.
 */
#ifndef OPM_MATERIAL_MATH_TOOLBOX_HPP
#define OPM_MATERIAL_MATH_TOOLBOX_HPP

#include <cmath>
#include <algorithm>
#include <type_traits>

namespace Opm {
/*
 * \brief A traits class which provides basic mathematical functions for arbitrary scalar
 *        floating point values.
 *
 * The reason why this is done in such a complicated way is to enable other approaches,
 * in particular automatic differentiation based ones.
 */
template <class ScalarT>
struct MathToolbox
{
    static_assert(std::is_floating_point<ScalarT>::value,
                  "This class expects floating point scalars! (specialization missing?)");
public:
    /*!
     * \brief The type used to represent scalar values
     */
    typedef ScalarT Scalar;

    /*!
     * \brief The type used to represent function evaluations
     *
     * In general, that is the value of the function plus a number of derivatives at the
     * evaluation point. In the case of the scalar toolbox, no derivatives will be
     * evaluated.
     */
    typedef ScalarT Evaluation;

    /*!
     * \brief Return the value of the function at a given evaluation point.
     *
     * For this toolbox, an evaluation is the value, so this method is the identity
     * function.
     */
    static Scalar value(const Evaluation& eval)
    { return eval; }

    /*!
     * \brief Given a scalar value, return an evaluation of a constant function.
     *
     * For this toolbox, an evaluation is the value, so this method is the identity
     * function. In general, this returns an evaluation object for which all derivatives
     * are zero.
     */
    static Scalar createConstant(Scalar value)
    { return value; }

    /*!
     * \brief Given a scalar value, return an evaluation of a linear function.
     *
     * i.e., Create an evaluation which represents f(x) = x and the derivatives with
     * regard to x. For scalar evaluations (which do not consider derivatives), this
     * method does nothing.
     */
    static Scalar createVariable(Scalar value, int /* varIdx */)
    { return value; }

    /*!
     * \brief Given a function evaluation, constrain it to its value (if necessary).
     *
     * If the left hand side is a scalar and the right hand side is an evaluation, the
     * scalar gets the value of the right hand side assigned. Also if both sides are
     * scalars, this method returns the identity. The final case (left hand side being an
     * evaluation, right hand side is a scalar) yields a compiler error.
     *
     * The purpose of this method is to be able to transparantly use evaluation objects
     * in scalar computations.
     */
    template <class LhsEval>
    static LhsEval toLhs(const Evaluation& eval)
    { return eval; }

    /*!
     * \brief Pass a value through if it is an evaluation, or create a constant
     *        evaluation if it is a scalar.
     *
     * In some sense, this method is the opposite of "toLhs()": If the function argument
     * is a Scalar, an Evaluation which represents a constant value is returned, if the
     * argument is an Evaluation, it is returned as is. This method makes it possible to
     * uniformly handle the cases where some condition is either given by a constant
     * value or as a value which depends on other variables. (E.g. boundary conditions.)
     */
    static Scalar passThroughOrCreateConstant(Scalar value)
    { return value; }

    /*!
     * \brief Returns true if two values are identical up to a specified tolerance
     */
    static bool isSame(Scalar a, Scalar b, Scalar tolerance)
    {
        Scalar valueDiff = a - b;
        Scalar denom = std::max<Scalar>(1.0, std::abs(a + b));

        return std::abs(valueDiff) < tolerance || std::abs(valueDiff)/denom < tolerance;
    }

    ////////////
    // arithmetic functions
    ////////////

    //! The maximum of two arguments
    static Scalar max(Scalar arg1, Scalar arg2)
    { return std::max(arg1, arg2); }

    //! The minimum of two arguments
    static Scalar min(Scalar arg1, Scalar arg2)
    { return std::min(arg1, arg2); }

    //! The absolute value
    static Scalar abs(Scalar arg)
    { return std::abs(arg); }

    //! The tangens of a value
    static Scalar tan(Scalar arg)
    { return std::tan(arg); }

    //! The arcus tangens of a value
    static Scalar atan(Scalar arg)
    { return std::atan(arg); }

    //! The arcus tangens of a value
    static Scalar atan2(Scalar arg1, Scalar arg2)
    { return std::atan2(arg1, arg2); }

    //! The sine of a value
    static Scalar sin(Scalar arg)
    { return std::sin(arg); }

    //! The arcus sine of a value
    static Scalar asin(Scalar arg)
    { return std::asin(arg); }

    //! The cosine of a value
    static Scalar cos(Scalar arg)
    { return std::cos(arg); }

    //! The arcus cosine of a value
    static Scalar acos(Scalar arg)
    { return std::acos(arg); }

    //! The square root of a value
    static Scalar sqrt(Scalar arg)
    { return std::sqrt(arg); }

    //! The natural exponentiation of a value
    static Scalar exp(Scalar arg)
    { return std::exp(arg); }

    //! The natural logarithm of a value
    static Scalar log(Scalar arg)
    { return std::log(arg); }

    //! Exponentiation to an arbitrary base
    static Scalar pow(Scalar base, Scalar exp)
    { return std::pow(base, exp); }
};

} // namespace Opm

#endif
