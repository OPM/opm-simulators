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

#include <opm/material/common/Unused.hpp>

#include <cmath>
#include <algorithm>
#include <type_traits>
#include <stdexcept>

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
     * \brief The type used to represent "primitive" scalar values
     */
    typedef ScalarT Scalar;

    /*!
     * \brief The type used to represent values
     *
     * In general, these objects represent the function value at a given point plus a
     * number of derivatives. In the case of the scalars, no derivatives will be
     * evaluated.
     */
    typedef ScalarT ValueType;

    /*!
     * \brief The toolbox for the type of value objects.
     *
     * For this class this is trivial because primitive floating point objects are
     * endpoints in the "nesting graph". This typedef makes sense if nested automatic
     * differentiation is used, though...
     */
    typedef MathToolbox<Scalar> InnerToolbox;

    /*!
     * \brief Return the value of the function at a given evaluation point.
     *
     * For this toolbox, there are no derivatives so this method is the identity
     * function.
     */
    static Scalar value(Scalar value)
    { return value; }

    /*!
     * \brief Return the primitive scalar value of a value object.
     *
     * Since this toolbox's value objects are primitive scalars, this method just passes
     * through the argument it was given.
     */
    static Scalar scalarValue(Scalar value)
    { return value; }

    /*!
     * \brief Given a scalar value, return a "compatible" object
     *
     * This basically boils down to creating an uninitialized object of sufficient size.
     * This is method only non-trivial for dynamically-sized Evaluation objects.
     */
    static Scalar createBlank(Scalar value OPM_UNUSED)
    { return Scalar(); }

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
     * \brief Given a scalar value, return an evaluation of a constant function that
     *        features a given number of derivatives
     *
     * For this toolbox, an evaluation is the value, so this method is the identity
     * function. In general, this returns an evaluation object for which all derivatives
     * are zero.
     */
    static Scalar createConstant(unsigned numDerivatives, Scalar value)
    {
        if (numDerivatives != 0)
            throw std::logic_error("Plain floating point objects cannot represent any derivatives");
        return value;
    }

    /*!
     * \brief Given a scalar value, return an evaluation of a constant function that is
     *        compatible to a "template" variable.
     *
     * For this toolbox, an evaluation is the value, so this method is the identity
     * function. In general, this returns an evaluation object for which all derivatives
     * are zero.
     */
    static Scalar createConstant(Scalar x OPM_UNUSED, Scalar value)
    { return value; }

    /*!
     * \brief Given a scalar value, return an evaluation of a linear function.
     *
     * i.e., Create an evaluation which represents f(x) = x and the derivatives with
     * regard to x. For scalars (which do not consider derivatives), this method does
     * nothing.
     */
    static Scalar createVariable(Scalar value OPM_UNUSED, unsigned varIdx OPM_UNUSED)
    { throw std::logic_error("Plain floating point objects cannot represent variables"); }

    /*!
     * \brief Given a scalar value, return an evaluation of a linear function that is
     *        compatible with a "template" evaluation.
     *
     * i.e., Create an evaluation which represents f(x) = x and the derivatives with
     * regard to x. For scalars (which do not consider derivatives), this method does
     * nothing.
     */
    static Scalar createVariable(Scalar x OPM_UNUSED, Scalar value OPM_UNUSED, unsigned varIdx OPM_UNUSED)
    { throw std::logic_error("Plain floating point objects cannot represent variables"); }

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
    static LhsEval decay(Scalar value)
    {
        static_assert(std::is_floating_point<LhsEval>::value,
                      "The left-hand side must be a primitive floating point type!");

        return value;
    }

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

    //! The 10 logarithm of a value
    static Scalar log10(Scalar arg)
    { return std::log10(arg); }

    //! The natural logarithm of a value
    static Scalar log(Scalar arg)
    { return std::log(arg); }

    //! Exponentiation to an arbitrary base
    static Scalar pow(Scalar base, Scalar exp)
    { return std::pow(base, exp); }

    //! Return true iff the argument's value and all its derivatives are finite values
    static bool isfinite(Scalar arg)
    { return std::isfinite(arg); }

    //! Return true iff the argument's value or any of its derivatives are NaN values
    static bool isnan(Scalar arg)
    { return std::isnan(arg); }
};

template <class Eval1, class Eval2>
struct ReturnEval_
{
    typedef typename std::remove_const< typename std::remove_reference<Eval1>::type >::type T;
    typedef typename std::remove_const< typename std::remove_reference<Eval2>::type >::type U;

    //static_assert(std::is_constructible<T, U>::value || std::is_constructible<U, T>::value,
    //              "One of the argument types must be constructible to the other");

    typedef typename std::conditional<std::is_constructible<T, U>::value,
                                      T,
                                      U>::type type;
};

// these are convenience functions for not having to type MathToolbox<Scalar>::foo()
template <class Evaluation>
Evaluation blank(const Evaluation& x)
{ return MathToolbox<Evaluation>::createBlank(x); }

template <class Evaluation, class Scalar>
Evaluation constant(const Scalar& value)
{ return MathToolbox<Evaluation>::createConstant(value); }

template <class Evaluation, class Scalar>
Evaluation constant(unsigned numDeriv, const Scalar& value)
{ return MathToolbox<Evaluation>::createConstant(numDeriv, value); }

template <class Evaluation, class Scalar>
Evaluation constant(const Evaluation& x, const Scalar& value)
{ return MathToolbox<Evaluation>::createConstant(x, value); }

template <class Evaluation, class Scalar>
Evaluation variable(unsigned numDeriv, const Scalar& value, unsigned idx)
{ return MathToolbox<Evaluation>::createVariable(numDeriv, value, idx); }

template <class Evaluation, class Scalar>
Evaluation variable(const Evaluation& x, const Scalar& value, unsigned idx)
{ return MathToolbox<Evaluation>::createVariable(x, value, idx); }

template <class Evaluation, class Scalar>
Evaluation variable(const Scalar& value, unsigned idx)
{ return MathToolbox<Evaluation>::createVariable(value, idx); }

template <class ResultEval, class Evaluation>
auto decay(const Evaluation& value)
    -> decltype(MathToolbox<Evaluation>::template decay<ResultEval>(value))
{ return MathToolbox<Evaluation>::template decay<ResultEval>(value); }

template <class Evaluation>
auto getValue(const Evaluation& val)
    -> decltype(MathToolbox<Evaluation>::value(val))
{ return MathToolbox<Evaluation>::value(val); }

template <class Evaluation>
auto scalarValue(const Evaluation& val)
    -> decltype(MathToolbox<Evaluation>::scalarValue(val))
{ return MathToolbox<Evaluation>::scalarValue(val); }

template <class Evaluation1, class Evaluation2>
typename ReturnEval_<Evaluation1, Evaluation2>::type
max(const Evaluation1& arg1, const Evaluation2& arg2)
{ return MathToolbox<typename ReturnEval_<Evaluation1, Evaluation2>::type>::max(arg1, arg2); }

template <class Evaluation1, class Evaluation2>
typename ReturnEval_<Evaluation1, Evaluation2>::type
min(const Evaluation1& arg1, const Evaluation2& arg2)
{ return MathToolbox<typename ReturnEval_<Evaluation1, Evaluation2>::type>::min(arg1, arg2); }

template <class Evaluation>
Evaluation abs(const Evaluation& value)
{ return MathToolbox<Evaluation>::abs(value); }

template <class Evaluation>
Evaluation tan(const Evaluation& value)
{ return MathToolbox<Evaluation>::tan(value); }

template <class Evaluation>
Evaluation atan(const Evaluation& value)
{ return MathToolbox<Evaluation>::atan(value); }

template <class Evaluation1, class Evaluation2>
typename ReturnEval_<Evaluation1, Evaluation2>::type
atan2(const Evaluation1& value1, const Evaluation2& value2)
{ return MathToolbox<typename ReturnEval_<Evaluation1, Evaluation2>::type>::atan2(value1, value2); }

template <class Evaluation>
Evaluation sin(const Evaluation& value)
{ return MathToolbox<Evaluation>::sin(value); }

template <class Evaluation>
Evaluation asin(const Evaluation& value)
{ return MathToolbox<Evaluation>::asin(value); }

template <class Evaluation>
Evaluation cos(const Evaluation& value)
{ return MathToolbox<Evaluation>::cos(value); }

template <class Evaluation>
Evaluation acos(const Evaluation& value)
{ return MathToolbox<Evaluation>::acos(value); }

template <class Evaluation>
Evaluation sqrt(const Evaluation& value)
{ return MathToolbox<Evaluation>::sqrt(value); }

template <class Evaluation>
Evaluation exp(const Evaluation& value)
{ return MathToolbox<Evaluation>::exp(value); }

template <class Evaluation>
Evaluation log(const Evaluation& value)
{ return MathToolbox<Evaluation>::log(value); }

template <class Evaluation>
Evaluation log10(const Evaluation& value)
{ return MathToolbox<Evaluation>::log10(value); }

template <class Evaluation1, class Evaluation2>
typename ReturnEval_<Evaluation1, Evaluation2>::type
pow(const Evaluation1& base, const Evaluation2& exp)
{ return MathToolbox<typename ReturnEval_<Evaluation1, Evaluation2>::type>::pow(base, exp); }

template <class Evaluation>
bool isfinite(const Evaluation& value)
{ return MathToolbox<Evaluation>::isfinite(value); }

template <class Evaluation>
bool isnan(const Evaluation& value)
{ return MathToolbox<Evaluation>::isnan(value); }

} // namespace Opm

#endif

