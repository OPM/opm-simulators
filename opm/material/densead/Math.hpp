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
 * \brief A number of commonly used algebraic functions for the localized OPM automatic
 *        differentiation (AD) framework.
 *
 * This file provides AD variants of the the most commonly used functions of the <cmath>
 * header file.
 */
#ifndef OPM_LOCAL_AD_MATH_HPP
#define OPM_LOCAL_AD_MATH_HPP

#include "Evaluation.hpp"

#include <opm/material/common/MathToolbox.hpp>

namespace Opm {
namespace DenseAd {
// forward declaration of the Evaluation template class
template <class ValueT, int numVars, unsigned staticSize>
class Evaluation;

// provide some algebraic functions
template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> abs(const Evaluation<ValueType, numVars, staticSize>& x)
{ return (x > 0.0)?x:-x; }

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> min(const Evaluation<ValueType, numVars, staticSize>& x1,
                                               const Evaluation<ValueType, numVars, staticSize>& x2)
{ return (x1 < x2)?x1:x2; }

template <class Arg1ValueType, class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> min(const Arg1ValueType& x1,
                                               const Evaluation<ValueType, numVars, staticSize>& x2)
{
    if (x1 < x2) {
        Evaluation<ValueType, numVars, staticSize> ret(x2);
        ret = x1;
        return ret;
    }
    else
        return x2;
}

template <class ValueType, int numVars, unsigned staticSize, class Arg2ValueType>
Evaluation<ValueType, numVars, staticSize> min(const Evaluation<ValueType, numVars, staticSize>& x1,
                                               const Arg2ValueType& x2)
{ return min(x2, x1); }

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> max(const Evaluation<ValueType, numVars, staticSize>& x1,
                                               const Evaluation<ValueType, numVars, staticSize>& x2)
{ return (x1 > x2)?x1:x2; }

template <class Arg1ValueType, class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> max(const Arg1ValueType& x1,
                                               const Evaluation<ValueType, numVars, staticSize>& x2)
{
    if (x1 > x2) {
        Evaluation<ValueType, numVars, staticSize> ret(x2);
        ret = x1;
        return ret;
    }
    else
        return x2;
}

template <class ValueType, int numVars, unsigned staticSize, class Arg2ValueType>
Evaluation<ValueType, numVars, staticSize> max(const Evaluation<ValueType, numVars, staticSize>& x1,
                                               const Arg2ValueType& x2)
{ return max(x2, x1); }

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> tan(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    const ValueType& tmp = ValueTypeToolbox::tan(x.value());
    result.setValue(tmp);

    // derivatives use the chain rule
    const ValueType& df_dx = 1 + tmp*tmp;
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> atan(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::atan(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1/(1 + x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> atan2(const Evaluation<ValueType, numVars, staticSize>& x,
                                                 const Evaluation<ValueType, numVars, staticSize>& y)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::atan2(x.value(), y.value()));

    // derivatives use the chain rule
    const ValueType& alpha = 1/(1 + (x.value()*x.value())/(y.value()*y.value()));
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx) {
        result.setDerivative(curVarIdx,
                             alpha/(y.value()*y.value())
                             *(x.derivative(curVarIdx)*y.value() - x.value()*y.derivative(curVarIdx)));
    }

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> atan2(const Evaluation<ValueType, numVars, staticSize>& x,
                                                 const ValueType& y)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::atan2(x.value(), y));

    // derivatives use the chain rule
    const ValueType& alpha = 1/(1 + (x.value()*x.value())/(y*y));
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx) {
        result.setDerivative(curVarIdx,
                             alpha/(y*y)
                             *(x.derivative(curVarIdx)*y));
    }

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> atan2(const ValueType& x,
                                                 const Evaluation<ValueType, numVars, staticSize>& y)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(y);

    result.setValue(ValueTypeToolbox::atan2(x, y.value()));

    // derivatives use the chain rule
    const ValueType& alpha = 1/(1 + (x.value()*x.value())/(y.value()*y.value()));
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx) {
        result.setDerivative(curVarIdx,
                             alpha/(y.value()*y.value())
                             *x*y.derivative(curVarIdx));
    }

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> sin(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::sin(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = ValueTypeToolbox::cos(x.value());
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> asin(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::asin(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1.0/ValueTypeToolbox::sqrt(1 - x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> cos(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::cos(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = -ValueTypeToolbox::sin(x.value());
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> acos(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::acos(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = - 1.0/ValueTypeToolbox::sqrt(1 - x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> sqrt(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    const ValueType& sqrt_x = ValueTypeToolbox::sqrt(x.value());
    result.setValue(sqrt_x);

    // derivatives use the chain rule
    ValueType df_dx = 0.5/sqrt_x;
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx) {
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));
    }

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> exp(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;
    Evaluation<ValueType, numVars, staticSize> result(x);

    const ValueType& exp_x = ValueTypeToolbox::exp(x.value());
    result.setValue(exp_x);

    // derivatives use the chain rule
    const ValueType& df_dx = exp_x;
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

// exponentiation of arbitrary base with a fixed constant
template <class ValueType, int numVars, unsigned staticSize, class ExpType>
Evaluation<ValueType, numVars, staticSize> pow(const Evaluation<ValueType, numVars, staticSize>& base,
                                               const ExpType& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;
    Evaluation<ValueType, numVars, staticSize> result(base);

    const ValueType& pow_x = ValueTypeToolbox::pow(base.value(), exp);
    result.setValue(pow_x);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result = 0.0;
    }
    else {
        // derivatives use the chain rule
        const ValueType& df_dx = pow_x/base.value()*exp;
        for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
            result.setDerivative(curVarIdx, df_dx*base.derivative(curVarIdx));
    }

    return result;
}

// exponentiation of constant base with an arbitrary exponent
template <class BaseType, class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> pow(const BaseType& base,
                                               const Evaluation<ValueType, numVars, staticSize>& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(exp);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result = 0.0;
    }
    else {
        const ValueType& lnBase = ValueTypeToolbox::log(base);
        result.setValue(ValueTypeToolbox::exp(lnBase*exp.value()));

        // derivatives use the chain rule
        const ValueType& df_dx = lnBase*result.value();
        for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
            result.setDerivative(curVarIdx, df_dx*exp.derivative(curVarIdx));
    }

    return result;
}

// this is the most expensive power function. Computationally it is pretty expensive, so
// one of the above two variants above should be preferred if possible.
template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> pow(const Evaluation<ValueType, numVars, staticSize>& base,
                                               const Evaluation<ValueType, numVars, staticSize>& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(base);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result = 0.0;
    }
    else {
        ValueType valuePow = ValueTypeToolbox::pow(base.value(), exp.value());
        result.setValue(valuePow);

        // use the chain rule for the derivatives. since both, the base and the exponent can
        // potentially depend on the variable set, calculating these is quite elaborate...
        const ValueType& f = base.value();
        const ValueType& g = exp.value();
        const ValueType& logF = ValueTypeToolbox::log(f);
        for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx) {
            const ValueType& fPrime = base.derivative(curVarIdx);
            const ValueType& gPrime = exp.derivative(curVarIdx);
            result.setDerivative(curVarIdx, (g*fPrime/f + logF*gPrime) * valuePow);
        }
    }

    return result;
}

template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> log(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::log(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1/x.value();
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}


template <class ValueType, int numVars, unsigned staticSize>
Evaluation<ValueType, numVars, staticSize> log10(const Evaluation<ValueType, numVars, staticSize>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    Evaluation<ValueType, numVars, staticSize> result(x);

    result.setValue(ValueTypeToolbox::log10(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1/x.value() * ValueTypeToolbox::log10(ValueTypeToolbox::exp(1.0));
    for (int curVarIdx = 0; curVarIdx < result.size(); ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

} // namespace DenseAd

// a kind of traits class for the automatic differentiation case. (The toolbox for the
// scalar case is provided by the MathToolbox.hpp header file.)
template <class ValueT, int numVars, unsigned staticSize>
struct MathToolbox<DenseAd::Evaluation<ValueT, numVars, staticSize> >
{
private:
public:
    typedef ValueT ValueType;
    typedef MathToolbox<ValueType> InnerToolbox;
    typedef typename InnerToolbox::Scalar Scalar;
    typedef DenseAd::Evaluation<ValueType, numVars, staticSize> Evaluation;

    static ValueType value(const Evaluation& eval)
    { return eval.value(); }

    static decltype(InnerToolbox::scalarValue(0.0)) scalarValue(const Evaluation& eval)
    { return InnerToolbox::scalarValue(eval.value()); }

    static Evaluation createBlank(const Evaluation& x)
    { return Evaluation::createBlank(x); }

    static Evaluation createConstantZero(const Evaluation& x)
    { return Evaluation::createConstantZero(x); }

    static Evaluation createConstantOne(const Evaluation& x)
    { return Evaluation::createConstantOne(x); }

    static Evaluation createConstant(ValueType value)
    { return Evaluation::createConstant(value); }

    static Evaluation createConstant(unsigned numDeriv, const ValueType value)
    { return Evaluation::createConstant(numDeriv, value); }

    static Evaluation createConstant(const Evaluation& x, const ValueType value)
    { return Evaluation::createConstant(x, value); }

    static Evaluation createVariable(ValueType value, int varIdx)
    { return Evaluation::createVariable(value, varIdx); }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<Evaluation, LhsEval>::value,
                                   LhsEval>::type
    decay(const Evaluation& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<Evaluation, LhsEval>::value,
                                   LhsEval>::type
    decay(const Evaluation&& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_floating_point<LhsEval>::value,
                                   LhsEval>::type
    decay(const Evaluation& eval)
    { return eval.value(); }

    // comparison
    static bool isSame(const Evaluation& a, const Evaluation& b, Scalar tolerance)
    {
        typedef MathToolbox<ValueType> ValueTypeToolbox;

        // make sure that the value of the evaluation is identical
        if (!ValueTypeToolbox::isSame(a.value(), b.value(), tolerance))
            return false;

        // make sure that the derivatives are identical
        for (int curVarIdx = 0; curVarIdx < numVars; ++curVarIdx)
            if (!ValueTypeToolbox::isSame(a.derivative(curVarIdx), b.derivative(curVarIdx), tolerance))
                return false;

        return true;
    }

    // arithmetic functions
    template <class Arg1Eval, class Arg2Eval>
    static Evaluation max(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return DenseAd::max(arg1, arg2); }

    template <class Arg1Eval, class Arg2Eval>
    static Evaluation min(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return DenseAd::min(arg1, arg2); }

    static Evaluation abs(const Evaluation& arg)
    { return DenseAd::abs(arg); }

    static Evaluation tan(const Evaluation& arg)
    { return DenseAd::tan(arg); }

    static Evaluation atan(const Evaluation& arg)
    { return DenseAd::atan(arg); }

    static Evaluation atan2(const Evaluation& arg1, const Evaluation& arg2)
    { return DenseAd::atan2(arg1, arg2); }

    template <class Eval2>
    static Evaluation atan2(const Evaluation& arg1, const Eval2& arg2)
    { return DenseAd::atan2(arg1, arg2); }

    template <class Eval1>
    static Evaluation atan2(const Eval1& arg1, const Evaluation& arg2)
    { return DenseAd::atan2(arg1, arg2); }

    static Evaluation sin(const Evaluation& arg)
    { return DenseAd::sin(arg); }

    static Evaluation asin(const Evaluation& arg)
    { return DenseAd::asin(arg); }

    static Evaluation cos(const Evaluation& arg)
    { return DenseAd::cos(arg); }

    static Evaluation acos(const Evaluation& arg)
    { return DenseAd::acos(arg); }

    static Evaluation sqrt(const Evaluation& arg)
    { return DenseAd::sqrt(arg); }

    static Evaluation exp(const Evaluation& arg)
    { return DenseAd::exp(arg); }

    static Evaluation log(const Evaluation& arg)
    { return DenseAd::log(arg); }

    static Evaluation log10(const Evaluation& arg)
    { return DenseAd::log10(arg); }

    template <class RhsValueType>
    static Evaluation pow(const Evaluation& arg1, const RhsValueType& arg2)
    { return DenseAd::pow(arg1, arg2); }

    template <class RhsValueType>
    static Evaluation pow(const RhsValueType& arg1, const Evaluation& arg2)
    { return DenseAd::pow(arg1, arg2); }

    static Evaluation pow(const Evaluation& arg1, const Evaluation& arg2)
    { return DenseAd::pow(arg1, arg2); }

    static bool isfinite(const Evaluation& arg)
    {
        if (!InnerToolbox::isfinite(arg.value()))
            return false;

        for (int i = 0; i < numVars; ++i)
            if (!InnerToolbox::isfinite(arg.derivative(i)))
                return false;

        return true;
    }

    static bool isnan(const Evaluation& arg)
    {
        if (InnerToolbox::isnan(arg.value()))
            return true;

        for (int i = 0; i < numVars; ++i)
            if (InnerToolbox::isnan(arg.derivative(i)))
                return true;

        return false;
    }
};

}

#endif
