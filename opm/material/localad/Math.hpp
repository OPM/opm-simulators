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
namespace LocalAd {
// provide some algebraic functions
template <class Scalar, int numVars>
Evaluation<Scalar, numVars> abs(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::abs(x.value);

    // derivatives use the chain rule
    if (x.value < 0.0) {
        for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
            result.derivatives[curVarIdx] = -x.derivatives[curVarIdx];
    }
    else {
        for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
            result.derivatives[curVarIdx] = x.derivatives[curVarIdx];
    }

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> min(const Evaluation<Scalar, numVars>& x1,
                                const Evaluation<Scalar, numVars>& x2)
{
    Evaluation<Scalar, numVars> result;

    if (x1.value < x2.value) {
        result.value = x1.value;

        std::copy(x1.derivatives.begin(),
                  x1.derivatives.end(),
                  result.derivatives.begin());
    }
    else  {
        result.value = x2.value;

        std::copy(x2.derivatives.begin(),
                  x2.derivatives.end(),
                  result.derivatives.begin());
    }

    return result;
}

template <class ScalarA, class Scalar, int numVars>
Evaluation<Scalar, numVars> min(ScalarA x1,
                                const Evaluation<Scalar, numVars>& x2)
{
    Evaluation<Scalar, numVars> result;

    if (x1 < x2.value) {
        result.value = x1;

        std::fill(result.derivatives.begin(),
                  result.derivatives.end(),
                  0.0);
    }
    else  {
        result.value = x2.value;

        std::copy(x2.derivatives.begin(),
                  x2.derivatives.end(),
                  result.derivatives.begin());
    }

    return result;
}

template <class Arg2Scalar, class Scalar, int numVars>
Evaluation<Scalar, numVars> min(const Evaluation<Scalar, numVars>& x1,
                                const Arg2Scalar& x2)
{ return min(x2, x1); }

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> max(const Evaluation<Scalar, numVars>& x1,
                                const Evaluation<Scalar, numVars>& x2)
{
    Evaluation<Scalar, numVars> result;

    if (x1.value > x2.value) {
        result.value = x1.value;

        std::copy(x1.derivatives.begin(),
                  x1.derivatives.end(),
                  result.derivatives.begin());
    }
    else  {
        result.value = x2.value;

        std::copy(x2.derivatives.begin(),
                  x2.derivatives.end(),
                  result.derivatives.begin());
    }

    return result;
}

template <class Arg1Scalar, class Scalar, int numVars>
Evaluation<Scalar, numVars> max(const Arg1Scalar& x1,
                                const Evaluation<Scalar, numVars>& x2)
{
    Evaluation<Scalar, numVars> result;

    if (x1 > x2.value) {
        result.value = x1;

        std::fill(result.derivatives.begin(),
                  result.derivatives.end(),
                  0.0);
    }
    else  {
        result.value = x2.value;

        std::copy(x2.derivatives.begin(),
                  x2.derivatives.end(),
                  result.derivatives.begin());
    }

    return result;
}

template <class Arg2Scalar, class Scalar, int numVars>
Evaluation<Scalar, numVars> max(const Evaluation<Scalar, numVars>& x1,
                                const Arg2Scalar& x2)
{ return max(x2, x1); }

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> tan(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    const Scalar& tmp = std::tan(x.value);
    result.value = tmp;

    // derivatives use the chain rule
    Scalar df_dx = 1 + tmp*tmp;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> atan(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::atan(x.value);

    // derivatives use the chain rule
    Scalar df_dx = 1/(1 + x.value*x.value);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> atan2(const Evaluation<Scalar, numVars>& x,
                                  const Evaluation<Scalar, numVars>& y)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::atan2(x.value, y.value);

    // derivatives use the chain rule
    Scalar alpha = 1/(1 + (x.value*x.value)/(y.value*y.value));
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx) {
        result.derivatives[curVarIdx] =
            alpha
            /(y.value*y.value)
            *(x.derivatives[curVarIdx]*y.value - x.value*y.derivatives[curVarIdx]);
    }

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> sin(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::sin(x.value);

    // derivatives use the chain rule
    Scalar df_dx = std::cos(x.value);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> asin(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::asin(x.value);

    // derivatives use the chain rule
    Scalar df_dx = 1.0/std::sqrt(1 - x.value*x.value);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> cos(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::cos(x.value);

    // derivatives use the chain rule
    Scalar df_dx = -std::sin(x.value);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> acos(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::acos(x.value);

    // derivatives use the chain rule
    Scalar df_dx = - 1.0/std::sqrt(1 - x.value*x.value);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> sqrt(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    Scalar sqrt_x = std::sqrt(x.value);
    result.value = sqrt_x;

    // derivatives use the chain rule
    Scalar df_dx = 0.5/sqrt_x;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx) {
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];
    }

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> exp(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    Scalar exp_x = std::exp(x.value);
    result.value = exp_x;

    // derivatives use the chain rule
    Scalar df_dx = exp_x;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

// exponentiation of arbitrary base with a fixed constant
template <class Scalar, int numVars, class ExpType>
Evaluation<Scalar, numVars> pow(const Evaluation<Scalar, numVars>& base,
                                const ExpType& exp)
{
    Evaluation<Scalar, numVars> result;

    Scalar pow_x = std::pow(base.value, exp);
    result.value = pow_x;

    // derivatives use the chain rule
    Scalar df_dx = pow_x/base.value*exp;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*base.derivatives[curVarIdx];

    return result;
}

// exponentiation of constant base with an arbitrary exponent
template <class BaseType, class Scalar, int numVars>
Evaluation<Scalar, numVars> pow(const BaseType& base,
                                const Evaluation<Scalar, numVars>& exp)
{
    Evaluation<Scalar, numVars> result;

    Scalar lnBase = std::log(base);
    result.value = std::exp(lnBase*exp.value);

    // derivatives use the chain rule
    Scalar df_dx = lnBase*result.value;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*exp.derivatives[curVarIdx];

    return result;
}

// this is the most expensive power function. Computationally it is pretty expensive, so
// one of the above two variants above should be preferred if possible.
template <class Scalar, int numVars>
Evaluation<Scalar, numVars> pow(const Evaluation<Scalar, numVars>& base,
                                const Evaluation<Scalar, numVars>& exp)
{
    Evaluation<Scalar, numVars> result;

    Scalar valuePow = std::pow(base.value, exp.value);
    result.value = valuePow;

    // use the chain rule for the derivatives. since both, the base and the exponent can
    // potentially depend on the variable set, calculating these is quite elaborate...
    Scalar f = base.value;
    Scalar g = exp.value;
    Scalar logF = std::log(f);
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx) {
        Scalar fPrime = base.derivatives[curVarIdx];
        Scalar gPrime = exp.derivatives[curVarIdx];
        result.derivatives[curVarIdx] = (g*fPrime/f + logF*gPrime) * valuePow;
    }

    return result;
}

template <class Scalar, int numVars>
Evaluation<Scalar, numVars> log(const Evaluation<Scalar, numVars>& x)
{
    Evaluation<Scalar, numVars> result;

    result.value = std::log(x.value);

    // derivatives use the chain rule
    Scalar df_dx = 1/x.value;
    for (unsigned curVarIdx = 0; curVarIdx < result.derivatives.size(); ++curVarIdx)
        result.derivatives[curVarIdx] = df_dx*x.derivatives[curVarIdx];

    return result;
}

} // namespace LocalAd

// a kind of traits class for the automatic differentiation case. (The toolbox for the
// scalar case is provided by the MathToolbox.hpp header file.)
template <class ScalarT, int numVars>
struct MathToolbox<Opm::LocalAd::Evaluation<ScalarT, numVars> >
{
private:
public:
    typedef ScalarT Scalar;
    typedef Opm::LocalAd::Evaluation<ScalarT, numVars> Evaluation;

    static Scalar value(const Evaluation& eval)
    { return eval.value; }

    static Evaluation createConstant(Scalar value)
    { return Evaluation::createConstant(value); }

    static Evaluation createVariable(Scalar value, int varIdx)
    { return Evaluation::createVariable(value, varIdx); }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<Evaluation, LhsEval>::value,
                                   LhsEval>::type
    toLhs(const Evaluation& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<Evaluation, LhsEval>::value,
                                   LhsEval>::type
    toLhs(const Evaluation&& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_floating_point<LhsEval>::value,
                                   LhsEval>::type
    toLhs(const Evaluation& eval)
    { return eval.value; }

    static const Evaluation passThroughOrCreateConstant(Scalar value)
    { return createConstant(value); }

    static const Evaluation& passThroughOrCreateConstant(const Evaluation& eval)
    { return eval; }


    // arithmetic functions
    template <class Arg1Eval, class Arg2Eval>
    static Evaluation max(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return Opm::LocalAd::max(arg1, arg2); }

    template <class Arg1Eval, class Arg2Eval>
    static Evaluation min(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return Opm::LocalAd::min(arg1, arg2); }

    static Evaluation abs(const Evaluation& arg)
    { return Opm::LocalAd::abs(arg); }

    static Evaluation tan(const Evaluation& arg)
    { return Opm::LocalAd::tan(arg); }

    static Evaluation atan(const Evaluation& arg)
    { return Opm::LocalAd::atan(arg); }

    static Evaluation atan2(const Evaluation& arg1, const Evaluation& arg2)
    { return Opm::LocalAd::atan2(arg1, arg2); }

    static Evaluation sin(const Evaluation& arg)
    { return Opm::LocalAd::sin(arg); }

    static Evaluation asin(const Evaluation& arg)
    { return Opm::LocalAd::asin(arg); }

    static Evaluation cos(const Evaluation& arg)
    { return Opm::LocalAd::cos(arg); }

    static Evaluation acos(const Evaluation& arg)
    { return Opm::LocalAd::acos(arg); }

    static Evaluation sqrt(const Evaluation& arg)
    { return Opm::LocalAd::sqrt(arg); }

    static Evaluation exp(const Evaluation& arg)
    { return Opm::LocalAd::exp(arg); }

    static Evaluation log(const Evaluation& arg)
    { return Opm::LocalAd::log(arg); }

    static Evaluation pow(const Evaluation& arg1, typename Evaluation::Scalar arg2)
    { return Opm::LocalAd::pow(arg1, arg2); }

    static Evaluation pow(typename Evaluation::Scalar arg1, const Evaluation& arg2)
    { return Opm::LocalAd::pow(arg1, arg2); }

    static Evaluation pow(const Evaluation& arg1, const Evaluation& arg2)
    { return Opm::LocalAd::pow(arg1, arg2); }
};

}

#endif
