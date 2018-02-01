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
 * \brief Low-level tests for the localized automatic differentiation (AD) framework.
 */
#include "config.h"

// for testing the "!=" and "==" operators, we need to disable the -Wfloat-equal to
// prevent clang from producing a warning with -Weverything
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/material/common/Unused.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <stdexcept>

//static const int numVars = 3;

template <class Scalar, int numVars>
struct TestEnv
{

void testOperators(const Scalar tolerance)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;

    // test the constructors of the Opm::DenseAd::Evaluation class
    const Scalar c = 1.234;
    const Scalar x = 4.567;
    const Scalar y = 8.910;
    const Eval cEval = Eval::createConstant(c);
    const Eval c2Eval OPM_UNUSED = c;
    const Eval xEval = Eval::createVariable(x, 0);
    const Eval yEval = Eval::createVariable(y, 1);

    Eval xyEval = xEval;
    Eval yxEval = yEval;

    xyEval.copyDerivatives(yxEval);
    yxEval.clearDerivatives();

    for (int i = 0; i < numVars; ++i) {
        if (i == 0 && xEval.derivative(i) != 1.0)
            throw std::logic_error("oops: createVariable");
        else if (i == 1 && yEval.derivative(i) != 1.0)
            throw std::logic_error("oops: createVariable");
        else if (i == 1 && xyEval.derivative(i) != 1.0)
            throw std::logic_error("oops: copyDerivatives");
        else continue;

        // the remaining derivatives must be zero
        if (xEval.derivative(i) != 0.0)
            throw std::logic_error("oops: createVariable");
        if (yEval.derivative(i) != 0.0)
            throw std::logic_error("oops: createVariable");
        if (cEval.derivative(i) != 0.0)
            throw std::logic_error("oops: createConstant");
        if (xyEval.derivative(i) != 0.0)
            throw std::logic_error("oops: copyDerivatives");
        if (xyEval.derivative(i) != 0.0)
            throw std::logic_error("oops: clearDerivatives");
    }

    // test the non-inplace operators
    {
        Eval a = xEval + yEval;
        if (std::abs(a.value() - (x + y)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval b = xEval + c;
        if (std::abs(b.value() - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval d = xEval + cEval;
        if (std::abs(d.value() - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");
    }

    {
        Eval a = xEval - yEval;
        if (std::abs(a.value() - (x - y)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval b = xEval - c;
        if (std::abs(b.value() - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval d = xEval - cEval;
        if (std::abs(d.value() - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");
    }

    {
        Eval a = xEval*yEval;
        if (std::abs(a.value() - (x*y)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval b = xEval*c;
        if (std::abs(b.value() - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval d = xEval*cEval;
        if (std::abs(d.value() - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");
    }

    {
        Eval a = xEval/yEval;
        if (std::abs(a.value() - (x/y)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval b = xEval/c;
        if (std::abs(b.value() - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval d = xEval/cEval;
        if (std::abs(d.value() - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");
    }

    // test the inplace operators
    {
        Eval a = xEval;
        a += yEval;
        if (std::abs(a.value() - (x + y)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval b = xEval;
        b += c;
        if (std::abs(b.value() - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval d = xEval;
        d += cEval;
        if (std::abs(d.value() - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");
    }

    {
        Eval a = xEval;
        a -= yEval;
        if (std::abs(a.value() - (x - y)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval b = xEval;
        b -= c;
        if (std::abs(b.value() - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval d = xEval;
        d -= cEval;
        if (std::abs(d.value() - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");
    }

    {
        Eval a = xEval;
        a *= yEval;
        if (std::abs(a.value() - (x*y)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval b = xEval;
        b *= c;
        if (std::abs(b.value() - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval d = xEval;
        d *= cEval;
        if (std::abs(d.value() - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");
    }

    {
        Eval a = xEval;
        a /= yEval;
        if (std::abs(a.value() - (x/y)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval b = xEval;
        b /= c;
        if (std::abs(b.value() - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval d = xEval;
        d /= cEval;
        if (std::abs(d.value() - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");
    }

    {
        Eval a = 1.0;
        Eval b = 2.0;
        if (a >= b)
            throw std::logic_error("oops: operator>=");
        if (a >= 2.0)
            throw std::logic_error("oops: operator>=");

        if (!(a >= a))
            throw std::logic_error("oops: operator>=");
        if (!(a >= 1.0))
            throw std::logic_error("oops: operator>=");
        if (!(1.0 <= a))
            throw std::logic_error("oops: operator<=");

        if (b <= a)
            throw std::logic_error("oops: operator<=");
        if (b <= 1.0)
            throw std::logic_error("oops: operator<=");

        if (!(b <= b))
            throw std::logic_error("oops: operator<=");
        if (!(b <= 2.0))
            throw std::logic_error("oops: operator<=");
        if (!(2.0 >= b))
            throw std::logic_error("oops: operator>=");

        if (a != a)
            throw std::logic_error("oops: operator!=");
        if (a != 1.0)
            throw std::logic_error("oops: operator!=");
        if (1.0 != a)
            throw std::logic_error("oops: operator!=");
    }
}

template <class AdFn, class ClassicFn>
void test1DFunction(AdFn* adFn, ClassicFn* classicFn, Scalar xMin = 1e-6, Scalar xMax = 1000)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar x = Scalar(i)/(n - 1)*(xMax - xMin) + xMin;

        const auto& xEval = Eval::createVariable(x, 0);
        const Eval& yEval = adFn(xEval);

        Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
        eps = std::max(eps, eps*std::abs(x));

        Scalar y = classicFn(x);
        Scalar yStar1 = classicFn(x - eps);
        Scalar yStar2 = classicFn(x + eps);
        Scalar yPrime = (yStar2 - yStar1)/(2*eps);

        if (std::abs(y-yEval.value()) > 5e-14)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(yPrime - yEval.derivative(0));
        Scalar deltaRel = std::abs(deltaAbs/yPrime);
        if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) yPrime) + " vs "
                                   + std::to_string((long double) yEval.derivative(0))
                                   + " delta: " + std::to_string((long double) std::abs(yPrime - yEval.derivative(0))));
    }
}

template <class AdFn,
          class ClassicFn>
void test2DFunction1(AdFn* adFn, ClassicFn* classicFn, Scalar xMin, Scalar xMax, Scalar y)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar x = Scalar(i)/(n - 1)*(xMax - xMin) + xMin;

        const auto& xEval = Eval::createVariable(x, 0);
        const auto& yEval = Eval::createConstant(y);
        const Eval& zEval = adFn(xEval, yEval);

        Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
        eps = std::max(eps, eps*std::abs(x));

        Scalar z = classicFn(x, y);
        Scalar zStar1 = classicFn(x - eps, y);
        Scalar zStar2 = classicFn(x + eps, y);
        Scalar zPrime = (zStar2 - zStar1)/(2.*eps);

        if (std::abs(z - zEval.value())/std::abs(z + zEval.value())
            > std::numeric_limits<Scalar>::epsilon()*1e2)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval.derivative(0));
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval.derivative(0))
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivative(0))));
    }
}

template <class AdFn,
          class ClassicFn>
void test2DFunction2(AdFn* adFn, ClassicFn* classicFn, Scalar x, Scalar yMin, Scalar yMax)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar y = Scalar(i)/(n - 1)*(yMax - yMin) + yMin;

        const auto& xEval = Eval::createConstant(x);
        const auto& yEval = Eval::createVariable(y, 1);
        const Eval& zEval = adFn(xEval, yEval);

        Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
        eps = std::max(eps, eps*std::abs(y));

        Scalar z = classicFn(x, y);
        Scalar zStar1 = classicFn(x, y - eps);
        Scalar zStar2 = classicFn(x, y + eps);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (std::abs(z - zEval.value())/std::abs(z + zEval.value())
            > std::numeric_limits<Scalar>::epsilon()*1e2)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval.derivative(1));
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval.derivative(1))
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivative(1))));
    }
}

void testPowBase(Scalar baseMin = 1e-2, Scalar baseMax = 100)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;
    typedef Opm::MathToolbox<Eval> EvalToolbox;

    Scalar exp = 1.234;
    const auto& expEval = Eval::createConstant(exp);

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar base = Scalar(i)/(n - 1)*(baseMax - baseMin) + baseMin;

        const auto& baseEval = Eval::createVariable(base, 0);
        const Eval& zEval1 = pow(baseEval, exp);
        const Eval& zEval2 = pow(baseEval, expEval);

        Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
        eps = std::max(eps, eps*std::abs(base));

        Scalar z = pow(base, exp);
        Scalar zStar1 = pow(base - eps, exp);
        Scalar zStar2 = pow(base + eps, exp);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (std::abs(z - zEval2.value())/std::abs(z + zEval2.value())
            > std::numeric_limits<Scalar>::epsilon()*1e2)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval1.derivative(0));
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) base)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval1.derivative(0))
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval1.derivative(0))));

        if (!EvalToolbox::isSame(zEval1, zEval2, /*tolerance=*/std::numeric_limits<Scalar>::epsilon()*1e3*zEval1.value()))
            throw std::logic_error("oops: pow(Eval, Scalar) != pow(Eval, Eval)");
    }
}

void testPowExp(Scalar expMin = -100, Scalar expMax = 100)
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;
    typedef Opm::MathToolbox<Eval> EvalToolbox;

    Scalar base = 1.234;
    const auto& baseEval = Eval::createConstant(base);

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar exp = Scalar(i)/(n - 1)*(expMax - expMin) + expMin;
        const auto& expEval = Eval::createVariable(exp, 1);

        const Eval& zEval1 = pow(base, expEval);
        const Eval& zEval2 = pow(baseEval, expEval);

        Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
        eps = std::max(eps, eps*std::abs(exp));

        Scalar z = pow(base, exp);
        Scalar zStar1 = pow(base, exp - eps);
        Scalar zStar2 = pow(base, exp + eps);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (std::abs(z - zEval2.value())/std::abs(z + zEval2.value())
            > std::numeric_limits<Scalar>::epsilon()*1e2)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval1.derivative(1));
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) base)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval1.derivative(1))
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval1.derivative(1))));

        if (!EvalToolbox::isSame(zEval1, zEval2, /*tolerance=*/std::numeric_limits<Scalar>::epsilon()*1e3*zEval1.value()))
            throw std::logic_error("oops: pow(Eval, Scalar) != pow(Eval, Eval)");
    }
}

void testAtan2()
{
    typedef Opm::DenseAd::Evaluation<Scalar, numVars> Eval;

    int n = 1000;
    Scalar maxVal = 10.0;
    for (int i = 1; i < n; ++ i) {
        Scalar x = 2*maxVal*Scalar(i)/n - maxVal;
        if (- 0.05 < x && x < 0.05)
            // avoid numerical problems
            continue;

        const Eval& xEval = Eval::createVariable(x, 0);

        for (int j = 1; j < n; ++ j) {
            Scalar y = 2*maxVal*Scalar(j)/n - maxVal;

            if (- 0.05 < y && y < 0.05)
                // avoid numerical problems
                continue;

            const Eval& yEval = Eval::createVariable(y, 0);
            const Eval& zEval = atan2(xEval, yEval);

            Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7;
            eps = std::max(eps, eps*std::abs(x));

            Scalar z = atan2(x, y);
            Scalar zStar1 = atan2(x - eps, y - eps);
            Scalar zStar2 = atan2(x + eps, y + eps);
            Scalar zPrime = (zStar2 - zStar1)/(2*eps);

            if (std::abs(z - zEval.value())/std::abs(z + zEval.value())
                > std::numeric_limits<Scalar>::epsilon()*1e2)
                throw std::logic_error("oops: value");

            Scalar deltaAbs = std::abs(zPrime - zEval.derivative(0));
            Scalar deltaRel = std::abs(deltaAbs/zPrime);
            if (deltaAbs > 1000*eps && deltaRel > 1000*eps)
                throw std::logic_error("oops: derivative @("+std::to_string((long double) x)+","+std::to_string((long double) y)+"): "
                                       + std::to_string((long double) zPrime) + " vs "
                                       + std::to_string((long double) zEval.derivative(0))
                                       + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivative(0))));

        }
    }
}

// prototypes
static double myScalarMin(double a, double b)
{ return std::min(a, b); }

static double myScalarMax(double a, double b)
{ return std::max(a, b); }

inline void testAll()
{
    // the following is commented out because it is supposed to produce a compiler
    // error. This is the case since the function does not calculate the derivatives
    // w.r.t. Pressure but they have been requested...
    //const auto& result2 = Opm::DenseAd::sqrt(TemperatureEval::createVariable<Pressure>(4.0));

    std::cout << "testing operators and constructors\n";
    const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e3;
    testOperators(eps);

    std::cout << "testing min()\n";
    test2DFunction1(Opm::DenseAd::min<Scalar, numVars>,
                            myScalarMin,
                            -1000, 1000,
                            /*p=*/1.234);

    test2DFunction2(Opm::DenseAd::min<Scalar, numVars>,
                            myScalarMin,
                            /*T=*/1.234,
                            -1000, 1000);

    std::cout << "testing max()\n";
    test2DFunction1(Opm::DenseAd::max<Scalar, numVars>,
                            myScalarMax,
                            -1000, 1000,
                            /*p=*/1.234);

    test2DFunction2(Opm::DenseAd::max<Scalar, numVars>,
                            myScalarMax,
                            /*T=*/1.234,
                            -1000, 1000);

    std::cout << "testing pow()\n";
    testPowBase();
    testPowExp();

    std::cout << "testing abs()\n";
    test1DFunction(Opm::DenseAd::abs<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::abs));

    std::cout << "testing sqrt()\n";
    test1DFunction(Opm::DenseAd::sqrt<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::sqrt));

    std::cout << "testing sin()\n";
    test1DFunction(Opm::DenseAd::sin<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::sin),
                           0, 2*M_PI);

    std::cout << "testing asin()\n";
    test1DFunction(Opm::DenseAd::asin<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::asin),
                           -1.0, 1.0);

    std::cout << "testing cos()\n";
    test1DFunction(Opm::DenseAd::cos<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::cos),
                           0, 2*M_PI);

    std::cout << "testing acos()\n";
    test1DFunction(Opm::DenseAd::acos<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::acos),
                           -1.0, 1.0);

    std::cout << "testing tan()\n";
    test1DFunction(Opm::DenseAd::tan<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::tan),
                           -M_PI / 2 * 0.95, M_PI / 2 * 0.95);

    std::cout << "testing atan()\n";
    test1DFunction(Opm::DenseAd::atan<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::atan),
                           -10*1000.0, 10*1000.0);

    std::cout << "testing atan2()\n";
    testAtan2();

    std::cout << "testing exp()\n";
    test1DFunction(Opm::DenseAd::exp<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::exp),
                           -100, 100);

    std::cout << "testing log()\n";
    test1DFunction(Opm::DenseAd::log<Scalar, numVars>,
                           static_cast<Scalar (*)(Scalar)>(std::log),
                           1e-6, 1e9);

    while (false) {
        Scalar val1 OPM_UNUSED = 0.0;
        Scalar val2 OPM_UNUSED = 1.0;
        Scalar resultVal OPM_UNUSED;
        typedef Opm::DenseAd::Evaluation<Scalar, numVars> TmpEval;
        TmpEval eval1 OPM_UNUSED = 1.0;
        TmpEval eval2 OPM_UNUSED = 2.0;
        TmpEval resultEval OPM_UNUSED;

        // make sure that the convenince functions work (i.e., that everything can be
        // accessed without the MathToolbox<Scalar> detour.)
        resultVal = Opm::constant<Scalar>(val1);
        resultVal = Opm::variable<Scalar>(val1, /*idx=*/0);
        resultVal = Opm::decay<Scalar>(val1);
        resultVal = Opm::scalarValue(val1);
        resultVal = Opm::getValue(val1);
        resultVal = Opm::min(val1, val2);
        resultVal = Opm::max(val1, val2);
        resultVal = Opm::atan2(val1, val2);
        resultVal = Opm::pow(val1, val2);
        resultVal = Opm::abs(val1);
        resultVal = Opm::atan(val1);
        resultVal = Opm::sin(val1);
        resultVal = Opm::asin(val1);
        resultVal = Opm::cos(val1);
        resultVal = Opm::acos(val1);
        resultVal = Opm::sqrt(val1);
        resultVal = Opm::exp(val1);
        resultVal = Opm::log(val1);

        resultEval = Opm::constant<TmpEval>(val1);
        resultEval = Opm::variable<TmpEval>(val1, /*idx=*/0);
        resultEval = Opm::decay<TmpEval>(eval1);
        resultVal = Opm::decay<Scalar>(eval1);
        resultVal = Opm::scalarValue(eval1);
        resultVal = Opm::getValue(eval1);
        resultEval = Opm::min(eval1, eval2);
        resultEval = Opm::min(eval1, val2);
        resultEval = Opm::max(eval1, eval2);
        resultEval = Opm::max(eval1, val2);
        resultEval = Opm::atan2(eval1, eval2);
        resultEval = Opm::atan2(eval1, val2);
        resultEval = Opm::pow(eval1, eval2);
        resultEval = Opm::pow(eval1, val2);
        resultEval = Opm::abs(eval1);
        resultEval = Opm::atan(eval1);
        resultEval = Opm::sin(eval1);
        resultEval = Opm::asin(eval1);
        resultEval = Opm::cos(eval1);
        resultEval = Opm::acos(eval1);
        resultEval = Opm::sqrt(eval1);
        resultEval = Opm::exp(eval1);
        resultEval = Opm::log(eval1);
    }
}

};//TestEnv


int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    TestEnv<double, 15>().testAll();
    TestEnv<double, 2>().testAll();
    TestEnv<float, 15>().testAll();
    TestEnv<float, 2>().testAll();

    return 0;
}
