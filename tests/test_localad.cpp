// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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

#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <opm/material/common/Unused.hpp>

#include <opm/material/localad/Evaluation.hpp>
#include <opm/material/localad/Math.hpp>

struct TestVariables
{
    static const int size = 3;

    static const int temperatureIdx = 0;
    static const int pressureIdx = 1;
    static const int saturationIdx = 2;
};

template <class Scalar>
struct Tolerance
{
    static constexpr Scalar eps = 1e-10;
};

template <>
struct Tolerance< float >
{
    static constexpr  float eps = 1e-6;
};


template <class Scalar, class VariablesDescriptor>
void testOperators(const Scalar tolerance )
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    // test the constructors of the Opm::LocalAd::Evaluation class
    const Scalar c = 1.234;
    const Scalar x = 4.567;
    const Scalar y = 8.910;
    const Eval cEval = Eval::createConstant(c);
    const Eval OPM_UNUSED c2Eval = c;
    const Eval xEval = Eval::createVariable(x, 0);
    const Eval yEval = Eval::createVariable(y, 1);

    // test the non-inplace operators
    {
        Eval a = xEval + yEval;
        if (std::abs(a.value - (x + y)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval b = xEval + c;
        if (std::abs(b.value - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval d = xEval + cEval;
        if (std::abs(d.value - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");
    }

    {
        Eval a = xEval - yEval;
        if (std::abs(a.value - (x - y)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval b = xEval - c;
        if (std::abs(b.value - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval d = xEval - cEval;
        if (std::abs(d.value - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");
    }

    {
        Eval a = xEval*yEval;
        if (std::abs(a.value - (x*y)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval b = xEval*c;
        if (std::abs(b.value - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval d = xEval*cEval;
        if (std::abs(d.value - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");
    }

    {
        Eval a = xEval/yEval;
        if (std::abs(a.value - (x/y)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval b = xEval/c;
        if (std::abs(b.value - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval d = xEval/cEval;
        if (std::abs(d.value - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");
    }

    // test the inplace operators
    {
        Eval a = xEval;
        a += yEval;
        if (std::abs(a.value - (x + y)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval b = xEval;
        b += c;
        if (std::abs(b.value - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");

        Eval d = xEval;
        d += cEval;
        if (std::abs(d.value - (x + c)) > tolerance)
            throw std::logic_error("oops: operator+");
    }

    {
        Eval a = xEval;
        a -= yEval;
        if (std::abs(a.value - (x - y)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval b = xEval;
        b -= c;
        if (std::abs(b.value - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");

        Eval d = xEval;
        d -= cEval;
        if (std::abs(d.value - (x - c)) > tolerance)
            throw std::logic_error("oops: operator-");
    }

    {
        Eval a = xEval;
        a *= yEval;
        if (std::abs(a.value - (x*y)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval b = xEval;
        b *= c;
        if (std::abs(b.value - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");

        Eval d = xEval;
        d *= cEval;
        if (std::abs(d.value - (x*c)) > tolerance)
            throw std::logic_error("oops: operator*");
    }

    {
        Eval a = xEval;
        a /= yEval;
        if (std::abs(a.value - (x/y)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval b = xEval;
        b /= c;
        if (std::abs(b.value - (x/c)) > tolerance)
            throw std::logic_error("oops: operator/");

        Eval d = xEval;
        d /= cEval;
        if (std::abs(d.value - (x/c)) > tolerance)
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

template <class Scalar, class VariablesDescriptor, class AdFn, class ClassicFn>
void test1DFunction(AdFn* adFn, ClassicFn* classicFn, Scalar xMin = 1e-6, Scalar xMax = 1000)
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar x = Scalar(i)/(n - 1)*(xMax - xMin) + xMin;

        const auto& xEval = Eval::createVariable(x, 0);
        const Eval& yEval = adFn(xEval);

        const Scalar eps = Tolerance< Scalar > :: eps;
        Scalar y = classicFn(x);
        Scalar yStar1 = classicFn(x - eps);
        Scalar yStar2 = classicFn(x + eps);
        Scalar yPrime = (yStar2 - yStar1)/(2*eps);

        if (std::abs(y-yEval.value) > 5e-14)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(yPrime - yEval.derivatives[0]);
        Scalar deltaRel = std::abs(deltaAbs/yPrime);
        if (deltaAbs > 1e-3 && deltaRel > 1e-3)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) yPrime) + " vs "
                                   + std::to_string((long double) yEval.derivatives[0])
                                   + " delta: " + std::to_string((long double) std::abs(yPrime - yEval.derivatives[0])));
    }
}

template <class Scalar,
          class VariablesDescriptor,
          class AdFn,
          class ClassicFn>
void test2DFunction1(AdFn* adFn, ClassicFn* classicFn, Scalar xMin, Scalar xMax, Scalar y)
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar x = Scalar(i)/(n - 1)*(xMax - xMin) + xMin;

        const auto& xEval = Eval::createVariable(x, 0);
        const auto& yEval = Eval::createConstant(y);
        const Eval& zEval = adFn(xEval, yEval);

        const Scalar eps = Tolerance< Scalar > :: eps;
        Scalar z = classicFn(x, y);
        Scalar zStar1 = classicFn(x - eps, y);
        Scalar zStar2 = classicFn(x + eps, y);
        Scalar zPrime = (zStar2 - zStar1)/(2.*eps);

        if (z != zEval.value)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval.derivatives[0]);
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1e-3 && deltaRel > 1e-3)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval.derivatives[0])
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivatives[0])));
    }
}

template <class Scalar,
          class VariablesDescriptor,
          class AdFn,
          class ClassicFn>
void test2DFunction2(AdFn* adFn, ClassicFn* classicFn, Scalar x, Scalar yMin, Scalar yMax)
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar y = Scalar(i)/(n - 1)*(yMax - yMin) + yMin;

        const auto& xEval = Eval::createConstant(x);
        const auto& yEval = Eval::createVariable(y, 1);
        const Eval& zEval = adFn(xEval, yEval);

        const Scalar eps = Tolerance< Scalar > :: eps;
        Scalar z = classicFn(x, y);
        Scalar zStar1 = classicFn(x, y - eps);
        Scalar zStar2 = classicFn(x, y + eps);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (z != zEval.value)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval.derivatives[1]);
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1e-3 && deltaRel > 1e-3)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) x)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval.derivatives[1])
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivatives[1])));
    }
}

template <class Scalar, class VariablesDescriptor>
void testPowBase(Scalar baseMin = 1e-2, Scalar baseMax = 100)
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    Scalar exp = 1.234;
    const auto& expEval = Eval::createConstant(exp);

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar base = Scalar(i)/(n - 1)*(baseMax - baseMin) + baseMin;

        const auto& baseEval = Eval::createVariable(base, 0);
        const Eval& zEval1 = pow(baseEval, exp);
        const Eval& zEval2 = pow(baseEval, expEval);

        const Scalar eps = 1e-5;
        Scalar z = pow(base, exp);
        Scalar zStar1 = pow(base - eps, exp);
        Scalar zStar2 = pow(base + eps, exp);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (z != zEval2.value)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval1.derivatives[0]);
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1e-3 && deltaRel > 1e-3)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) base)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval1.derivatives[0])
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval1.derivatives[0])));

        if (!zEval1.isSame(zEval2, /*tolerance=*/1e-9))
            throw std::logic_error("oops: pow(Eval, Scalar) != pow(Eval, Eval)");
    }
}

template <class Scalar, class VariablesDescriptor>
void testPowExp(Scalar expMin = -100, Scalar expMax = 100)
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

    Scalar base = 1.234;
    const auto& baseEval = Eval::createConstant(base);

    int n = 100*1000;
    for (int i = 0; i < n; ++ i) {
        Scalar exp = Scalar(i)/(n - 1)*(expMax - expMin) + expMin;
        const auto& expEval = Eval::createVariable(exp, 1);

        const Eval& zEval1 = pow(base, expEval);
        const Eval& zEval2 = pow(baseEval, expEval);

        const Scalar eps = 1e-8;
        Scalar z = pow(base, exp);
        Scalar zStar1 = pow(base, exp - eps);
        Scalar zStar2 = pow(base, exp + eps);
        Scalar zPrime = (zStar2 - zStar1)/(2*eps);

        if (z != zEval2.value)
            throw std::logic_error("oops: value");

        Scalar deltaAbs = std::abs(zPrime - zEval1.derivatives[1]);
        Scalar deltaRel = std::abs(deltaAbs/zPrime);
        if (deltaAbs > 1e-3 && deltaRel > 1e-3)
            throw std::logic_error("oops: derivative @"+std::to_string((long double) base)+": "
                                   + std::to_string((long double) zPrime) + " vs "
                                   + std::to_string((long double) zEval1.derivatives[1])
                                   + " delta: " + std::to_string((long double) std::abs(zPrime - zEval1.derivatives[1])));

        if (!zEval1.isSame(zEval2, /*tolerance=*/1e-5))
            throw std::logic_error("oops: pow(Eval, Scalar) != pow(Eval, Eval)");
    }
}

template <class Scalar, class VariablesDescriptor>
void testAtan2()
{
    typedef Opm::LocalAd::Evaluation<Scalar, VariablesDescriptor, VariablesDescriptor::size> Eval;

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

            const Scalar eps = 1e-8;
            Scalar z = atan2(x, y);
            Scalar zStar1 = atan2(x - eps, y - eps);
            Scalar zStar2 = atan2(x + eps, y + eps);
            Scalar zPrime = (zStar2 - zStar1)/(2*eps);

            if (z != zEval.value)
                throw std::logic_error("oops: value");

            Scalar deltaAbs = std::abs(zPrime - zEval.derivatives[0]);
            Scalar deltaRel = std::abs(deltaAbs/zPrime);
            if (deltaAbs > 1e-3 && deltaRel > 1e-3)
                throw std::logic_error("oops: derivative @("+std::to_string((long double) x)+","+std::to_string((long double) y)+"): "
                                       + std::to_string((long double) zPrime) + " vs "
                                       + std::to_string((long double) zEval.derivatives[0])
                                       + " delta: " + std::to_string((long double) std::abs(zPrime - zEval.derivatives[0])));

        }
    }
}

// prototypes
double myScalarMin(double a, double b);
double myScalarMax(double a, double b);

double myScalarMin(double a, double b)
{ return std::min(a, b); }

double myScalarMax(double a, double b)
{ return std::max(a, b); }

template <class Scalar>
inline void testAll()
{
    typedef TestVariables VarsDescriptor;

    // the following is commented out because it is supposed to produce a compiler
    // error. This is the case since the function does not calculate the derivatives
    // w.r.t. Pressure but they have been requested...
    //const auto& result2 = Opm::LocalAd::sqrt(TemperatureEval::createVariable<Pressure>(4.0));

    std::cout << "testing operators and constructors\n";
    testOperators<Scalar, VarsDescriptor>( Tolerance< Scalar > :: eps );

    std::cout << "testing min()\n";
    test2DFunction1<Scalar, VarsDescriptor>(Opm::LocalAd::min<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                            myScalarMin,
                                            -1000, 1000,
                                            /*p=*/1.234);

    test2DFunction2<Scalar, VarsDescriptor>(Opm::LocalAd::min<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                            myScalarMin,
                                            /*T=*/1.234,
                                            -1000, 1000);

    std::cout << "testing max()\n";
    test2DFunction1<Scalar, VarsDescriptor>(Opm::LocalAd::max<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                            myScalarMax,
                                            -1000, 1000,
                                            /*p=*/1.234);

    test2DFunction2<Scalar, VarsDescriptor>(Opm::LocalAd::max<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                            myScalarMax,
                                            /*T=*/1.234,
                                            -1000, 1000);

    std::cout << "testing pow()\n";
    testPowBase<Scalar, VarsDescriptor>();
    testPowExp<Scalar, VarsDescriptor>();

    std::cout << "testing abs()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::abs<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::abs));

    std::cout << "testing sqrt()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::sqrt<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::sqrt));

    std::cout << "testing sin()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::sin<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::sin),
                                           0, 2*M_PI);

    std::cout << "testing asin()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::asin<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::asin),
                                           -1.0, 1.0);

    std::cout << "testing cos()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::cos<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::cos),
                                           0, 2*M_PI);

    std::cout << "testing acos()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::acos<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::acos),
                                           -1.0, 1.0);

    std::cout << "testing tan()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::tan<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::tan),
                                           -M_PI / 2 * 0.95, M_PI / 2 * 0.95);

    std::cout << "testing atan()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::atan<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::atan),
                                           -10*1000.0, 10*1000.0);

    std::cout << "testing atan2()\n";
    testAtan2<Scalar, VarsDescriptor>();

    std::cout << "testing exp()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::exp<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::exp),
                                           -100, 100);

    std::cout << "testing log()\n";
    test1DFunction<Scalar, VarsDescriptor>(Opm::LocalAd::log<Scalar, VarsDescriptor, VarsDescriptor::size>,
                                           static_cast<Scalar (*)(Scalar)>(std::log),
                                           1e-6, 1e9);
}

int main()
{
    testAll< double >();
    // testAll< float  >();
    return 0;
}
