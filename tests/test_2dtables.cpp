/*
  Copyright (C) 2014 by Andreas Lauser

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
 * \brief This is the unit test for the 2D tabulation classes.
 *
 * I.e., for the UniformTabulated2DFunction and UniformXTabulated2DFunction classes.
 */
#include "config.h"

#include <opm/material/UniformXTabulated2DFunction.hpp>
#include <opm/material/UniformTabulated2DFunction.hpp>

#include <memory>
#include <cmath>
#include <iostream>

typedef double Scalar;

Scalar testFn1(Scalar x, Scalar y)
{ return x; }

Scalar testFn2(Scalar x, Scalar y)
{ return y; }

Scalar testFn3(Scalar x, Scalar y)
{ return x*y; }

template <class Fn>
std::shared_ptr<Opm::UniformTabulated2DFunction<Scalar> >
createUniformTabulatedFunction(Fn &f)
{
    Scalar xMin = -2.0;
    Scalar xMax = 3.0;
    Scalar m = 50;

    Scalar yMin = -1/2.0;
    Scalar yMax = 1/3.0;
    Scalar n = 40;

    auto tab = std::make_shared<Opm::UniformTabulated2DFunction<Scalar>>(
        xMin, xMax, m,
        yMin, yMax, n);
    for (int i = 0; i < m; ++i) {
        Scalar x = xMin + Scalar(i)/(m - 1) * (xMax - xMin);
        for (int j = 0; j < n; ++j) {
            Scalar y = yMin + Scalar(j)/(n - 1) * (yMax - yMin);
            tab->setSamplePoint(i, j, f(x, y));
        }
    }

    return tab;
}

template <class Fn>
std::shared_ptr<Opm::UniformXTabulated2DFunction<Scalar> >
createUniformXTabulatedFunction(Fn &f)
{
    Scalar xMin = -2.0;
    Scalar xMax = 3.0;
    Scalar m = 50;

    Scalar yMin = -1/2.0;
    Scalar yMax = 1/3.0;
    Scalar n = 40;

    auto tab = std::make_shared<Opm::UniformXTabulated2DFunction<Scalar>>();
    for (int i = 0; i < m; ++i) {
        Scalar x = xMin + Scalar(i)/(m - 1) * (xMax - xMin);
        tab->appendXPos(x);
        for (int j = 0; j < n; ++j) {
            Scalar y = yMin + Scalar(j)/(n -1) * (yMax - yMin);
            tab->appendSamplePoint(i, y, f(x, y));
        }
    }

    return tab;
}


template <class Fn>
std::shared_ptr<Opm::UniformXTabulated2DFunction<Scalar> >
createUniformXTabulatedFunction2(Fn &f)
{
    Scalar xMin = -2.0;
    Scalar xMax = 3.0;
    Scalar m = 50;


    auto tab = std::make_shared<Opm::UniformXTabulated2DFunction<Scalar>>();
    for (int i = 0; i < m; ++i) {
        Scalar x = xMin + Scalar(i)/(m - 1) * (xMax - xMin);
        tab->appendXPos(x);

        Scalar n = i + 10;
        Scalar yMin = - (x + 1);
        Scalar yMax = (x + 1);

        for (int j = 0; j < n; ++j) {
            Scalar y = yMin + Scalar(j)/(n -1) * (yMax - yMin);
            tab->appendSamplePoint(i, y, f(x, y));
        }
    }

    return tab;
}

template <class Fn, class Table>
bool compareTableWithAnalyticFn(const Table &table,
                                Scalar xMin,
                                Scalar xMax,
                                int numX,

                                Scalar yMin,
                                Scalar yMax,
                                int numY,

                                Fn &f,
                                Scalar tolerance = 1e-8)
{
    // make sure that the tabulated function evaluates to the same thing as the analytic
    // one (modulo tolerance)
    for (int i = 1; i <= numX; ++i) {
        Scalar x = xMin + Scalar(i)/numX*(xMax - xMin);

        for (int j = 0; j < numY; ++j) {
            Scalar y = yMin + Scalar(j)/numY*(yMax - yMin);
            if (std::abs(table->eval(x, y) - f(x, y)) > tolerance) {
                std::cerr << __FILE__ << ":" << __LINE__ << ": table->eval("<<x<<","<<y<<") != f("<<x<<","<<y<<"): " << table->eval(x,y) << " != " << f(x,y) << "\n";
                return false;
            }
        }
    }

    return true;
}

template <class UniformTablePtr, class UniformXTablePtr, class Fn>
bool compareTables(const UniformTablePtr uTable,
                   const UniformXTablePtr uXTable,
                   Fn &f,
                   Scalar tolerance = 1e-8)
{
    // make sure the uniform and the non-uniform tables exhibit the same dimensions
    if (std::abs(uTable->xMin() - uXTable->xMin()) > 1e-8) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->xMin() != uXTable->xMin(): " << uTable->xMin() << " != " << uXTable->xMin() << "\n";
        return false;
    }
    if (std::abs(uTable->xMax() - uXTable->xMax()) > 1e-8) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->xMax() != uXTable->xMax(): " << uTable->xMax() << " != " << uXTable->xMax() << "\n";
        return false;
    }
    if (uTable->numX() != uXTable->numX()) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->numX() != uXTable->numX(): " << uTable->numX() << " != " << uXTable->numX() << "\n";
        return false;
    }

    for (int i = 0; i < uTable->numX(); ++i) {
        if (std::abs(uTable->yMin() - uXTable->yMin(i)) > 1e-8) {
            std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->yMin() != uXTable->yMin("<<i<<"): " << uTable->yMin() << " != " << uXTable->yMin(i) << "\n";
            return false;
        }

        if (std::abs(uTable->yMax() - uXTable->yMax(i)) > 1e-8) {
            std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->yMax() != uXTable->yMax("<<i<<"): " << uTable->yMax() << " != " << uXTable->yMax(i) << "\n";
            return false;
        }

        if (uTable->numY() != uXTable->numY(i)) {
            std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->numY() != uXTable->numY("<<i<<"): " << uTable->numY() << " != " << uXTable->numY(i) << "\n";
            return false;
        }
    }

    // make sure that the x and y values are identical
    for (int i = 0; i < uTable->numX(); ++i) {
        if (std::abs(uTable->iToX(i) - uXTable->iToX(i)) > 1e-8) {
            std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->iToX("<<i<<") != uXTable->iToX("<<i<<"): " << uTable->iToX(i) << " != " << uXTable->iToX(i) << "\n";
            return false;
        }

        for (int j = 0; j < uTable->numY(); ++j) {
            if (std::abs(uTable->jToY(j) - uXTable->jToY(i, j)) > 1e-8) {
                std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->jToY("<<j<<") != uXTable->jToY("<<i<<","<<j<<"): " << uTable->jToY(i) << " != " << uXTable->jToY(i, j) << "\n";
                return false;
            }
        }
    }

    // check that the appicable range is correct. Note that due to rounding errors it is
    // undefined whether the table applies to the boundary of the tabulated domain or not
    Scalar xMin = uTable->xMin();
    Scalar yMin = uTable->yMin();
    Scalar xMax = uTable->xMax();
    Scalar yMax = uTable->yMax();

    Scalar x = xMin - 1e-8;
    Scalar y = yMin - 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMin - 1e-8;
    y = yMin + 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMin + 1e-8;
    y = yMin - 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMin + 1e-8;
    y = yMin + 1e-8;
    if (!uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": !uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (!uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": !uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMax + 1e-8;
    y = yMax + 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMax - 1e-8;
    y = yMax + 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMax + 1e-8;
    y = yMax - 1e-8;
    if (uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    x = xMax - 1e-8;
    y = yMax - 1e-8;
    if (!uTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": !uTable->applies("<<x<<","<<y<<")\n";
        return false;
    }
    if (!uXTable->applies(x, y)) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": !uXTable->applies("<<x<<","<<y<<")\n";
        return false;
    }

    // make sure that the function values at the sampling points are identical and that
    // they correspond to the analytic function
    int m2 = uTable->numX()*5;
    int n2 = uTable->numY()*5;
    if (!compareTableWithAnalyticFn(uTable,
                                    xMin, xMax, m2,
                                    yMin, yMax, n2,
                                    f,
                                    tolerance))
        return false;
    if (!compareTableWithAnalyticFn(uXTable,
                                    xMin, xMax, m2,
                                    yMin, yMax, n2,
                                    f,
                                    tolerance))
        return false;

    return true;
}

int main()
{
    auto uniformTab = createUniformTabulatedFunction(testFn1);
    auto uniformXTab = createUniformXTabulatedFunction(testFn1);
    if (!compareTables(uniformTab, uniformXTab, testFn1, /*tolerance=*/1e-12))
        return 1;

    uniformTab = createUniformTabulatedFunction(testFn2);
    uniformXTab = createUniformXTabulatedFunction(testFn2);
    if (!compareTables(uniformTab, uniformXTab, testFn2, /*tolerance=*/1e-12))
        return 1;

    uniformTab = createUniformTabulatedFunction(testFn3);
    uniformXTab = createUniformXTabulatedFunction(testFn3);
    if (!compareTables(uniformTab, uniformXTab, testFn3, /*tolerance=*/1e-2))
        return 1;

    uniformXTab = createUniformXTabulatedFunction2(testFn3);
    if (!compareTableWithAnalyticFn(uniformXTab,
                                    -10, 10, 100,
                                    -10, 10, 100,
                                    testFn3,
                                    /*tolerance=*/1e-2))
        return 1;

    // CSV output for debugging
#if 0
    int m = 100;
    int n = 100;
    Scalar xMin = -3.0;
    Scalar xMax = 4.0;

    Scalar yMin = -1;
    Scalar yMax = 1;
    for (int i = 0; i < m; ++i) {
        Scalar x = xMin + Scalar(i)/m * (xMax - xMin);

        for (int j = 0; j < n; ++j) {
            Scalar y = yMin + Scalar(j)/n * (yMax - yMin);

            std::cout << x << " "
                      << y << " "
                      << uniformXTab->eval(x,y,true) << "\n";
        }
        std::cout << "\n";
    }
#endif

    return 0;
}
