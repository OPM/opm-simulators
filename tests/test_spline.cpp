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
 * \brief This is a program to test the polynomial spline interpolation.
 *
 * It just prints some function to stdout. You can look at the result
 * using the following commands:
 *
----------- snip -----------
./test_spline > spline.csv
gnuplot
gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
              "spline.csv" using 1:3 w l ti "Derivative", \
              "spline.csv" using 1:4 w p ti "Monotonical"
----------- snap -----------
*/
#include "config.h"

#include <opm/material/common/Spline.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>

template <class Spline>
void testCommon(const Spline& sp,
                const double* x,
                const double* y)
{
    static double eps = 1e-10;
    static double epsFD = 1e-7;
    size_t n = sp.numSamples();
    for (size_t i = 0; i < n; ++i) {
        // sure that we hit all sampling points
        double y0 = (i>0)?sp.eval(x[i]-eps):y[0];
        double y1 = sp.eval(x[i]);
        double y2 = (i<n-1)?sp.eval(x[i]+eps):y[n-1];
        if (std::abs(y0 - y[i]) > 100*eps || std::abs(y2 - y[i]) > 100*eps)
            throw std::runtime_error("Spline seems to be discontinuous at sampling point "+std::to_string(i)+"!");
        if (std::abs(y1 - y[i]) > eps)
            throw std::runtime_error("Spline does not capture sampling point "+std::to_string(i)+"!");
        // make sure the derivative is continuous (assuming that the
        // second derivative is smaller than 1000)
        double d1 = sp.evalDerivative(x[i]);
        double d0 = (i>0)?sp.evalDerivative(x[i]-eps):d1;
        double d2 = (i<n-1)?sp.evalDerivative(x[i]+eps):d1;
        if (std::abs(d1 - d0) > 1000*eps || std::abs(d2 - d0) > 1000*eps)
            throw std::runtime_error("Spline seems to exhibit a discontinuous derivative at sampling point "+std::to_string(i)+"!");
    }
    // make sure the derivatives are consistent with the curve
    size_t np = 3*n;
    for (size_t i = 0; i < np; ++i) {
        double xMin = sp.xAt(0);
        double xMax = sp.xAt(sp.numSamples() - 1);
        double xval = xMin + (xMax - xMin)*i/np;
        // first derivative
        double y1 = sp.eval(xval+epsFD);
        double y0 = sp.eval(xval);
        double mFD = (y1 - y0)/epsFD;
        double m = sp.evalDerivative(xval);
        if (std::abs( mFD - m ) > 1000*epsFD)
            throw std::runtime_error("Derivative of spline seems to be inconsistent with cuve"
                                     " ("+std::to_string(mFD)+" - "+std::to_string(m)+" = "
                                     +std::to_string(mFD - m)+")!");
        // second derivative
        y1 = sp.evalDerivative(xval+epsFD);
        y0 = sp.evalDerivative(xval);
        mFD = (y1 - y0)/epsFD;
        m = sp.evalSecondDerivative(xval);
        if (std::abs( mFD - m ) > 1000*epsFD)
            throw std::runtime_error("Second derivative of spline seems to be inconsistent with cuve"
                                     " ("+std::to_string(mFD)+" - "+std::to_string(m)+" = "+std::to_string(mFD - m)+")!");
        // Third derivative
        y1 = sp.evalSecondDerivative(xval+epsFD);
        y0 = sp.evalSecondDerivative(xval);
        mFD = (y1 - y0)/epsFD;
        m = sp.evalThirdDerivative(xval);
        if (std::abs( mFD - m ) > 1000*epsFD)
            throw std::runtime_error("Third derivative of spline seems to be inconsistent with cuve"
                                     " ("+std::to_string(mFD)+" - "+std::to_string(m)+" = "+std::to_string(mFD - m)+")!");
    }
}

template <class Spline>
void testFull(const Spline& sp,
              const double* x,
              const double* y,
              double m0,
              double m1)
{
    // test the common properties of splines
    testCommon(sp, x, y);
    static double eps = 1e-5;
    size_t n = sp.numSamples();
    // make sure the derivative at both end points is correct
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[n-1]);
    if (std::abs(d0 - m0) > eps)
        throw std::runtime_error("Invalid derivative at beginning of interval: is "
                                 +std::to_string(d0)+" ought to be "+std::to_string(m0));
    if (std::abs(d1 - m1) > eps)
        throw std::runtime_error("Invalid derivative at end of interval: is "
                                 +std::to_string(d1)+" ought to be "+std::to_string(m1));
}

template <class Spline>
void testNatural(const Spline& sp,
                 const double* x,
                 const double* y)
{
    // test the common properties of splines
    testCommon(sp, x, y);
    static double eps = 1e-5;
    size_t n = sp.numSamples();
    // make sure the second derivatives at both end points are 0
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[0] + eps);
    double d2 = sp.evalDerivative(x[n-1] - eps);
    double d3 = sp.evalDerivative(x[n-1]);
    if (std::abs(d1 - d0)/eps > 1000*eps)
        throw std::runtime_error("Invalid second derivative at beginning of interval: is "
                                 +std::to_string((d1 - d0)/eps)+" ought to be 0");
    if (std::abs(d3 - d2)/eps > 1000*eps)
        throw std::runtime_error("Invalid second derivative at end of interval: is "
                                 +std::to_string((d3 - d2)/eps)+" ought to be 0");
}

template <class Spline>
void testMonotonic(const Spline& sp,
                   const double* x,
                   const double* y)
{
    // test the common properties of splines
    testCommon(sp, x, y);
    size_t n = sp.numSamples();
    for (size_t i = 0; i < n - 1; ++ i) {
        // make sure that the spline is monotonic for each interval
        // between sampling points
        if (!sp.monotonic(x[i], x[i + 1]))
            throw std::runtime_error("Spline says it is not monotonic in interval "
                                     +std::to_string(i)+" where it should be");
        // test the intersection methods
        double d = (y[i] + y[i+1])/2;
        double interX = sp.template intersectInterval<double>(x[i], x[i+1],
                                                              /*a=*/0, /*b=*/0, /*c=*/0, d);
        double interY = sp.eval(interX);
        if (std::abs(interY - d) > 1e-5)
            throw std::runtime_error("Spline::intersectInterval() seems to be broken: "
                                     +std::to_string(sp.eval(interX))+" - "+std::to_string(d)+
                                     " = "+std::to_string(sp.eval(interX) - d)+"!");
    }
    // make sure the spline says to be monotonic on the (extrapolated)
    // left and right sides
    if (!sp.monotonic(x[0] - 1.0, (x[0] + x[1])/2, /*extrapolate=*/true))
        throw std::runtime_error("Spline says it is not monotonic on left side where it should be");
    if (!sp.monotonic((x[n - 2]+ x[n - 1])/2, x[n-1] + 1.0, /*extrapolate=*/true))
        throw std::runtime_error("Spline says it is not monotonic on right side where it should be");
    for (size_t i = 0; i < n - 2; ++ i) {
        // make sure that the spline says that it is non-monotonic for
        // if extrema are within the queried interval
        if (sp.monotonic((x[i] + x[i + 1])/2, (x[i + 1] + x[i + 2])/2))
            throw std::runtime_error("Spline says it is monotonic in interval "
                                     +std::to_string(i)+" where it should not be");
    }
}

// function prototype to prevent some compilers producing a warning
void testAll();
void testAll()
{
    double x[] = { 0, 5, 7.5, 8.75, 9.375 };
    double y[] = { 10, 0, 10, 0, 10 };
    double m0 = 10;
    double m1 = -10;
    double points[][2] =
        {
            {x[0], y[0]},
            {x[1], y[1]},
            {x[2], y[2]},
            {x[3], y[3]},
            {x[4], y[4]},
        };

    std::initializer_list<const std::pair<double, double> > pointsInitList =
        {
            {x[0], y[0]},
            {x[1], y[1]},
            {x[2], y[2]},
            {x[3], y[3]},
            {x[4], y[4]},
        };

    std::vector<double> xVec;
    std::vector<double> yVec;
    std::vector<double*> pointVec;
    for (int i = 0; i < 5; ++i) {
        xVec.push_back(x[i]);
        yVec.push_back(y[i]);
        pointVec.push_back(points[i]);
    }
    /////////
    // test spline with two sampling points
    /////////
    // full spline
    { Opm::Spline<double> sp(x[0], x[1], y[0], y[1], m0, m1); sp.set(x[0],x[1],y[0],y[1],m0, m1); testFull(sp, x, y, m0, m1); };
    { Opm::Spline<double> sp(2, x, y, m0, m1); sp.setXYArrays(2, x, y, m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp(static_cast<size_t>(2), points, m0, m1); sp.setArrayOfPoints(2, points, m0, m1); testFull(sp, x, y, m0, m1); };
    /////////
    // test variable length splines
    /////////
    // full spline
    { Opm::Spline<double> sp(5, x, y, m0, m1); sp.setXYArrays(5,x,y,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp(xVec, yVec, m0, m1); sp.setXYContainers(xVec,yVec,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp; sp.setArrayOfPoints(5,points,m0, m1); testFull(sp, x, y, m0, m1); };
    { Opm::Spline<double> sp; sp.setContainerOfPoints(pointVec,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp; sp.setContainerOfTuples(pointsInitList,m0, m1); testFull(sp, x, y, m0, m1); };
    // natural spline
    { Opm::Spline<double> sp(5, x, y); sp.setXYArrays(5,x,y); testNatural(sp, x, y);  };
    { Opm::Spline<double> sp(xVec, yVec); sp.setXYContainers(xVec,yVec); testNatural(sp, x, y); };
    { Opm::Spline<double> sp; sp.setArrayOfPoints(5,points); testNatural(sp, x, y); };
    { Opm::Spline<double> sp; sp.setContainerOfPoints(pointVec); testNatural(sp, x, y); };
    { Opm::Spline<double> sp; sp.setContainerOfTuples(pointsInitList); testNatural(sp, x, y); };
}

// function prototype to prevent some compilers producing a warning
void plot();
void plot()
{
    const int numSamples = 5;
    const int n = numSamples - 1;
    typedef std::array<double, numSamples> FV;
    double x_[] = { 0, 5, 7.5, 8.75, 10 };
    double y_[] = { 10, 0, 10, 0, 10 };
    double m1 = 10;
    double m2 = -10;
    FV& xs = *reinterpret_cast<FV*>(x_);
    FV& ys = *reinterpret_cast<FV*>(y_);
    Opm::Spline<double> spFull(xs, ys, m1, m2);
    Opm::Spline<double> spNatural(xs, ys);
    Opm::Spline<double> spPeriodic(xs, ys, /*type=*/Opm::Spline<double>::Periodic);
    Opm::Spline<double> spMonotonic(xs, ys, /*type=*/Opm::Spline<double>::Monotonic);
    testMonotonic(spMonotonic, x_, y_);
    spFull.printCSV(x_[0] - 1.00001,
                    x_[n] + 1.00001,
                    1000);
    std::cout << "\n";
    spNatural.printCSV(x_[0] - 1.00001,
                       x_[n] + 1.00001,
                       1000);
    std::cout << "\n";
    spPeriodic.printCSV(x_[0] - 1.00001,
                        x_[n] + 1.00001,
                       1000);
    std::cout << "\n";
    spMonotonic.printCSV(x_[0] - 1.00001,
                         x_[n] + 1.00001,
                         1000);
    std::cout << "\n";
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    try {
        testAll();
        plot();
    }
    catch (const std::exception& e) {
        std::cout << "Caught OPM exception: " << e.what() << "\n";
    }

    return 0;
}
