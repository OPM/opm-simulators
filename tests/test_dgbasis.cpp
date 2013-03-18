/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE DGBasisTest
#include <boost/test/unit_test.hpp>

#include <opm/core/transport/reorder/DGBasis.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>
#include <cmath>

using namespace Opm;


namespace
{


    bool aequal(double a, double b)
    {
        const double eps = 1e-15;
        return std::fabs(a - b) < eps;
    }

} // anonymous namespace


namespace cart2d
{

    static void test()
    {
        // Set up 2d 1-cell cartesian case.
        GridManager g(1, 1);
        const UnstructuredGrid& grid = *g.c_grid();

        // Test DGBasisBoundedTotalDegree, degree 0.
        {
            DGBasisBoundedTotalDegree b(grid, 0);
            BOOST_CHECK_EQUAL(b.numBasisFunc(), 1);
            std::vector<double> bx(b.numBasisFunc(), 0.0);
            b.eval(0, grid.cell_centroids, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            double x[2] = { 0.123, 0.456 };
            b.eval(0, x, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            std::vector<double> c(b.numBasisFunc(), 0.0);
            b.addConstant(0.789, &c[0]);
            BOOST_CHECK(aequal(c[0], 0.789));
            b.multiplyGradient(1.234, &c[0]);
            BOOST_CHECK(aequal(c[0], 0.789));
        }
        // Test DGBasisBoundedTotalDegree, degree 1.
        {
            DGBasisBoundedTotalDegree b(grid, 1);
            BOOST_CHECK_EQUAL(b.numBasisFunc(), 3);
            std::vector<double> bx(b.numBasisFunc(), 0.0);
            b.eval(0, grid.cell_centroids, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            BOOST_CHECK(aequal(bx[1], 0.0));
            BOOST_CHECK(aequal(bx[2], 0.0));
            double x[2] = { 0.123, 0.456 };
            b.eval(0, x, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            BOOST_CHECK(aequal(bx[1], 0.123 - 0.5));
            BOOST_CHECK(aequal(bx[2], 0.456 - 0.5));
            std::vector<double> c(b.numBasisFunc(), 0.0);
            c[0] = 1.0; c[1] = 2.0; c[2] = 3.0;
            b.addConstant(0.789, &c[0]);
            BOOST_CHECK(aequal(c[0], 1.789));
            BOOST_CHECK(aequal(c[1], 2.0));
            BOOST_CHECK(aequal(c[2], 3.0));
            const double fx = c[0]*bx[0] + c[1]*bx[1] + c[2]*bx[2];
            b.multiplyGradient(1.234, &c[0]);
            BOOST_CHECK(aequal(c[0], 1.789));
            BOOST_CHECK(aequal(c[1], 2.0*1.234));
            BOOST_CHECK(aequal(c[2], 3.0*1.234));
            const double fx2 = c[0]*bx[0] + c[1]*bx[1] + c[2]*bx[2];
            BOOST_CHECK(aequal(fx2 - c[0], 1.234*(fx - c[0])));
        }
        // Test DGBasisMultilin, degree 0.
        {
            DGBasisMultilin b(grid, 0);
            BOOST_CHECK_EQUAL(b.numBasisFunc(), 1);
            std::vector<double> bx(b.numBasisFunc(), 0.0);
            b.eval(0, grid.cell_centroids, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            double x[2] = { 0.123, 0.456 };
            b.eval(0, x, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 1.0));
            std::vector<double> c(b.numBasisFunc(), 0.0);
            b.addConstant(0.789, &c[0]);
            BOOST_CHECK(aequal(c[0], 0.789));
            b.multiplyGradient(1.234, &c[0]);
            BOOST_CHECK(aequal(c[0], 0.789));
        }
        // Test DGBasisMultilin, degree 1.
        {
            DGBasisMultilin b(grid, 1);
            BOOST_CHECK_EQUAL(b.numBasisFunc(), 4);
            std::vector<double> bx(b.numBasisFunc(), 0.0);
            b.eval(0, grid.cell_centroids, &bx[0]);
            BOOST_CHECK(aequal(bx[0], 0.25));
            BOOST_CHECK(aequal(bx[1], 0.25));
            BOOST_CHECK(aequal(bx[2], 0.25));
            BOOST_CHECK(aequal(bx[3], 0.25));
            double x[2] = { 0.123, 0.456 };
            b.eval(0, x, &bx[0]);
            const double xm[2] = { 1.0 - x[0], x[0] };
            const double ym[2] = { 1.0 - x[1], x[1] };
            BOOST_CHECK(aequal(bx[0], xm[0]*ym[0]));
            BOOST_CHECK(aequal(bx[1], xm[0]*ym[1]));
            BOOST_CHECK(aequal(bx[2], xm[1]*ym[0]));
            BOOST_CHECK(aequal(bx[3], xm[1]*ym[1]));
            std::vector<double> c(b.numBasisFunc(), 0.0);
            c[0] = -1.567; c[1] = 1.42; c[2] = 0.59; c[3] = 3.225;
            std::vector<double> corig = c;
            b.addConstant(0.789, &c[0]);
            BOOST_CHECK(aequal(c[0], corig[0] + 0.25*0.789));
            BOOST_CHECK(aequal(c[1], corig[1] + 0.25*0.789));
            BOOST_CHECK(aequal(c[2], corig[2] + 0.25*0.789));
            BOOST_CHECK(aequal(c[3], corig[3] + 0.25*0.789));
            const double fx = c[0]*bx[0] + c[1]*bx[1] + c[2]*bx[2] + c[3]*bx[3];
            const double fc = 0.25*(c[0] + c[1] + c[2] + c[3]);
            b.multiplyGradient(1.234, &c[0]);
            const double fx2 = c[0]*bx[0] + c[1]*bx[1] + c[2]*bx[2] + c[3]*bx[3];
            BOOST_CHECK(aequal(fx2 - fc, 1.234*(fx - fc)));
        }

    }
} // namespace cart2d


BOOST_AUTO_TEST_CASE(test_dgbasis)
{
    cart2d::test();
}
