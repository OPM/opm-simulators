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

#define BOOST_TEST_MODULE AutoDiffBlockTest

#include <opm/autodiff/AutoDiffBlock.hpp>

#include <boost/test/unit_test.hpp>

#include <Eigen/Eigen>
#include <Eigen/Sparse>



BOOST_AUTO_TEST_CASE(ScalarMultiplication)
{
    typedef AutoDiff::ForwardBlock<double> ADB;
    std::vector<int> blocksizes = { 3, 1, 2 };

    ADB::V vx(3);
    vx << 1.0, 2.0, 10.0;
    
    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADB x = ADB::variable(FirstVar, vx, blocksizes);

    ADB::V const_vector(3);
    const_vector << 3.14, 3.14, 3.14;
    // std::cout << "const_vector:\n" << const_vector << std::endl;
    
    const ADB x2 = x * const_vector;
    const ADB x3 = const_vector * x;
    BOOST_CHECK_EQUAL( x2.value().cwiseNotEqual( x3.value() ).count(), 0 );

    // The new operator:

    const ADB y2 = x * 3.14;
    BOOST_CHECK_EQUAL( x2.value().cwiseNotEqual( y2.value() ).count(), 0 );

    const ADB y3 = 3.14 * x;
    BOOST_CHECK_EQUAL( x3.value().cwiseNotEqual( y3.value() ).count(), 0 );
}



