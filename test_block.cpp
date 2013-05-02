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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE AutoDiffBlockTest

#include "AutoDiffBlock.hpp"

#include <boost/test/unit_test.hpp>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <iostream>

BOOST_AUTO_TEST_CASE(ConstantInitialisation)
{
    typedef AutoDiff::ForwardBlock<double> ADV;

    std::vector<int> blocksizes = { 3, 1, 2 };

    ADV::V v(3);
    v << 0.2, 1.2, 13.4;

    ADV a = ADV::constant(v, blocksizes);
    BOOST_REQUIRE(a.value().matrix() == v.matrix());

    const std::vector<ADV::M>& J = a.derivative();
    for (std::vector<ADV::M>::const_iterator
             b = J.begin(), e = J.end(); b != e; ++b) {
        BOOST_REQUIRE(b->nonZeros() == 0);
    }
}

BOOST_AUTO_TEST_CASE(VariableInitialisation)
{
    typedef AutoDiff::ForwardBlock<double> ADV;

    std::vector<int> blocksizes = { 3, 1, 2 };

    ADV::V v(3);
    v << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADV x = ADV::variable(FirstVar, v, blocksizes);

    BOOST_REQUIRE(x.value().matrix() == v.matrix());

    const std::vector<ADV::M>& J = x.derivative();
    BOOST_REQUIRE(J[0].nonZeros() == v.size());

    const Eigen::Diagonal<const ADV::M, 0>& d = J[0].diagonal();
    BOOST_REQUIRE((d.array() == 1.0).all());

    for (std::vector<ADV::M>::const_iterator
             b = J.begin() + 1, e = J.end(); b != e; ++b) {
        BOOST_REQUIRE(b->nonZeros() == 0);
    }
}

BOOST_AUTO_TEST_CASE(FunctionInitialisation)
{
    typedef AutoDiff::ForwardBlock<double> ADV;

    std::vector<int>            blocksizes = { 3, 1, 2 };
    std::vector<int>::size_type num_blocks = blocksizes.size();

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADV::V v(3);
    v << 1.0, 2.2, 3.4;

    std::vector<ADV::M> jacs(num_blocks);
    for (std::vector<int>::size_type j = 0; j < num_blocks; ++j) {
        jacs[j] = ADV::M(blocksizes[FirstVar], blocksizes[j]);
        jacs[j].insert(0,0) = -1.0;
    }

    ADV f = ADV::function(v, jacs);

    BOOST_REQUIRE(f.value().matrix() == v.matrix());

    const std::vector<ADV::M>& J = f.derivative();
    for (std::vector<ADV::M>::const_iterator
             bf = J.begin(), ef = J.end(), bj = jacs.begin();
         bf != ef; ++bf, ++bj) {
        BOOST_REQUIRE(bf->nonZeros() == bj->nonZeros());

        BOOST_REQUIRE(bf->outerSize() == bj->outerSize());
        BOOST_REQUIRE(bf->innerSize() == bj->innerSize());

        for (ADV::M::Index k = 0; k < bf->outerSize(); ++k) {
            for (ADV::M::InnerIterator
                     ileft(*bf, k), iright(*bj, k);
                 ileft && iright; ++ileft, ++iright) {
                BOOST_REQUIRE(ileft.row()   == iright.row()  );
                BOOST_REQUIRE(ileft.col()   == iright.col()  );
                BOOST_REQUIRE(ileft.index() == iright.index());
                BOOST_REQUIRE(ileft.value() == iright.value());
            }
        }
    }
}

#if 0
#include <iostream>

int main()
{
    typedef AutoDiff::ForwardBlock<double> ADV;
    std::vector<int> blocksizes = { 3, 1, 2 };
    int num_blocks = blocksizes.size();
    ADV::V v1(3);
    v1 << 0.2, 1.2, 13.4;
    ADV::V v2(3);
    v2 << 1.0, 2.2, 3.4;
    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };
    ADV a = ADV::constant(v1, blocksizes);
    ADV x = ADV::variable(FirstVar, v2, blocksizes);
    std::vector<ADV::M> jacs(num_blocks);
    for (int i = 0; i < num_blocks; ++i) {
        jacs[i] = ADV::M(blocksizes[FirstVar], blocksizes[i]);
        jacs[i].insert(0,0) = -1.0;
    }
    ADV f = ADV::function(v2, jacs);

    ADV xpx = x + x;
    std::cout << xpx;
    ADV xpxpa = x + x + a;
    std::cout << xpxpa;


    std::cout << xpxpa - xpx;

    ADV sqx = x * x;

    std::cout << sqx;

    ADV sqxdx = sqx / x;

    std::cout << sqxdx;

    ADV::M m(2,3);
    m.insert(0,0) = 4;
    m.insert(0,1) = 3;
    m.insert(1,1) = 1;
    std::cout << m*sqx;
}
#endif
