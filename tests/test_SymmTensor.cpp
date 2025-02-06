/*
  Copyright 2025 Equinor ASA

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

#include <opm/simulators/utils/SymmTensor.hpp>

#define BOOST_TEST_MODULE SymmTensorTest
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>

namespace {
    #if FLOW_INSTANTIATE_FLOAT
    using Types = boost::mpl::list<float,double>;
    #else
    using Types = boost::mpl::list<double>;
    #endif
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Basics, Scalar, Types)
{
    Opm::SymmTensor<Scalar> tensor;
    constexpr auto& indices = Opm::SymmTensor<Scalar>::unique_indices;
    for (const auto i : indices) {
        BOOST_CHECK_EQUAL(tensor[i], Scalar{0});
    }

    tensor = Scalar{1};
    for (const auto i : indices) {
        BOOST_CHECK_EQUAL(tensor[i], Scalar{1});
    }

    const Opm::SymmTensor<Scalar> tensor2{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::XX], 1.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::YY], 2.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::ZZ], 3.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::YZ], 4.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::YZ], tensor2[Opm::VoigtIndex::ZY]);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::XZ], 5.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::XZ], tensor2[Opm::VoigtIndex::ZX]);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::XY], 6.0);
    BOOST_CHECK_EQUAL(tensor2[Opm::VoigtIndex::XY], tensor2[Opm::VoigtIndex::YX]);

    const Opm::SymmTensor<Scalar> tensor2s{1.0, 2.0, 3.0};
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::XX], 1.0);
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::YY], 2.0);
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::ZZ], 3.0);
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::YZ], 0.0);
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::XZ], 0.0);
    BOOST_CHECK_EQUAL(tensor2s[Opm::VoigtIndex::XY], 0.0);

    const auto tensor3 = Scalar{2} * tensor + tensor2;
    auto tensor4 = tensor;
    tensor4 *= Scalar{2};
    tensor += tensor2;
    constexpr auto& uindices = Opm::SymmTensor<Scalar>::unique_indices;
    for (const auto i : uindices) {
        BOOST_CHECK_EQUAL(tensor[i], tensor2[i] + 1.0);
        BOOST_CHECK_EQUAL(tensor3[i], tensor2[i] + 2.0);
        BOOST_CHECK_EQUAL(tensor4[i], 2.0);
    }
    BOOST_CHECK_EQUAL(tensor.trace(), Scalar{2} + Scalar{3} + Scalar{4});
    BOOST_CHECK_EQUAL(tensor.traction({1.0, 0.0, 0.0}), 2.0);
    BOOST_CHECK_EQUAL(tensor.traction({0.0, 1.0, 0.0}), 3.0);
    BOOST_CHECK_EQUAL(tensor.traction({0.0, 0.0, 1.0}), 4.0);
    BOOST_CHECK_EQUAL(tensor.traction({1.0, 1.0, 0.0}), 2.0 + 3.0 + 2 * 5.0);
    BOOST_CHECK_EQUAL(tensor.traction({1.0, 0.0, 1.0}), 2.0 + 4.0 + 2 * 6.0);
    BOOST_CHECK_EQUAL(tensor.traction({0.0, 1.0, 1.0}), 3.0 + 4.0 + 2 * 7.0);
}
