/*
  Copyright 2021 Equinor.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <stdexcept>
#include <opm/simulators/wells/ALQState.hpp>


#define BOOST_TEST_MODULE GroupStateTest
#include <boost/test/unit_test.hpp>

using namespace Opm;



BOOST_AUTO_TEST_CASE(ALQStateCreate) {
    ALQState alq_state;


    alq_state.update_default("W1", 100);
    alq_state.update_default("W2", 200);
    alq_state.set("W1", 1);

    BOOST_CHECK_THROW( alq_state.get("W3"), std::exception);
    BOOST_CHECK_EQUAL( alq_state.get("W2"), 200);
    BOOST_CHECK_EQUAL( alq_state.get("W1"), 1);

    alq_state.update_default("W1", 100);
    BOOST_CHECK_EQUAL( alq_state.get("W1"), 1);

    alq_state.update_default("W1", 0);
    BOOST_CHECK_EQUAL( alq_state.get("W1"), 0);

    BOOST_CHECK(!alq_state.oscillation("WX"));

    alq_state.update_count("WX", true);
    BOOST_CHECK(!alq_state.oscillation("WX"));

    alq_state.update_count("WX", false);
    BOOST_CHECK(alq_state.oscillation("WX"));

    alq_state.reset_count();
    BOOST_CHECK_EQUAL(alq_state.get_increment_count("NO_SUCH_WELL"), 0);
    BOOST_CHECK_EQUAL(alq_state.get_increment_count("WX"), 0);
    alq_state.update_count("WX", true);
    BOOST_CHECK_EQUAL(alq_state.get_increment_count("WX"), 1);

    BOOST_CHECK_EQUAL(alq_state.get_decrement_count("WX"), 0);
    alq_state.update_count("WX", false);
    BOOST_CHECK_EQUAL(alq_state.get_decrement_count("WX"), 1);

    alq_state.set("W1", 1);
    alq_state.set("W2", 2);
    BOOST_CHECK_EQUAL(alq_state.pack_size(), 2);
    std::vector<double> data(2);

    BOOST_CHECK_EQUAL( alq_state.get("W1"), 1);
    BOOST_CHECK_EQUAL( alq_state.get("W2"), 2);
    BOOST_CHECK_EQUAL(alq_state.pack_data(data.data()), 2);
    alq_state.set("W1", 99);
    alq_state.set("W2", 101);

    BOOST_CHECK_EQUAL(alq_state.unpack_data(data.data()), 2);
    BOOST_CHECK_EQUAL( alq_state.get("W1"), 1);
    BOOST_CHECK_EQUAL( alq_state.get("W2"), 2);
}
