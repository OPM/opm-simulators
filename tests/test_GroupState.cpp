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

#include <opm/simulators/wells/GroupState.hpp>
#include <opm/json/JsonObject.hpp>


#define BOOST_TEST_MODULE GroupStateTest
#include <boost/test/unit_test.hpp>

using namespace Opm;

class TestCommunicator {
public:
    void sum(const double *, std::size_t) const {}
};



BOOST_AUTO_TEST_CASE(GroupStateCreate) {
    std::size_t num_phases{3};
    GroupState gs(num_phases);

    BOOST_CHECK(!gs.has_production_rates("AGROUP"));
    BOOST_CHECK_THROW( gs.update_production_rates("AGROUP", {0}), std::exception);

    std::vector<double> rates{0,1,2};
    gs.update_production_rates("AGROUP", rates);
    BOOST_CHECK(gs.has_production_rates("AGROUP"));

    BOOST_CHECK_THROW( gs.production_rates("NO_SUCH_GROUP"), std::exception );
    auto r2 = gs.production_rates("AGROUP");
    BOOST_CHECK( r2 == rates );

    gs.update_injection_rein_rates("CGROUP", rates);


    BOOST_CHECK(!gs.has_production_control("NO_SUCH_GROUP"));
    BOOST_CHECK(!gs.has_production_control("AGROUP"));

    gs.production_control("AGROUP", Group::ProductionCMode::GRAT);
    BOOST_CHECK(gs.has_production_control("AGROUP"));
    BOOST_CHECK(gs.production_control("AGROUP") == Group::ProductionCMode::GRAT);
    BOOST_CHECK_THROW(gs.production_control("BGROUP"), std::exception);


    auto gs2 = gs;
    BOOST_CHECK(gs2 == gs);
    TestCommunicator comm;
    gs.communicate_rates(comm);
    BOOST_CHECK(gs2 == gs);
}


