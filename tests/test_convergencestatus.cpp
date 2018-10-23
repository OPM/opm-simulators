/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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
#define BOOST_TEST_MODULE ConvergenceStatusTest
#include <boost/test/unit_test.hpp>

#include <opm/simulators/timestepping/ConvergenceStatus.hpp>

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    Opm::ConvergenceStatus s;
    BOOST_CHECK(s.converged());
    BOOST_CHECK(!s.reservoirFailed());
    BOOST_CHECK(!s.wellFailed());
}

BOOST_AUTO_TEST_CASE(Failures)
{
    using CS = Opm::ConvergenceStatus;
    Opm::ConvergenceStatus s1;
    s1.setReservoirFailed({CS::ReservoirFailure::Type::Cnv, 2, 100});
    {
        BOOST_CHECK(!s1.converged());
        BOOST_CHECK(s1.reservoirFailed());
        BOOST_CHECK(!s1.wellFailed());
        BOOST_REQUIRE(s1.reservoirFailures().size() == 1);
        const auto f = s1.reservoirFailures()[0];
        BOOST_CHECK(f.type == CS::ReservoirFailure::Type::Cnv);
        BOOST_CHECK(f.phase == 2);
        BOOST_CHECK(f.cell_index == 100);
        BOOST_CHECK(s1.wellFailures().empty());
    }

    Opm::ConvergenceStatus s2;
    s2.setWellFailed({CS::WellFailure::Type::Ctrl, -1, "PRODUCER-123"});
    s2.setWellFailed({CS::WellFailure::Type::Mb, 2, "INJECTOR-XYZ"});
    {
        BOOST_CHECK(!s2.converged());
        BOOST_CHECK(!s2.reservoirFailed());
        BOOST_CHECK(s2.wellFailed());
        BOOST_CHECK(s2.reservoirFailures().empty());
        BOOST_REQUIRE(s2.wellFailures().size() == 2);
        const auto f0 = s2.wellFailures()[0];
        BOOST_CHECK(f0.type == CS::WellFailure::Type::Ctrl);
        BOOST_CHECK(f0.phase == -1);
        BOOST_CHECK(f0.well_name == "PRODUCER-123");
        const auto f1 = s2.wellFailures()[1];
        BOOST_CHECK(f1.type == CS::WellFailure::Type::Mb);
        BOOST_CHECK(f1.phase == 2);
        BOOST_CHECK(f1.well_name == "INJECTOR-XYZ");
    }

    s1 += s2;
    {
        BOOST_CHECK(!s1.converged());
        BOOST_CHECK(s1.reservoirFailed());
        BOOST_CHECK(s1.wellFailed());
        BOOST_REQUIRE(s1.reservoirFailures().size() == 1);
        const auto f = s1.reservoirFailures()[0];
        BOOST_CHECK(f.type == CS::ReservoirFailure::Type::Cnv);
        BOOST_CHECK(f.phase == 2);
        BOOST_CHECK(f.cell_index == 100);
        BOOST_REQUIRE(s1.wellFailures().size() == 2);
        const auto f0 = s1.wellFailures()[0];
        BOOST_CHECK(f0.type == CS::WellFailure::Type::Ctrl);
        BOOST_CHECK(f0.phase == -1);
        BOOST_CHECK(f0.well_name == "PRODUCER-123");
        const auto f1 = s1.wellFailures()[1];
        BOOST_CHECK(f1.type == CS::WellFailure::Type::Mb);
        BOOST_CHECK(f1.phase == 2);
        BOOST_CHECK(f1.well_name == "INJECTOR-XYZ");
    }

    s1.clear();
    {
        BOOST_CHECK(s1.converged());
        BOOST_CHECK(!s1.reservoirFailed());
        BOOST_CHECK(!s1.wellFailed());
    }

    s1 += s2;
    {
        BOOST_CHECK(!s1.converged());
        BOOST_CHECK(!s1.reservoirFailed());
        BOOST_CHECK(s1.wellFailed());
        BOOST_CHECK(s1.reservoirFailures().empty());
        BOOST_REQUIRE(s1.wellFailures().size() == 2);
        const auto f0 = s1.wellFailures()[0];
        BOOST_CHECK(f0.type == CS::WellFailure::Type::Ctrl);
        BOOST_CHECK(f0.phase == -1);
        BOOST_CHECK(f0.well_name == "PRODUCER-123");
        const auto f1 = s1.wellFailures()[1];
        BOOST_CHECK(f1.type == CS::WellFailure::Type::Mb);
        BOOST_CHECK(f1.phase == 2);
        BOOST_CHECK(f1.well_name == "INJECTOR-XYZ");
    }
}

