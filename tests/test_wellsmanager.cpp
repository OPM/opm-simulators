/*
  Copyright 2013 Statoil ASA.

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

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE WellsManagerTests
#include <boost/test/unit_test.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>

BOOST_AUTO_TEST_CASE(Constructor_Old_Works) {
    Opm::EclipseGridParser oldParser("wells_manager_data.data");
    Opm::GridManager gridManager(oldParser);
    Opm::WellsManager wellsManager(oldParser, *gridManager.c_grid(), NULL);
    const Wells* oldWells = wellsManager.c_wells();
    BOOST_CHECK_EQUAL(oldWells->number_of_wells, 2);
}
