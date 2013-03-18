/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#define BOOST_TEST_MODULE WellsModuleTest
#include <boost/test/unit_test.hpp>

#include <opm/core/wells.h>

#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

BOOST_AUTO_TEST_CASE(Construction)
{
    const int nphases = 2;
    const int nwells  = 2;
    const int nperfs  = 2;

    boost::shared_ptr<Wells> W(create_wells(nphases, nwells, nperfs),
                                destroy_wells);

    if (W) {
        int          cells[] = { 0, 9 };
        double       WI      = 1.0;
        const double ifrac[] = { 1.0, 0.0 };

        const bool ok0 = add_well(INJECTOR, 0.0, 1, &ifrac[0], &cells[0],
                                  &WI, "INJECTOR", W.get());

        const double pfrac[] = { 0.0, 0.0 };
        const bool ok1 = add_well(PRODUCER, 0.0, 1, &pfrac[0], &cells[1],
                                  &WI, "PRODUCER", W.get());

        if (ok0 && ok1) {
            BOOST_CHECK_EQUAL(W->number_of_phases, nphases);
            BOOST_CHECK_EQUAL(W->number_of_wells , nwells );

            BOOST_CHECK_EQUAL(W->well_connpos[0], 0);
            BOOST_CHECK_EQUAL(W->well_connpos[1], 1);
            BOOST_CHECK_EQUAL(W->well_connpos[W->number_of_wells], nperfs);

            BOOST_CHECK_EQUAL(W->well_cells[W->well_connpos[0]], cells[0]);
            BOOST_CHECK_EQUAL(W->well_cells[W->well_connpos[1]], cells[1]);

            BOOST_CHECK_EQUAL(W->WI[W->well_connpos[0]], WI);
            BOOST_CHECK_EQUAL(W->WI[W->well_connpos[1]], WI);

            using std::string;
            BOOST_CHECK_EQUAL(string(W->name[0]), string("INJECTOR"));
            BOOST_CHECK_EQUAL(string(W->name[1]), string("PRODUCER"));
        }
    }
}


BOOST_AUTO_TEST_CASE(Controls)
{
    const int nphases = 2;
    const int nwells  = 1;
    const int nperfs  = 2;

    boost::shared_ptr<Wells> W(create_wells(nphases, nwells, nperfs),
                               destroy_wells);

    if (W) {
        int          cells[] = { 0  , 9   };
        double       WI   [] = { 1.0, 1.0 };
        const double ifrac[] = { 1.0, 0.0 };

        const bool ok = add_well(INJECTOR, 0.0, nperfs, &ifrac[0], &cells[0],
                                 &WI[0], "INJECTOR", W.get());

        if (ok) {
            const double distr[] = { 1.0, 0.0 };
            const bool   ok1     = append_well_controls(BHP, 1, &distr[0],
                                                        0, W.get());
            const bool   ok2     = append_well_controls(SURFACE_RATE, 1,
                                                        &distr[0], 0, W.get());

            if (ok1 && ok2) {
                WellControls* ctrls = W->ctrls[0];

                BOOST_CHECK_EQUAL(ctrls->num    ,  2);
                BOOST_CHECK_EQUAL(ctrls->current, -1);

                set_current_control(0, 0, W.get());
                BOOST_CHECK_EQUAL(ctrls->current, 0);

                set_current_control(0, 1, W.get());
                BOOST_CHECK_EQUAL(ctrls->current, 1);

                BOOST_CHECK_EQUAL(ctrls->type[0], BHP);
                BOOST_CHECK_EQUAL(ctrls->type[1], SURFACE_RATE);

                BOOST_CHECK_EQUAL(ctrls->target[0], 1.0);
                BOOST_CHECK_EQUAL(ctrls->target[1], 1.0);
            }
        }
    }
}


BOOST_AUTO_TEST_CASE(Copy)
{
    const int nphases = 2;
    const int nwells  = 2;
    const int nperfs  = 2;

    boost::shared_ptr<Wells> W1(create_wells(nphases, nwells, nperfs),
                                destroy_wells);
    boost::shared_ptr<Wells> W2;

    if (W1) {
        int          cells[] = { 0, 9 };
        const double WI      = 1.0;
        const double ifrac[] = { 1.0, 0.0 };

        const bool ok0 = add_well(INJECTOR, 0.0, 1, &ifrac[0], &cells[0],
                                  &WI, "INJECTOR", W1.get());

        const double pfrac[] = { 0.0, 0.0 };
        const bool ok1 = add_well(PRODUCER, 0.0, 1, &pfrac[0], &cells[1],
                                  &WI, "PRODUCER", W1.get());

        bool ok = ok0 && ok1;
        for (int w = 0; ok && (w < W1->number_of_wells); ++w) {
            const double distr[] = { 1.0, 0.0 };
            const bool   okc1     = append_well_controls(BHP, 1, &distr[0],
                                                         w, W1.get());
            const bool   okc2     = append_well_controls(SURFACE_RATE, 1,
                                                         &distr[0], w,
                                                         W1.get());

            ok = okc1 && okc2;
        }

        if (ok) {
            W2.reset(clone_wells(W1.get()), destroy_wells);
        }
    }

    if (W2) {
        BOOST_CHECK_EQUAL(W2->number_of_phases, W1->number_of_phases);
        BOOST_CHECK_EQUAL(W2->number_of_wells , W1->number_of_wells );
        BOOST_CHECK_EQUAL(W2->well_connpos[0] , W1->well_connpos[0] );

        for (int w = 0; w < W1->number_of_wells; ++w) {
            using std::string;
            BOOST_CHECK_EQUAL(string(W2->name[w]), string(W1->name[w]));
            BOOST_CHECK_EQUAL(       W2->type[w] ,        W1->type[w] );

            BOOST_CHECK_EQUAL(W2->well_connpos[w + 1],
                              W1->well_connpos[w + 1]);

            for (int j = W1->well_connpos[w];
                j < W1->well_connpos[w + 1]; ++j) {
                BOOST_CHECK_EQUAL(W2->well_cells[j], W1->well_cells[j]);
                BOOST_CHECK_EQUAL(W2->WI        [j], W1->WI        [j]);
            }

            BOOST_CHECK(W1->ctrls[w] != 0);
            BOOST_CHECK(W2->ctrls[w] != 0);

            WellControls* c1 = W1->ctrls[w];
            WellControls* c2 = W2->ctrls[w];

            BOOST_CHECK_EQUAL(c2->num    , c1->num    );
            BOOST_CHECK_EQUAL(c2->current, c1->current);

            for (int c = 0; c < c1->num; ++c) {
                BOOST_CHECK_EQUAL(c2->type  [c], c1->type  [c]);
                BOOST_CHECK_EQUAL(c2->target[c], c1->target[c]);

                for (int p = 0; p < W1->number_of_phases; ++p) {
                    BOOST_CHECK_EQUAL(c2->distr[c*W1->number_of_phases + p],
                                      c1->distr[c*W1->number_of_phases + p]);
                }
            }
        }
    }
}
