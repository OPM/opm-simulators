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

#ifndef OPM_NEWWELLS_H_INCLUDED
#define OPM_NEWWELLS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

enum well_type         { INJECTOR, PRODUCER };
enum control_type      { BHP     , RATE, GRUP     };
enum surface_component { WATER = 0, OIL = 1, GAS = 2 };


/* Control for single well */
struct WellControls {
    int                  num;
    int                  cpty;
    enum control_type   *type;
    double              *target;
    int                  current;
};

struct Wells {
    int                  number_of_wells;
    int                  well_cpty;
    int                  perf_cpty;

    enum well_type      *type;
    double              *depth_ref;
    double              *zfrac;        /* Surface volume fraction
                                        * (3*number_of_wells) */

    int                 *well_connpos;
    int                 *well_cells;
    double              *WI;           /* Well index */

    struct WellControls **ctrls;       /* One struct for each well */
};

struct CompletionData {
    double *gpot;               /* Gravity potential */
    double *A;                  /* RB^{-1} */
    double *phasemob;           /* Phase mobility */
};

struct Wells *
wells_create(int nwells, int nperf);

int
wells_add(enum well_type type     ,
          double         depth_ref,
          int            nperf    ,
          const double  *zfrac    , /* Injection fraction or NULL */
          const int     *cells    ,
          const double  *WI       , /* Well index per perf (or NULL) */
          struct Wells  *W        );

void
wells_destroy(struct Wells *W);

int
well_controls_append(enum control_type    type  ,
                     double               target,
                     struct WellControls *ctrl  );

void
well_controls_clear(struct WellControls *ctrl);

#ifdef __cplusplus
}
#endif

#endif /* OPM_NEWWELLS_H_INCLUDED */
