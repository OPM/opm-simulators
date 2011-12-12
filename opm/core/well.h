/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELL_HEADER_INCLUDED
#define OPM_WELL_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif


enum well_type    { INJECTOR, PRODUCER };
enum well_control { BHP     , RATE     };
enum surface_component { WATER = 0, OIL = 1, GAS = 2 };

struct WellCompletions {
    int  number_of_wells;
    int *well_connpos;
    int *well_cells;
};

struct WellControls {
    enum well_type    *type;
    enum well_control *ctrl;
    double            *target;
    double            *zfrac;   /* Surface volume fraction (3*nwells) */
};

struct completion_data {
    double *WI;                 /* Productivity index */
    double *gpot;               /* Gravity potential */
    double *A;                  /* RB^{-1} */
    double *phasemob;           /* Phase mobility */
};

typedef struct WellCompletions well_t;
typedef struct WellControls    well_control_t;

int
allocate_cell_wells(int nc, well_t *W, int **cwpos, int **cwells);

void
deallocate_cell_wells(int *cvpos, int *cwells);

void
derive_cell_wells(int nc, well_t *W, int *cwpos, int *cwells);


#ifdef __cplusplus
}
#endif

#endif /* OPM_WELL_HEADER_INCLUDED */
