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

#ifndef OPM_WELL_CONTROLS_H_INCLUDED
#define OPM_WELL_CONTROLS_H_INCLUDED

#include <stdbool.h>



#ifdef __cplusplus
extern "C" {
#endif

enum WellControlType  {
    BHP,              /**< Well constrained by BHP target */
    RESERVOIR_RATE,   /**< Well constrained by reservoir volume flow rate */
    SURFACE_RATE      /**< Well constrained by surface volume flow rate */
};

struct WellControls;

bool 
well_controls_equal(const struct WellControls *ctrls1, const struct WellControls *ctrls2);

//int 
//well_controls_reserve(int nctrl, int nphases, struct WellControls *ctrl);

struct WellControls * 
well_controls_create(void);

void
well_controls_destroy(struct WellControls *ctrl);


int 
well_controls_get_num(const struct WellControls *ctrl);

int 
well_controls_get_current( const struct WellControls * ctrl);

void
well_controls_set_current( struct WellControls * ctrl, int current);

void
well_controls_invert_current( struct WellControls * ctrl );

int
well_controls_add_new(enum WellControlType type , double target , const double * distr , struct WellControls * ctrl);

enum WellControlType 
well_controls_iget_type(const struct WellControls * ctrl, int control_index);

enum WellControlType 
well_controls_get_current_type(const struct WellControls * ctrl);

void
well_controls_iset_type( struct WellControls * ctrls , int control_index , enum WellControlType type);

void
well_controls_iset_target(struct WellControls * ctrl, int control_index , double target);

double
well_controls_iget_target(const struct WellControls * ctrl, int control_index);

double
well_controls_get_current_target(const struct WellControls * ctrl);

const double *
well_controls_iget_distr(const struct WellControls * ctrl, int control_index);

void 
well_controls_iset_distr(const struct WellControls * ctrl, int control_index, const double * distr);

const double *
well_controls_get_current_distr(const struct WellControls * ctrl);

void 
well_controls_assert_number_of_phases(struct WellControls * ctrl , int number_of_phases);

void 
well_controls_clear(struct WellControls * ctrl);


#ifdef __cplusplus
}
#endif

#endif /* OPM_WELL_CONTROLS_H_INCLUDED */
