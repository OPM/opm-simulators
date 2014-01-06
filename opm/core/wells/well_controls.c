/*===========================================================================
//
// File: newwells.c
//
// Created: 2012-02-03 11:28:40+0100
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include "config.h"

#define HAVE_WELLCONTROLS
#include <opm/core/well_controls.h>

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>




/* ---------------------------------------------------------------------- */
void
well_controls_destroy(struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    if (ctrl != NULL) {
        free             (ctrl->distr);
        free             (ctrl->target);
        free             (ctrl->type);
    }

    free(ctrl);
}


/* ---------------------------------------------------------------------- */
struct WellControls *
well_controls_create(void)
/* ---------------------------------------------------------------------- */
{
    struct WellControls *ctrl;

    ctrl = malloc(1 * sizeof *ctrl);

    if (ctrl != NULL) {
        /* Initialise empty control set */
        ctrl->num               = 0;
        ctrl->number_of_phases  = 0;
        ctrl->type              = NULL;
        ctrl->target            = NULL;
        ctrl->distr             = NULL;
        ctrl->current           = -1;
        ctrl->cpty              = 0;         
    }

    return ctrl;
}


/* ---------------------------------------------------------------------- */
static int
well_controls_reserve(int nctrl, struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    int   c, p, ok;
    void *type, *target, *distr;

    type   = realloc(ctrl->type  , nctrl * 1                      * sizeof *ctrl->type  );
    target = realloc(ctrl->target, nctrl * 1                      * sizeof *ctrl->target);
    distr  = realloc(ctrl->distr , nctrl * ctrl->number_of_phases * sizeof *ctrl->distr );

    ok = 0;
    if (type   != NULL) { ctrl->type   = type  ; ok++; }
    if (target != NULL) { ctrl->target = target; ok++; }
    if (distr  != NULL) { ctrl->distr  = distr ; ok++; }

    if (ok == 3) {
        for (c = ctrl->cpty; c < nctrl; c++) {
            ctrl->type  [c] =  BHP;
            ctrl->target[c] = -1.0;
        }

        for (p = ctrl->cpty * ctrl->number_of_phases; p < nctrl * ctrl->number_of_phases; ++p) {
            ctrl->distr[ p ] = 0.0;
        }

        ctrl->cpty = nctrl;
    }

    return ok == 3;
}


int well_controls_get_num(const struct WellControls *ctrl) {
  return ctrl->num;
}


int well_controls_get_cpty(const struct WellControls *ctrl) {
  return ctrl->cpty;
}


int well_controls_get_current( const struct WellControls * ctrl) {
    return ctrl->current;
}

void
well_controls_set_current( struct WellControls * ctrl, int current) {
    ctrl->current = current;
}

void 
well_controls_invert_current( struct WellControls * ctrl ) {
    ctrl->current = ~ctrl->current;
}


enum WellControlType 
well_controls_iget_type(const struct WellControls * ctrl, int control_index) {
    return ctrl->type[control_index];
}

void
well_controls_iset_type( struct WellControls * ctrls , int control_index , enum WellControlType type) {
    ctrls->type[control_index] = type;
}


double
well_controls_iget_target(const struct WellControls * ctrl, int control_index) {
    return ctrl->target[control_index];
}

void
well_controls_iset_target(struct WellControls * ctrl, int control_index , double target) {
    ctrl->target[control_index] = target;
}


const double *
well_controls_iget_distr(const struct WellControls * ctrl, int control_index) {
    int offset = control_index * ctrl->number_of_phases;
    return &ctrl->distr[offset];
}


void 
well_controls_iset_distr(const struct WellControls * ctrl, int control_index, const double * distr) {
    int offset = control_index * ctrl->number_of_phases;
    for (int p=0; p < ctrl->number_of_phases; p++)
        ctrl->distr[offset + p] = distr[p];
}


void 
well_controls_assert_number_of_phases(struct WellControls * ctrl , int number_of_phases) {
    if (ctrl->num == 0)
        ctrl->number_of_phases = number_of_phases;

    assert( ctrl->number_of_phases == number_of_phases );
}

void 
well_controls_clear(struct WellControls * ctrl) {
    ctrl->num = 0;
    ctrl->number_of_phases = 0;
}



int
well_controls_add_new(enum WellControlType type , double target , const double * distr , struct WellControls * ctrl) { 
    if (ctrl->num == ctrl->cpty) {
        int new_cpty = 2*ctrl->cpty;
        if (new_cpty == ctrl->num)
            new_cpty += 1;

        if (!well_controls_reserve( new_cpty , ctrl))
            return 0;
    }

    ctrl->type  [ctrl->num] = type  ;
    ctrl->target[ctrl->num] = target;
    
    if (distr != NULL) {
        int offset = ctrl->num * ctrl->number_of_phases;
        memcpy(&ctrl->distr[offset] , distr, ctrl->number_of_phases * sizeof * ctrl->distr);
    }
    
    ctrl->num += 1;
    return 1;
}


bool
well_controls_equal(const struct WellControls *ctrls1, const struct WellControls *ctrls2)
/* ---------------------------------------------------------------------- */
{
    bool are_equal = true;
    are_equal = (ctrls1->num == ctrls2->num);
    are_equal &= (ctrls1->number_of_phases == ctrls2->number_of_phases);
    if (!are_equal) {
        return are_equal;
    }

    are_equal &= (memcmp(ctrls1->type, ctrls2->type, ctrls1->num * sizeof *ctrls1->type ) == 0);
    are_equal &= (memcmp(ctrls1->target, ctrls2->target, ctrls1->num * sizeof *ctrls1->target ) == 0);
    are_equal &= (memcmp(ctrls1->distr, ctrls2->distr, ctrls1->num * ctrls1->number_of_phases * sizeof *ctrls1->distr ) == 0);
    are_equal &= (ctrls1->cpty == ctrls2->cpty);

    return are_equal;
}

