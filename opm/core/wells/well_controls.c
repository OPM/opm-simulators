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

#include <opm/core/well_controls.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Controls for a single well.
 * Each control specifies a well rate or bottom-hole pressure. Only
 * one control can be active at a time, indicated by current. The
 * meaning of each control's target value depends on the control type:
 *
 *  - BHP            -> target pressure in Pascal.
 *  - RESERVOIR_RATE -> target reservoir volume rate in cubic(meter)/second
 *  - SURFACE_RATE   -> target surface volume rate in cubic(meter)/second
 *
 * The sign convention for RATE targets is as follows:
 *
 *  - (+) Fluid flowing into reservoir, i.e. injecting.
 *  - (-) Fluid flowing out of reservoir, i.e. producing.
 *
 * For *_RATE controls, the distribution of phases used for the control
 * is also needed. For example, a total rate control should have 1.0
 * for each phase, whereas a control on oil rate should have 1.0 for
 * the oil phase and 0.0 for the rest. For BHP controls, this is unused.
 * The active control acts as an equality constraint, whereas the
 * non-active controls should be interpreted as inequality
 * constraints (upper or lower bounds).  For instance, a PRODUCER's
 * BHP constraint defines a minimum acceptable bottom-hole pressure
 * value for the well.
 */

struct WellControls
{
    /**
     * Number of controls.
     */
    int num;

    int number_of_phases;

    /**
     * Array of control types.
     */
    enum WellControlType *type;

    /**
     * Array of control targets.
     */
    double *target;

    /**
     * Array of rate control distributions,
     * <CODE>number_of_phases</CODE> numbers for each control
     */
    double *distr;

    /**
     * Index of current active control.
     */
    int current;

    bool well_is_open;

    /* 
       The capacity allocated.
    */
    int cpty;
};


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
        /* Initialise empty control set; the well is created open. */
        ctrl->num               = 0;
        ctrl->number_of_phases  = 0;
        ctrl->type              = NULL;
        ctrl->target            = NULL;
        ctrl->distr             = NULL;
        ctrl->current           = -1;
        ctrl->cpty              = 0;         
        ctrl->well_is_open      = true;  
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


/* ---------------------------------------------------------------------- */
struct WellControls *
well_controls_clone(const struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    int                   ok, i, n;
    double                target;
    const double         *distr;
    struct WellControls  *new;
    enum WellControlType  type;

    new = well_controls_create();

    if (new != NULL) {
        /* Assign appropriate number of phases */
        well_controls_assert_number_of_phases(new, ctrl->number_of_phases);

        n  = well_controls_get_num(ctrl);
        ok = well_controls_reserve(n, new);

        if (! ok) {
            well_controls_destroy(new);
            new = NULL;
        }
        else {
            for (i = 0; ok && (i < n); i++) {
                type   = well_controls_iget_type  (ctrl, i);
                distr  = well_controls_iget_distr (ctrl, i);
                target = well_controls_iget_target(ctrl, i);

                ok = well_controls_add_new(type, target, distr, new);
            }

            if (i < n) {
                assert (!ok);
                well_controls_destroy(new);

                new = NULL;
            }
            else {
                i = well_controls_get_current(ctrl);
                well_controls_set_current(new, i);

                if (well_controls_well_is_open(ctrl)) {
                    well_controls_open_well(new);
                }
                else {
                    well_controls_shut_well(new);
                }
            }
        }
    }

    assert (well_controls_equal(ctrl, new, true));

    return new;
}


int well_controls_get_num(const struct WellControls *ctrl) {
  return ctrl->num;
}


int well_controls_get_current( const struct WellControls * ctrl) {
    return ctrl->current;
}

void
well_controls_set_current( struct WellControls * ctrl, int current) {
    ctrl->current = current;
}

bool well_controls_well_is_shut(const struct WellControls * ctrl) {
    return !ctrl->well_is_open;
}

bool well_controls_well_is_open(const struct WellControls * ctrl) {
    return ctrl->well_is_open;
}

void well_controls_open_well( struct WellControls * ctrl) {
    ctrl->well_is_open = true;
}

void well_controls_shut_well( struct WellControls * ctrl) {
    ctrl->well_is_open = false;
}



enum WellControlType 
well_controls_iget_type(const struct WellControls * ctrl, int control_index) {
    return ctrl->type[control_index];
}


enum WellControlType 
well_controls_get_current_type(const struct WellControls * ctrl) {
    return well_controls_iget_type( ctrl , ctrl->current);
}


void
well_controls_iset_type( struct WellControls * ctrls , int control_index , enum WellControlType type) {
    ctrls->type[control_index] = type;
}


double
well_controls_iget_target(const struct WellControls * ctrl, int control_index) {
    return ctrl->target[control_index];
}

double
well_controls_get_current_target(const struct WellControls * ctrl) {
    return ctrl->target[ctrl->current];
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


const double *
well_controls_get_current_distr(const struct WellControls * ctrl) {
    return well_controls_iget_distr( ctrl , ctrl->current );
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

    well_controls_iset_type( ctrl , ctrl->num , type);
    well_controls_iset_target( ctrl , ctrl->num , target);
    
    if (distr != NULL) 
        well_controls_iset_distr( ctrl , ctrl->num , distr);

    ctrl->num += 1;
    return 1;
}


bool
well_controls_equal(const struct WellControls *ctrls1, const struct WellControls *ctrls2 , bool verbose)
/* ---------------------------------------------------------------------- */
{
    bool are_equal = true;

    if (ctrls1->num !=  ctrls2->num) {
        are_equal = false;
        if (verbose)
            printf("ctrls1->num:%d    ctrls2->num:%d \n",ctrls1->num , ctrls2->num);
    }

    if (ctrls1->number_of_phases !=  ctrls2->number_of_phases) {
        are_equal = false;
        if (verbose)
            printf("ctrls1->number_of_phases:%d    ctrls2->number_of_phases:%d \n",ctrls1->number_of_phases , ctrls2->number_of_phases);
    }

    if (!are_equal) {
        return are_equal;
    }

    if (memcmp(ctrls1->type, ctrls2->type, ctrls1->num * sizeof *ctrls1->type ) != 0) {
        are_equal = false;
        if (verbose) 
            printf("The ->type vectors are different \n");
    }

    if (memcmp(ctrls1->target, ctrls2->target, ctrls1->num * sizeof *ctrls1->target ) != 0) {
        are_equal = false;
        if (verbose) 
            printf("The ->target vectors are different \n");
    }

    return are_equal;
}

