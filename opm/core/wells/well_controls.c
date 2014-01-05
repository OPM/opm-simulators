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

#include <string.h>
#include <stdlib.h>


static void
destroy_ctrl_mgmt(struct WellControlMgmt *m)
{
    free(m);
}


static struct WellControlMgmt *
create_ctrl_mgmt(void)
{
    struct WellControlMgmt *m;

    m = malloc(1 * sizeof *m);

    if (m != NULL) {
        m->cpty = 0;
    }

    return m;
}



/* ---------------------------------------------------------------------- */
void
well_controls_destroy(struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    if (ctrl != NULL) {
        destroy_ctrl_mgmt(ctrl->data);
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

        ctrl->data              = create_ctrl_mgmt();

        if (ctrl->data == NULL) {
            well_controls_destroy(ctrl);
            ctrl = NULL;
        }
    }

    return ctrl;
}


/* ---------------------------------------------------------------------- */
int
well_controls_reserve(int nctrl, int nphases, struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    int   c, p, ok;
    void *type, *target, *distr;

    struct WellControlMgmt *m;

    type   = realloc(ctrl->type  , nctrl * 1       * sizeof *ctrl->type  );
    target = realloc(ctrl->target, nctrl * 1       * sizeof *ctrl->target);
    distr  = realloc(ctrl->distr , nctrl * nphases * sizeof *ctrl->distr );

    ok = 0;
    if (type   != NULL) { ctrl->type   = type  ; ok++; }
    if (target != NULL) { ctrl->target = target; ok++; }
    if (distr  != NULL) { ctrl->distr  = distr ; ok++; }

    ctrl->number_of_phases = nphases;

    if (ok == 3) {
        m = ctrl->data;
        for (c = m->cpty; c < nctrl; c++) {
            ctrl->type  [c] =  BHP;
            ctrl->target[c] = -1.0;
        }

        for (p = m->cpty * ctrl->number_of_phases; p < nctrl * ctrl->number_of_phases; ++p) {
            ctrl->distr[ p ] = 0.0;
        }

        m->cpty = nctrl;
    }

    return ok == 3;
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

    struct WellControlMgmt* mgmt1 = (struct WellControlMgmt*)(ctrls1->data);
    struct WellControlMgmt* mgmt2 = (struct WellControlMgmt*)(ctrls2->data);
    are_equal &= (mgmt1->cpty == mgmt2->cpty);

    return are_equal;
}
