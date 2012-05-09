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

#include <opm/core/newwells.h>

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

struct WellControlMgmt {
    int cpty;
};

struct WellMgmt {
    int well_cpty;
    int perf_cpty;
};


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


static void
destroy_well_mgmt(struct WellMgmt *m)
{
    free(m);
}


static struct WellMgmt *
create_well_mgmt(void)
{
    struct WellMgmt *m;

    m = malloc(1 * sizeof *m);

    if (m != NULL) {
        m->well_cpty = 0;
        m->perf_cpty = 0;
    }

    return m;
}


/* ---------------------------------------------------------------------- */
static void
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
static struct WellControls *
well_controls_create(void)
/* ---------------------------------------------------------------------- */
{
    struct WellControls *ctrl;

    ctrl = malloc(1 * sizeof *ctrl);

    if (ctrl != NULL) {
        /* Initialise empty control set */
        ctrl->num     = 0;
        ctrl->type    = NULL;
        ctrl->target  = NULL;
        ctrl->distr   = NULL;
        ctrl->current = -1;

        ctrl->data    = create_ctrl_mgmt();

        if (ctrl->data == NULL) {
            well_controls_destroy(ctrl);
            ctrl = NULL;
        }
    }

    return ctrl;
}


/* ---------------------------------------------------------------------- */
static int
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

    if (ok == 3) {
        m = ctrl->data;
        for (c = m->cpty; c < nctrl; c++) {
            ctrl->type  [c] =  BHP;
            ctrl->target[c] = -1.0;
        }

        for (p = m->cpty * nphases; p < nctrl * nphases; ++p) {
            ctrl->distr[ p ] = 0.0;
        }

        m->cpty = nctrl;
    }

    return ok == 3;
}


/* ---------------------------------------------------------------------- */
static int
wells_allocate(int nwells, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int   ok, np;
    void *type, *depth_ref, *comp_frac;
    void *well_connpos;
    void *ctrls;

    np = W->number_of_phases;

    type      = realloc(W->type     ,  1 * nwells * sizeof *W->type);
    depth_ref = realloc(W->depth_ref,  1 * nwells * sizeof *W->depth_ref);
    comp_frac = realloc(W->comp_frac, np * nwells * sizeof *W->comp_frac);
    ctrls     = realloc(W->ctrls    ,  1 * nwells * sizeof *W->ctrls);

    well_connpos = realloc(W->well_connpos,
                           (nwells + 1) * sizeof *W->well_connpos);

    ok = 0;
    if (type         != NULL) { W->type         = type        ; ok++; }
    if (depth_ref    != NULL) { W->depth_ref    = depth_ref   ; ok++; }
    if (comp_frac    != NULL) { W->comp_frac    = comp_frac   ; ok++; }
    if (well_connpos != NULL) { W->well_connpos = well_connpos; ok++; }
    if (ctrls        != NULL) { W->ctrls        = ctrls       ; ok++; }

    return ok == 5;
}


/* ---------------------------------------------------------------------- */
static int
perfs_allocate(int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int   ok;
    void *well_cells, *WI;

    well_cells = realloc(W->well_cells, nperf * sizeof *W->well_cells);
    WI         = realloc(W->WI        , nperf * sizeof *W->WI        );

    ok = 0;
    if (well_cells != NULL) { W->well_cells = well_cells; ok++; }
    if (WI         != NULL) { W->WI         = WI        ; ok++; }

    return ok == 2;
}


/* ---------------------------------------------------------------------- */
static int
initialise_new_wells(int nwells, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int ok, w, p;

    struct WellMgmt *m;

    m = W->data;

    for (w = m->well_cpty; w < nwells; w++) {
        W->type     [w]        = PRODUCER;
        W->depth_ref[w]        = -1.0;

        for (p = 0; p < W->number_of_phases; ++p) {
            W->comp_frac[W->number_of_phases*w + p] = 0.0;
        }

        W->well_connpos[w + 1] = W->well_connpos[w];
    }

    for (w = m->well_cpty, ok = 1; ok && (w < nwells); w++) {
        W->ctrls[w] = well_controls_create();

        ok = W->ctrls[w] != NULL;
    }

    if (! ok) {
        for (; w < nwells; w++) {
            W->ctrls[w] = NULL;
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
initialise_new_perfs(int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int k;

    struct WellMgmt *m;

    m = W->data;

    for (k = m->perf_cpty; k < nperf; k++) {
        W->well_cells[k] = -1 ;
        W->WI        [k] = 0.0;
    }
}


/* ---------------------------------------------------------------------- */
static int
wells_reserve(int nwells, int nperf, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int ok;

    struct WellMgmt *m;

    m = W->data;

    assert (nwells >= m->well_cpty);
    assert (nperf  >= m->perf_cpty);

    ok = 1;

    if (nwells > m->well_cpty) {
        ok = wells_allocate(nwells, W);

        if (ok) {
            ok = initialise_new_wells(nwells, W);
        }

        if (ok) {
            m->well_cpty = nwells;
        }
    }

    if (ok && (nperf > m->perf_cpty)) {
        ok = perfs_allocate(nperf, W);

        if (ok) {
            initialise_new_perfs(nperf, W);
            m->perf_cpty = nperf;
        }
    }

    return ok;
}


/* ======================================================================
 * Public entry points below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct Wells *
create_wells(int nphases, int nwells, int nperf)
/* ---------------------------------------------------------------------- */
{
    int           ok;
    struct Wells *W;

    W = malloc(1 * sizeof *W);

    if (W != NULL) {
        W->number_of_wells = 0;
        W->number_of_phases = nphases;

        W->type            = NULL;
        W->depth_ref       = NULL;
        W->comp_frac       = NULL;

        W->well_connpos    = malloc(1 * sizeof *W->well_connpos);
        W->well_cells      = NULL;
        W->WI              = NULL;

        W->ctrls           = NULL;

        W->data            = create_well_mgmt();

        ok = (W->well_connpos != NULL) && (W->data != NULL);
        if (ok) {
            W->well_connpos[0] = 0;

            if ((nwells > 0) || (nperf > 0)) {
                ok = wells_reserve(nwells, nperf, W);
            }
        }

        if (! ok) {
            destroy_wells(W);
            W = NULL;
        }
    }

    return W;
}


/* ---------------------------------------------------------------------- */
void
destroy_wells(struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int w;

    struct WellMgmt *m;

    if (W != NULL) {
        m = W->data;

        for (w = 0; w < m->well_cpty; w++) {
            well_controls_destroy(W->ctrls[w]);
        }

        destroy_well_mgmt(m);

        free(W->ctrls);
        free(W->WI);
        free(W->well_cells);
        free(W->well_connpos);
        free(W->comp_frac);
        free(W->depth_ref);
        free(W->type);
    }

    free(W);
}


/* ---------------------------------------------------------------------- */
static int
alloc_size(int n, int a, int cpty)
/* ---------------------------------------------------------------------- */
{
    if (cpty < n + a) {
        cpty *= 2;              /* log_2(n) allocations */

        if (cpty < n + a) {     /* Typically for the first few allocs */
            cpty = n + a;
        }
    }

    return cpty;
}


/* ---------------------------------------------------------------------- */
int
add_well(enum WellType  type     ,
         double         depth_ref,
         int            nperf    ,
         const double  *comp_frac, /* Injection fraction or NULL */
         const int     *cells    ,
         const double  *WI       , /* Well index per perf (or NULL) */
         struct Wells  *W        )
/* ---------------------------------------------------------------------- */
{
    int ok, nw, np, nperf_tot, off;
    int nwalloc, nperfalloc;

    struct WellMgmt *m;

    nw        = W->number_of_wells;
    nperf_tot = W->well_connpos[nw];

    m = W->data;

    ok = (nw < m->well_cpty) && (nperf_tot + nperf <= m->perf_cpty);

    if (! ok) {
        nwalloc    = alloc_size(nw       , 1    , m->well_cpty);
        nperfalloc = alloc_size(nperf_tot, nperf, m->perf_cpty);

        ok = wells_reserve(nwalloc, nperfalloc, W);
    }

    off = W->well_connpos[nw];

    if (ok && (nperf > 0)) {
        assert (cells != NULL);

        memcpy(W->well_cells + off,
               cells, nperf * sizeof *W->well_cells);

        if (WI != NULL) {
            memcpy(W->WI + off, WI, nperf * sizeof *W->WI);
        }
    }

    if (ok) {
        W->type     [nw] = type     ;
        W->depth_ref[nw] = depth_ref;

        np = W->number_of_phases;
        if (comp_frac != NULL) {
            memcpy(W->comp_frac + np*nw, comp_frac, np * sizeof *W->comp_frac);
        }

        W->well_connpos[nw + 1]  = off + nperf;
        W->number_of_wells      += 1;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
int
append_well_controls(enum WellControlType type,
                     double               target,
                     const double        *distr,
                     int                  well_index,
                     struct Wells        *W)
/* ---------------------------------------------------------------------- */
{
    int ok, alloc, p, np;
    struct WellControls    *ctrl;
    struct WellControlMgmt *m;

    assert (W != NULL);

    ctrl = W->ctrls[well_index];
    np = W->number_of_phases;

    assert (ctrl != NULL);

    m  = ctrl->data;
    ok = ctrl->num < m->cpty;

    if (! ok) {
        alloc = alloc_size(ctrl->num, 1, m->cpty);
        ok    = well_controls_reserve(alloc, np, ctrl);
    }

    if (ok) {
        ctrl->type  [ctrl->num] = type  ;
        ctrl->target[ctrl->num] = target;

        if (distr != NULL) {
            memcpy(ctrl->distr + (ctrl->num * np), distr,
                   np * sizeof *ctrl->distr);
        }

        ctrl->num += 1;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
void
set_current_control(int well_index, int current_control, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    assert (W->ctrls[well_index] != NULL);

    W->ctrls[well_index]->current = current_control;
}


/* ---------------------------------------------------------------------- */
void
clear_well_controls(int well_index, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    if (W->ctrls[well_index] != NULL) {
        W->ctrls[well_index]->num = 0;
    }
}
