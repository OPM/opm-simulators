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

/* ---------------------------------------------------------------------- */
static void
well_controls_destroy(struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    if (ctrl != NULL) {
        free(ctrl->target);
        free(ctrl->type);
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
        ctrl->cpty    = 0;
        ctrl->type    = NULL;
        ctrl->target  = NULL;
        ctrl->current = -1;
    }

    return ctrl;
}


/* ---------------------------------------------------------------------- */
static int
well_controls_reserve(int nctrl, struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    int   c, ok;
    void *type, *target;

    type   = realloc(ctrl->type  , nctrl * sizeof *ctrl->type  );
    target = realloc(ctrl->target, nctrl * sizeof *ctrl->target);

    ok = 0;
    if (type   != NULL) { ctrl->type   = type  ; ok++; }
    if (target != NULL) { ctrl->target = target; ok++; }

    if (ok == 2) {
        for (c = ctrl->cpty; c < nctrl; c++) {
            ctrl->type  [c] =  BHP;
            ctrl->target[c] = -1.0;
        }

        ctrl->cpty = nctrl;
    }

    return ok == 2;
}


/* ---------------------------------------------------------------------- */
static int
wells_allocate(int nwells, struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int   ok;
    void *type, *depth_ref, *zfrac;
    void *well_connpos;
    void *ctrls;

    type      = realloc(W->type     , 1 * nwells * sizeof *W->type);
    depth_ref = realloc(W->depth_ref, 1 * nwells * sizeof *W->depth_ref);
    zfrac     = realloc(W->zfrac    , 3 * nwells * sizeof *W->zfrac);
    ctrls     = realloc(W->ctrls    , 1 * nwells * sizeof *W->ctrls);

    well_connpos = realloc(W->well_connpos,
                           (nwells + 1) * sizeof *W->well_connpos);

    ok = 0;
    if (type         != NULL) { W->type         = type        ; ok++; }
    if (depth_ref    != NULL) { W->depth_ref    = depth_ref   ; ok++; }
    if (zfrac        != NULL) { W->zfrac        = zfrac       ; ok++; }
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
    int ok, w;

    for (w = W->well_cpty; w < nwells; w++) {
        W->type     [w]        = PRODUCER;
        W->depth_ref[w]        = -1.0;

        W->zfrac[3*w + WATER]  = 0.0;
        W->zfrac[3*w + OIL  ]  = 0.0;
        W->zfrac[3*w + GAS  ]  = 0.0;

        W->well_connpos[w + 1] = W->well_connpos[w];
    }

    for (w = W->well_cpty, ok = 1; ok && (w < nwells); w++) {
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

    for (k = W->perf_cpty; k < nperf; k++) {
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

    assert (nwells >= W->well_cpty);
    assert (nperf  >= W->perf_cpty);

    ok = 1;

    if (nwells > W->well_cpty) {
        ok = wells_allocate(nwells, W);

        if (ok) {
            ok = initialise_new_wells(nwells, W);
        }

        if (ok) {
            W->well_cpty = nwells;
        }
    }

    if (ok && (nperf > W->perf_cpty)) {
        ok = perfs_allocate(nperf, W);

        if (ok) {
            initialise_new_perfs(nperf, W);
            W->perf_cpty = nperf;
        }
    }

    return ok;
}


/* ======================================================================
 * Public entry points below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct Wells *
wells_create(int nwells, int nperf)
/* ---------------------------------------------------------------------- */
{
    int           ok;
    struct Wells *W;

    W = malloc(1 * sizeof *W);

    if (W != NULL) {
        W->number_of_wells = 0;
        W->well_cpty       = 0;
        W->perf_cpty       = 0;

        W->type            = NULL;
        W->depth_ref       = NULL;
        W->zfrac           = NULL;

        W->well_connpos    = malloc(1 * sizeof *W->well_connpos);
        W->well_cells      = NULL;
        W->WI              = NULL;

        W->ctrls           = NULL;

        ok = W->well_connpos != NULL;
        if (ok) {
            W->well_connpos[0] = 0;

            if ((nwells > 0) || (nperf > 0)) {
                ok = wells_reserve(nwells, nperf, W);
            }
        }

        if (! ok) {
            wells_destroy(W);
            W = NULL;
        }
    }

    return W;
}


/* ---------------------------------------------------------------------- */
void
wells_destroy(struct Wells *W)
/* ---------------------------------------------------------------------- */
{
    int w;

    if (W != NULL) {
        for (w = 0; w < W->well_cpty; w++) {
            well_controls_destroy(W->ctrls[w]);
        }

        free(W->ctrls);
        free(W->WI);
        free(W->well_cells);
        free(W->well_connpos);
        free(W->zfrac);
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
wells_add(enum well_type type     ,
          double         depth_ref,
          int            nperf    ,
          const double  *zfrac    , /* Injection fraction or NULL */
          const int     *cells    ,
          const double  *WI       , /* Well index per perf (or NULL) */
          struct Wells  *W        )
/* ---------------------------------------------------------------------- */
{
    int ok, nw, nperf_tot, off;
    int nwalloc, nperfalloc;

    nw        = W->number_of_wells;
    nperf_tot = W->well_connpos[nw];

    ok = (nw < W->well_cpty) && (nperf_tot + nperf <= W->perf_cpty);

    if (! ok) {
        nwalloc    = alloc_size(nw       , 1    , W->well_cpty);
        nperfalloc = alloc_size(nperf_tot, nperf, W->perf_cpty);

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

        if (zfrac != NULL) {
            memcpy(W->zfrac + 3*nw, zfrac, 3 * sizeof *W->zfrac);
        }

        W->well_connpos[nw + 1]  = off + nperf;
        W->number_of_wells      += 1;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
int
well_controls_append(enum control_type    type  ,
                     double               target,
                     struct WellControls *ctrl  )
/* ---------------------------------------------------------------------- */
{
    int ok, alloc;

    assert (ctrl != NULL);

    ok = ctrl->num < ctrl->cpty;

    if (! ok) {
        alloc = alloc_size(ctrl->num, 1, ctrl->cpty);
        ok    = well_controls_reserve(alloc, ctrl);
    }

    if (ok) {
        ctrl->type  [ctrl->num] = type  ;
        ctrl->target[ctrl->num] = target;

        ctrl->num += 1;
        
        /* TODO: Review this: */
        ctrl->current = 0;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
void
well_controls_clear(struct WellControls *ctrl)
/* ---------------------------------------------------------------------- */
{
    if (ctrl != NULL) {
        ctrl->num = 0;
    }
}
