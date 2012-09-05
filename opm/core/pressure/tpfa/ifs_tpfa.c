#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/newwells.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>


/* ---------------------------------------------------------------------- */
static void
mult_csr_matrix(const struct CSRMatrix *A,
                const double           *u,
                double                 *v)
/* ---------------------------------------------------------------------- */
{
    int    j;
    size_t i;

    for (i = 0, j = 0; i < A->m; i++) {
        v[i] = 0.0;

        for (; j < A->ia[i + 1]; j++) {
            v[i] += A->sa[j] * u[ A->ja[j] ];
        }
    }
}


struct ifs_tpfa_impl {
    double *fgrav;              /* Accumulated grav contrib/face */
    double *work;

    /* Linear storage */
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
impl_deallocate(struct ifs_tpfa_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        free(pimpl->ddata);
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
static struct ifs_tpfa_impl *
impl_allocate(struct UnstructuredGrid *G,
              struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_impl *new;

    size_t nnu;
    size_t ddata_sz;

    nnu = G->number_of_cells;
    if (W != NULL) {
        nnu += W->number_of_wells;
    }

    ddata_sz  = 2 * nnu;                 /* b, x */
    ddata_sz += 1 * G->number_of_faces;  /* fgrav */
    ddata_sz += 1 * nnu;                 /* work */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->ddata = malloc(ddata_sz * sizeof *new->ddata);

        if (new->ddata == NULL) {
            impl_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct CSRMatrix *
ifs_tpfa_construct_matrix(struct UnstructuredGrid *G,
                          struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2, w, i, nc, nnu;
    size_t nnz;

    struct CSRMatrix *A;

    nc = nnu = G->number_of_cells;
    if (W != NULL) {
        nnu += W->number_of_wells;
    }

    A = csrmatrix_new_count_nnz(nnu);

    if (A != NULL) {
        /* Self connections */
        for (c1 = 0; c1 < nnu; c1++) {
            A->ia[ c1 + 1 ] = 1;
        }

        /* Other connections */
        for (f = 0; f < G->number_of_faces; f++) {
            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            if ((c1 >= 0) && (c2 >= 0)) {
                A->ia[ c1 + 1 ] += 1;
                A->ia[ c2 + 1 ] += 1;
            }
        }

        if (W != NULL) {
            /* Well <-> cell connections */
            for (w = i = 0; w < W->number_of_wells; w++) {
                for (; i < W->well_connpos[w + 1]; i++) {
                    c1 = W->well_cells[i];

                    A->ia[ 0  + c1 + 1 ] += 1; /* c -> w */
                    A->ia[ nc + w  + 1 ] += 1; /* w -> c */
                }
            }
        }

        nnz = csrmatrix_new_elms_pushback(A);
        if (nnz == 0) {
            csrmatrix_delete(A);
            A = NULL;
        }
    }

    if (A != NULL) {
        /* Fill self connections */
        for (i = 0; i < nnu; i++) {
            A->ja[ A->ia[ i + 1 ] ++ ] = i;
        }

        /* Fill other connections */
        for (f = 0; f < G->number_of_faces; f++) {
            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            if ((c1 >= 0) && (c2 >= 0)) {
                A->ja[ A->ia[ c1 + 1 ] ++ ] = c2;
                A->ja[ A->ia[ c2 + 1 ] ++ ] = c1;
            }
        }

        if (W != NULL) {
            /* Fill well <-> cell connections */
            for (w = i = 0; w < W->number_of_wells; w++) {
                for (; i < W->well_connpos[w + 1]; i++) {
                    c1 = W->well_cells[i];

                    A->ja[ A->ia[ 0  + c1 + 1 ] ++ ] = nc + w;
                    A->ja[ A->ia[ nc + w  + 1 ] ++ ] = c1    ;
                }
            }
        }

        assert ((size_t) A->ia[ nnu ] == nnz);

        /* Guarantee sorted rows */
        csrmatrix_sortrows(A);
    }

    return A;
}


/* ---------------------------------------------------------------------- */
/* fgrav = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]) */
/* ---------------------------------------------------------------------- */
static void
compute_grav_term(struct UnstructuredGrid *G, const double *gpress,
                  double *fgrav)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f, c1, c2;
    double s;

    vector_zero(G->number_of_faces, fgrav);

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f  = G->cell_faces[i];

            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            s  = 2.0*(c1 == c) - 1.0;

            if ((c1 >= 0) && (c2 >= 0)) {
                fgrav[f] += s * gpress[i];
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
assemble_bhp_well(int nc, int w,
                  const struct Wells   *W  ,
                  const double         *mt ,
                  const double         *wdp,
                  struct ifs_tpfa_data *h  )
/* ---------------------------------------------------------------------- */
{
    int    c, i, wdof;
    size_t jc, jw;
    double trans, bhp;

    struct WellControls *ctrls;

    ctrls = W->ctrls[ w ];
    wdof  = nc + w;
    bhp   = ctrls->target[ ctrls->current ];

    jw    = csrmatrix_elm_index(wdof, wdof, h->A);

    for (i = W->well_connpos[w]; i < W->well_connpos[w + 1]; i++) {

        c     = W->well_cells  [ i ];
        trans = mt[ c ] * W->WI[ i ];

        jc = csrmatrix_elm_index(c, c, h->A);

        /* c<->c diagonal contribution from well */
        h->A->sa[ jc   ] += trans;
        h->b    [ c    ] += trans * (bhp + wdp[ i ]);

        /* w<->w diagonal contribution from well, trivial eqn. */
        h->A->sa[ jw   ] += trans;
        h->b    [ wdof ] += trans * bhp;
    }
}


/* ---------------------------------------------------------------------- */
static void
assemble_rate_well(int nc, int w,
                   const struct Wells   *W  ,
                   const double         *mt ,
                   const double         *wdp,
                   struct ifs_tpfa_data *h  )
/* ---------------------------------------------------------------------- */
{
    int    c, i, wdof;
    size_t jcc, jcw, jwc, jww;
    double trans, resv;

    struct WellControls *ctrls;

    ctrls = W->ctrls[ w ];
    wdof  = nc + w;
    resv  = ctrls->target[ ctrls->current ];

    jww   = csrmatrix_elm_index(wdof, wdof, h->A);

    for (i = W->well_connpos[w]; i < W->well_connpos[w + 1]; i++) {

        c   = W->well_cells[ i ];

        jcc = csrmatrix_elm_index(c   , c   , h->A);
        jcw = csrmatrix_elm_index(c   , wdof, h->A);
        jwc = csrmatrix_elm_index(wdof, c   , h->A);

        /* Connection transmissibility */
        trans = mt[ c ] * W->WI[ i ];

        /* c->w connection */
        h->A->sa[ jcc  ] += trans;
        h->A->sa[ jcw  ] -= trans;
        h->b    [ c    ] += trans * wdp[ i ];

        /* w->c connection */
        h->A->sa[ jwc  ] -= trans;
        h->A->sa[ jww  ] += trans;
        h->b    [ wdof ] -= trans * wdp[ i ];
    }

    h->b[ wdof ] += resv;
}


/* ---------------------------------------------------------------------- */
static void
assemble_shut_well(int nc, int w,
                   const struct Wells   *W  ,
                   const double         *mt ,
                   struct ifs_tpfa_data *h  )
/* ---------------------------------------------------------------------- */
{
    int    c, i, wdof;
    size_t jw;
    double trans;

    /* The equation added for a shut well w is
     *   \sum_{c~w} T_{w,c} * mt_c p_w = 0.0
     * where c~w indicates the cell perforated by w, T are production
     * indices, mt is total mobility and p_w is the bottom-hole
     * pressure.
     * Note: post-processing in ifs_tpfa_press_flux() may change the
     * well bhp away from zero, the equation used here is intended
     * just to
     *    1) give the same solution as if the well was not there for
     *       all other degrees of freedom,
     *    2) not disturb the conditioning of the matrix.
     */

    wdof  = nc + w;

    jw    = csrmatrix_elm_index(wdof, wdof, h->A);

    for (i = W->well_connpos[w]; i < W->well_connpos[w + 1]; i++) {

        c     = W->well_cells  [ i ];
        trans = mt[ c ] * W->WI[ i ];

        /* w<->w diagonal contribution from well, trivial eqn. */
        h->A->sa[ jw ] += trans;
    }
    h->b[ wdof ] = 0.0;
}


/* ---------------------------------------------------------------------- */
static void
assemble_well_contrib(int                   nc ,
                      const struct Wells   *W  ,
                      const double         *mt ,
                      const double         *wdp,
                      struct ifs_tpfa_data *h   ,
                      int                  *all_rate,
                      int                  *ok)
/* ---------------------------------------------------------------------- */
{
    int w, p, np;

    struct WellControls *ctrls;

    np = W->number_of_phases;

    *all_rate = 1;
    *ok = 1;

    for (w = 0; w < W->number_of_wells; w++) {
        ctrls = W->ctrls[ w ];

        if (ctrls->current < 0) {

            /* Threat this well as a shut well, isolated from the domain. */

            assemble_shut_well(nc, w, W, mt, h);

        } else {

            assert (ctrls->current <  ctrls->num);

            switch (ctrls->type[ ctrls->current ]) {
            case BHP:
                *all_rate = 0;
                assemble_bhp_well (nc, w, W, mt, wdp, h);
                break;

            case RESERVOIR_RATE:
                for (p = 0; p < np; p++) {
                    if (ctrls->distr[np * ctrls->current + p] != 1.0) {
                        *ok = 0;
                        break;
                    }
                }
                if (*ok) {
                    assemble_rate_well(nc, w, W, mt, wdp, h);
                }
                break;

            case SURFACE_RATE:
                /* We cannot handle this case, since we do
                 * not have access to formation volume factors,
                 * needed to convert between reservoir and
                 * surface rates
                 */
                fprintf(stderr, "ifs_tpfa cannot handle SURFACE_RATE well controls.\n");
                *ok = 0;
                break;
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
assemble_bc_contrib(struct UnstructuredGrid             *G    ,
                    const struct FlowBoundaryConditions *bc   ,
                    const double                        *trans,
                    struct ifs_tpfa_data                *h    )
/* ---------------------------------------------------------------------- */
{
    int is_neumann, is_outflow;
    int f, c1, c2;

    size_t i, j, ix;
    double s, t;

    is_neumann = 1;

    for (i = 0; i < bc->nbc; i++) {
        if (bc->type[ i ] == BC_PRESSURE) {
            is_neumann = 0;

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[i + 1]; j++) {
                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0)); /* BCs on ext. faces only */

                is_outflow = c1 >= 0;

                t  = trans[ f ];
                s  = 2.0*is_outflow - 1.0;
                c1 = is_outflow ? c1 : c2;
                ix = csrmatrix_elm_index(c1, c1, h->A);

                h->A->sa[ ix ] += t;
                h->b    [ c1 ] += t * bc->value[ i ];
                h->b    [ c1 ] -= s * t * h->pimpl->fgrav[ f ];
            }
        }

        else if (bc->type[ i ] == BC_FLUX_TOTVOL) {
            /* We currently support individual flux faces only. */
            assert (bc->cond_pos[i + 1] - bc->cond_pos[i] == 1);

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[i + 1]; j++) {
                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0)); /* BCs on ext. faces only */

                c1 = (c1 >= 0) ? c1 : c2;

                /* Interpret BC as flow *INTO* cell */
                h->b[ c1 ] += bc->value[ i ];
            }
        }

        /* Other types currently not handled */
    }

    return is_neumann;
}


/* ---------------------------------------------------------------------- */
static void
boundary_fluxes(struct UnstructuredGrid             *G     ,
                const struct FlowBoundaryConditions *bc    ,
                const double                        *trans ,
                const double                        *cpress,
                const struct ifs_tpfa_data          *h     ,
                double                              *fflux )
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2;
    size_t i, j;
    double s, dh;

    for (i = 0; i < bc->nbc; i++) {
        if (bc->type[ i ] == BC_PRESSURE) {
            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[ i + 1 ]; j++) {

                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0));

                if (c1 < 0) {   /* Environment -> c2 */
                    dh = bc->value[ i ] - cpress[c2];
                }
                else {          /* c1 -> environment */
                    dh = cpress[c1] - bc->value[ i ];
                }

                fflux[f] = trans[f] * (dh + h->pimpl->fgrav[f]);
            }
        }

        else if (bc->type[ i ] == BC_FLUX_TOTVOL) {
            assert (bc->cond_pos[i+1] - bc->cond_pos[i] == 1);

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[ i + 1 ]; j++) {

                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0));

                /* BC flux is positive into reservoir. */
                s = 2.0*(c1 < 0) - 1.0;

                fflux[f] = s * bc->value[ i ];
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
well_solution(const struct UnstructuredGrid *G   ,
              const struct ifs_tpfa_forces  *F   ,
              const struct ifs_tpfa_data    *h   ,
              struct ifs_tpfa_solution      *soln)
/* ---------------------------------------------------------------------- */
{
    int           c, w, i, shut;
    double        bhp, dp, trans;
    const double *mt, *WI, *wdp;

    if (soln->well_press != NULL) {
        /* Extract BHP directly from solution vector for non-shut wells */
        for (w = 0; w < F->W->number_of_wells; w++) {
            if (F->W->ctrls[w]->current >= 0) {
                soln->well_press[w] = h->x[G->number_of_cells + w];
            }
        }
    }

    if (soln->well_flux != NULL) {
        WI  = F->W->WI;
        mt  = F->totmob;
        wdp = F->wdp;

        for (w = i = 0; w < F->W->number_of_wells; w++) {
            bhp = h->x[G->number_of_cells + w];
            shut = (F->W->ctrls[w]->current < 0);

            for (; i < F->W->well_connpos[ w + 1 ]; i++) {

                c = F->W->well_cells[ i ];

                trans = mt[ c ] * WI[ i ];
                dp    = bhp + wdp[ i ] - soln->cell_press[ c ];

                soln->well_flux[i] = shut ? 0.0 : trans * dp;
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
assemble_incompressible(struct UnstructuredGrid      *G     ,
                        const struct ifs_tpfa_forces *F     ,
                        const double                 *trans ,
                        const double                 *gpress,
                        struct ifs_tpfa_data         *h     ,
                        int                          *singular,
                        int                          *ok    )
/* ---------------------------------------------------------------------- */
{
    int c1, c2, c, i, f, j1, j2;

    int res_is_neumann, wells_are_rate;

    double s;

    *ok = 1;
    csrmatrix_zero(         h->A);
    vector_zero   (h->A->m, h->b);

    compute_grav_term(G, gpress, h->pimpl->fgrav);

    for (c = i = 0; c < G->number_of_cells; c++) {
        j1 = csrmatrix_elm_index(c, c, h->A);

        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            s  = 2.0*(c1 == c) - 1.0;
            c2 = (c1 == c) ? c2 : c1;

            h->b[c] -= trans[f] * (s * h->pimpl->fgrav[f]);

            if (c2 >= 0) {
                j2 = csrmatrix_elm_index(c, c2, h->A);

                h->A->sa[j1] += trans[f];
                h->A->sa[j2] -= trans[f];
            }
        }
    }


    /* Assemble contributions from driving forces other than gravity */
    res_is_neumann = 1;
    wells_are_rate = 1;
    if (F != NULL) {
        if ((F->W != NULL) && (F->totmob != NULL) && (F->wdp != NULL)) {
            /* Contributions from wells */

            /* Ensure modicum of internal consistency. */
            assert (h->A->m ==
                    (size_t) G->number_of_cells +
                    (size_t) F->W->number_of_wells);

            assemble_well_contrib(G->number_of_cells, F->W,
                                  F->totmob, F->wdp, h,
                                  &wells_are_rate, ok);
        }

        if (F->bc != NULL) {
            /* Contributions from boundary conditions */
            res_is_neumann = assemble_bc_contrib(G, F->bc, trans, h);
        }

        if (F->src != NULL) {
            /* Contributions from explicit source terms. */
            for (c = 0; c < G->number_of_cells; c++) {
                h->b[c] += F->src[c];
            }
        }
    }

    *singular = res_is_neumann && wells_are_rate;
}

/* ======================================================================
 * Public interface below separator.
 * ====================================================================== */

/* ---------------------------------------------------------------------- */
struct ifs_tpfa_data *
ifs_tpfa_construct(struct UnstructuredGrid *G,
                   struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->pimpl = impl_allocate(G, W);
        new->A     = ifs_tpfa_construct_matrix(G, W);

        if ((new->pimpl == NULL) || (new->A == NULL)) {
            ifs_tpfa_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        new->b = new->pimpl->ddata;
        new->x = new->b                       + new->A->m;

        new->pimpl->fgrav = new->x            + new->A->m;
        new->pimpl->work  = new->pimpl->fgrav + G->number_of_faces;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
int
ifs_tpfa_assemble(struct UnstructuredGrid      *G     ,
                  const struct ifs_tpfa_forces *F     ,
                  const double                 *trans ,
                  const double                 *gpress,
                  struct ifs_tpfa_data         *h     )
/* ---------------------------------------------------------------------- */
{
    int system_singular, ok;

    assemble_incompressible(G, F, trans, gpress, h, &system_singular, &ok);

    if (ok && system_singular) {
        /* Remove zero eigenvalue associated to constant pressure */
        h->A->sa[0] *= 2.0;
    }
    return ok;
}


/* ---------------------------------------------------------------------- */
int
ifs_tpfa_assemble_comprock(struct UnstructuredGrid      *G        ,
                           const struct ifs_tpfa_forces *F        ,
                           const double                 *trans    ,
                           const double                 *gpress   ,
                           const double                 *porevol  ,
                           const double                 *rock_comp,
                           const double                  dt       ,
                           const double                 *pressure ,
                           struct ifs_tpfa_data         *h        )
/* ---------------------------------------------------------------------- */
{
    int    c, system_singular, ok;
    size_t j;
    double d;

    assemble_incompressible(G, F, trans, gpress, h, &system_singular, &ok);

    /*
     * The extra term of the equation is
     *
     *     porevol*rock_comp*(p - p0)/dt.
     *
     * The p part goes on the diagonal, the p0 on the rhs.
     * We don't care if the system was singular before this modification,
     * after it will always be nonsingular.
     */
    if (ok) {
        for (c = 0; c < G->number_of_cells; c++) {
            j = csrmatrix_elm_index(c, c, h->A);

            d = porevol[c] * rock_comp[c] / dt;

            h->A->sa[j] += d;
            h->b[c]     += d * pressure[c];
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
int
ifs_tpfa_assemble_comprock_increment(struct UnstructuredGrid      *G        ,
                                     const struct ifs_tpfa_forces *F        ,
                                     const double                 *trans    ,
                                     const double                 *gpress   ,
                                     const double                 *porevol  ,
                                     const double                 *rock_comp,
                                     const double                  dt       ,
                                     const double                 *prev_pressure,
                                     const double                 *initial_porevolume,
                                     struct ifs_tpfa_data         *h        )
/* ---------------------------------------------------------------------- */
{
    int     c, w, wdof, system_singular, ok;
    size_t  j;
    double *v, dpvdt;

    ok = 1;
    assemble_incompressible(G, F, trans, gpress, h, &system_singular, &ok);

    /* We want to solve a Newton step for the residual
     * (porevol(pressure)-porevol(initial_pressure))/dt + residual_for_incompressible
     *
     */

    if (ok) {
        v = h->pimpl->work;
        mult_csr_matrix(h->A, prev_pressure, v);

        for (c = 0; c < G->number_of_cells; c++) {
            j = csrmatrix_elm_index(c, c, h->A);

            dpvdt = (porevol[c] - initial_porevolume[c]) / dt;

            h->A->sa[j] += porevol[c] * rock_comp[c] / dt;
            h->b[c]     -= dpvdt + v[c];
        }

        if (F->W != NULL) {
            wdof = G->number_of_cells;

            for (w = 0; w < F->W->number_of_wells; w++, wdof++) {
                h->b[wdof] -= v[wdof];
            }
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_press_flux(struct UnstructuredGrid      *G    ,
                    const struct ifs_tpfa_forces *F    ,
                    const double                 *trans,
                    struct ifs_tpfa_data         *h    ,
                    struct ifs_tpfa_solution     *soln )
/* ---------------------------------------------------------------------- */
{
    int    c1, c2, f;
    double dh;

    double *cpress, *fflux;

    assert (soln             != NULL);
    assert (soln->cell_press != NULL);
    assert (soln->face_flux  != NULL);

    cpress = soln->cell_press;
    fflux  = soln->face_flux ;

    /* Assign cell pressure directly from solution vector */
    memcpy(cpress, h->x, G->number_of_cells * sizeof *cpress);

    for (f = 0; f < G->number_of_faces; f++) {
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if ((c1 >= 0) && (c2 >= 0)) {
            dh = cpress[c1] - cpress[c2] + h->pimpl->fgrav[f];
            fflux[f] = trans[f] * dh;
        } else {
            fflux[f] = 0.0;
        }
    }

    if (F != NULL) {
        if (F->bc != NULL) {
            boundary_fluxes(G, F->bc, trans, cpress, h, fflux);
        }

        if (F->W != NULL) {
            well_solution(G, F, h, soln);
        }
    }
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_destroy(struct ifs_tpfa_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        csrmatrix_delete(h->A);
        impl_deallocate (h->pimpl);
    }

    free(h);
}
