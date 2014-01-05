#include "config.h"
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define HAVE_WELLCONTROLS
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/pressure/tpfa/compr_quant_general.h>
#include <opm/core/pressure/tpfa/compr_source.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/pressure/tpfa/cfs_tpfa_residual.h>

#if defined(MAX)
#undef MAX
#endif

#define MAX(a,b) (((a) > (b)) ? (a) : (b))



struct densrat_util {
    MAT_SIZE_T *ipiv;

    double      residual;
    double     *lu;
    double     *t1;
    double     *t2;
    double     *mat_row;
    double     *coeff;
    double     *linsolve_buffer;
};


struct cfs_tpfa_res_impl {
    int                  is_incomp;

    /* One entry per component per face */
    double              *compflux_f;       /* A_{ij} v_{ij} */
    double              *compflux_deriv_f; /* A_{ij} \partial_{p} v_{ij} */

    /* One entry per component per perforation */
    double              *compflux_p;       /* A_{wi} q_{wi} */
    double              *compflux_deriv_p; /* A_{wi} \partial_{p} q_{wi} */

    double              *flux_work;

    /* Scratch array for face pressure calculation */
    double              *scratch_f;

    struct densrat_util *ratio;

    /* Linear storage */
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
deallocate_densrat(struct densrat_util *ratio)
/* ---------------------------------------------------------------------- */
{
    if (ratio != NULL) {
        free(ratio->lu);
        free(ratio->ipiv);
    }

    free(ratio);
}


/* ---------------------------------------------------------------------- */
static struct densrat_util *
allocate_densrat(size_t max_conn, int np)
/* ---------------------------------------------------------------------- */
{
    size_t               alloc_sz, n_buffer_col;
    struct densrat_util *ratio;

    ratio = malloc(1 * sizeof *ratio);

    if (ratio != NULL) {
        n_buffer_col  = 1;            /* z */
        n_buffer_col += 1 * max_conn; /* A_{ij} v_{ij} */
        n_buffer_col += 2 * max_conn; /* A_{ij} \partial_{p} v_{ij} */

        alloc_sz  = np             * np; /* lu */
        alloc_sz += 2              * np; /* t1, t2 */
        alloc_sz += (max_conn + 1) * 1 ; /* mat_row */
        alloc_sz += (max_conn + 1) * 1 ; /* coeff */
        alloc_sz += n_buffer_col   * np; /* linsolve_buffer */

        ratio->ipiv = malloc(np       * sizeof *ratio->ipiv);
        ratio->lu   = malloc(alloc_sz * sizeof *ratio->lu  );

        if ((ratio->ipiv == NULL) || (ratio->lu == NULL)) {
            deallocate_densrat(ratio);
            ratio = NULL;
        } else {
            ratio->t1              = ratio->lu      + (np             * np);
            ratio->t2              = ratio->t1      + (1              * np);
            ratio->mat_row         = ratio->t2      + (1              * np);
            ratio->coeff           = ratio->mat_row + ((max_conn + 1) * 1 );
            ratio->linsolve_buffer = ratio->coeff   + ((max_conn + 1) * 1 );
        }
    }

    return ratio;
}


/* ---------------------------------------------------------------------- */
static void
impl_deallocate(struct cfs_tpfa_res_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        free              (pimpl->ddata);
        deallocate_densrat(pimpl->ratio);
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
static struct cfs_tpfa_res_impl *
impl_allocate(struct UnstructuredGrid   *G       ,
              struct cfs_tpfa_res_wells *wells   ,
              size_t                     max_conn,
              int                        np      )
/* ---------------------------------------------------------------------- */
{
    size_t                nnu, nwperf;
    struct cfs_tpfa_res_impl *new;

    size_t ddata_sz;

    nnu    = G->number_of_cells;
    nwperf = 0;

    if ((wells != NULL) && (wells->W != NULL)) {
        nnu    += wells->W->number_of_wells;
        nwperf  = wells->W->well_connpos[ wells->W->number_of_wells ];
    }

    /* Linear system */
    ddata_sz  = nnu;            /* b */

    /* Reservoir */
    ddata_sz += np *      G->number_of_faces ; /* compflux_f */
    ddata_sz += np * (2 * G->number_of_faces); /* compflux_deriv_f */

    /* Well perforations */
    ddata_sz += np *      nwperf ;             /* compflux_p */
    ddata_sz += np * (2 * nwperf);             /* compflux_deriv_p */

    ddata_sz += np * (1 + 2)                 ; /* flux_work */

    ddata_sz += 1  *      G->number_of_faces ; /* scratch_f */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->ddata = malloc(ddata_sz * sizeof *new->ddata);
        new->ratio = allocate_densrat(max_conn, np);

        if (new->ddata == NULL || new->ratio == NULL) {
            impl_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct CSRMatrix *
construct_matrix(struct UnstructuredGrid   *G    ,
                 struct cfs_tpfa_res_wells *wells)
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2, w, i, nc, nnu;
    int    wells_present;
    size_t nnz;

    struct CSRMatrix *A;

    nc = nnu = G->number_of_cells;

    wells_present = (wells != NULL) && (wells->W != NULL);
    if (wells_present) {
        nnu += wells->W->number_of_wells;
    }

    A = csrmatrix_new_count_nnz(nnu);

    if (A != NULL) {
        /* Self connections */
        for (i = 0; i < nnu; i++) {
            A->ia[ i + 1 ] = 1;
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

        if (wells_present) {
            /* Well <-> cell connections */
            struct Wells *W = wells->W;

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

        if (wells_present) {
            /* Fill well <-> cell connections */
            struct Wells *W = wells->W;

            for (w = i = 0; w < W->number_of_wells; w++) {
                for (; i < W->well_connpos[w + 1]; i++) {
                    c1 = W->well_cells[i];

                    A->ja[ A->ia[ 0  + c1 + 1 ] ++ ] = nc + w;
                    A->ja[ A->ia[ nc + w  + 1 ] ++ ] = c1    ;
                }
            }
        }

        assert ((size_t) A->ia[ nnu ] == nnz);

        /* Enforce sorted connection structure per row */
        csrmatrix_sortrows(A);
    }

    return A;
}


static void
factorise_fluid_matrix(int np, const double *A, struct densrat_util *ratio)
{
    int        np2;
    MAT_SIZE_T m, n, ld, info;

    m = n = ld = np;
    np2 = np * np;

    memcpy (ratio->lu, A, np2 * sizeof *ratio->lu);
    dgetrf_(&m, &n, ratio->lu, &ld, ratio->ipiv, &info);

    assert (info == 0);
}


static void
solve_linear_systems(int                  np   ,
                     MAT_SIZE_T           nrhs ,
                     struct densrat_util *ratio,
                     double              *b    )
{
    MAT_SIZE_T n, ldA, ldB, info;

    n = ldA = ldB = np;

    dgetrs_("No Transpose", &n,
            &nrhs, ratio->lu, &ldA, ratio->ipiv,
            b               , &ldB, &info);

    assert (info == 0);
}


static void
matvec(int nrow, int ncol, const double *A, const double *x, double *y)
{
    MAT_SIZE_T m, n, ld, incx, incy;
    double     a1, a2;

    m    = ld = nrow;
    n    = ncol;
    incx = incy = 1;
    a1   = 1.0;
    a2   = 0.0;

    dgemv_("No Transpose", &m, &n,
           &a1, A, &ld, x, &incx,
           &a2,         y, &incy);
}


static void
matmat(int np, int ncol, const double *A, const double *B, double *C)
{
    MAT_SIZE_T m, n, k, ldA, ldB, ldC;
    double     a1, a2;

    m  = k = ldA = ldB = ldC = np;
    n  = ncol;
    a1 = 1.0;
    a2 = 0.0;

    dgemm_("No Transpose", "No Transpose", &m, &n, &k,
           &a1, A, &ldA, B, &ldB, &a2, C, &ldC);
}


static void
compute_darcyflux_and_deriv(int           np,
                            double        trans,
                            double        dp,
                            const double *pmobf,
                            const double *gcapf,
                            double       *dflux,
                            double       *dflux_deriv)
{
    int    p;
    double a;

    for (p = 0; p < np; p++) {
        a = trans * pmobf[p];

        dflux      [       p] =   a * (dp + gcapf[p]);
        dflux_deriv[0*np + p] =   a; /* ignore gravity... */
        dflux_deriv[1*np + p] = - a;
    }
}


static void
compute_compflux_and_deriv(struct UnstructuredGrid  *G     ,
                           int                       np    ,
                           const double             *cpress,
                           const double             *trans ,
                           const double             *pmobf ,
                           const double             *gcapf ,
                           const double             *Af    ,
                           struct cfs_tpfa_res_impl *pimpl )
{
    int     c1, c2, f, np2;
    double  dp;
    double *cflux, *dcflux;

    np2    = np * np;

    cflux  = pimpl->compflux_f;
    dcflux = pimpl->compflux_deriv_f;

    for (f = 0; f < G->number_of_faces;
         f++, pmobf += np, gcapf += np, Af += np2,
             cflux += np, dcflux += 2 * np) {

        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if ((c1 >= 0) && (c2 >= 0)) {
            dp = cpress[c1] - cpress[c2];

            compute_darcyflux_and_deriv(np, trans[f], dp, pmobf, gcapf,
                                        pimpl->flux_work,
                                        pimpl->flux_work + np);

            /* Component flux = Af * v*/
            matvec(np, np, Af, pimpl->flux_work     , cflux );

            /* Derivative = Af * (dv/dp) */
            matmat(np, 2 , Af, pimpl->flux_work + np, dcflux);
        }

        /* Boundary connections excluded */
    }
}


static void
compute_well_compflux_and_deriv(struct cfs_tpfa_res_wells *wells ,
                                int                        np    ,
                                const double              *cpress,
                                const double              *wpress,
                                struct cfs_tpfa_res_impl  *pimpl )
{
    int           c, i, w, np2;
    double        pw, dp;
    const double *WI, *wdp, *Ap, *pmobp;
    double       *pflux, *dpflux, gpot[3] = { 0.0 };

    struct Wells *W;

    assert (wells->W != NULL);
    assert (wells->W->number_of_phases <= 3);

    W = wells->W;

    assert (W->ctrls != NULL);
    assert (W->data  != NULL);

    WI     = W->WI;
    wdp    = wells->data->wdp;
    Ap     = wells->data->A;
    pmobp  = wells->data->phasemob;

    np2    = np * np;

    pflux  = pimpl->compflux_p;
    dpflux = pimpl->compflux_deriv_p;

    for (w = i = 0; w < W->number_of_wells; w++) {
        pw = wpress[w];

        for (; i < W->well_connpos[w + 1]; i++,
                 Ap += np2, pmobp += np,
                 pflux += np, dpflux += 2 * np) {

            c  = W->well_cells[i];
            dp = pw + wdp[i]- cpress[c];

            compute_darcyflux_and_deriv(np, WI[i], dp, pmobp, gpot,
                                        pimpl->flux_work,
                                        pimpl->flux_work + np);

            /* Component flux = Ap * q*/
            matvec(np, np, Ap, pimpl->flux_work     , pflux );

            /* Derivative = Ap * (dq/dp) */
            matmat(np, 2 , Ap, pimpl->flux_work + np, dpflux);
        }
    }
}



static int
count_internal_conn(struct UnstructuredGrid *G, int c)
{
    int c1, c2, f, i, nconn;

    nconn = 0;

    for (i = G->cell_facepos[c]; i < G->cell_facepos[c + 1]; i++) {
        f  = G->cell_faces[i];
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        nconn += (c1 >= 0) && (c2 >= 0);
    }

    return nconn;
}


static int
init_cell_contrib(struct UnstructuredGrid  *G    ,
                  int                       c    ,
                  int                       np   ,
                  double                    pvol ,
                  double                    dt   ,
                  const double             *z    ,
                  struct cfs_tpfa_res_impl *pimpl)
{
    int     c1, c2, f, i, conn, nconn;
    double *cflx, *dcflx;

    nconn = count_internal_conn(G, c);

    memcpy(pimpl->ratio->linsolve_buffer, z, np * sizeof *z);

    pimpl->ratio->coeff[0] = -pvol;
    conn = 1;

    cflx  = pimpl->ratio->linsolve_buffer + (1 * np);
    dcflx = cflx + (nconn * np);

    for (i = G->cell_facepos[c]; i < G->cell_facepos[c + 1]; i++) {
        f  = G->cell_faces[i];
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if ((c1 >= 0) && (c2 >= 0)) {
            memcpy(cflx, pimpl->compflux_f + (f*np + 0),
                   np * sizeof *cflx);

            memcpy(dcflx, pimpl->compflux_deriv_f + (f*(2 * np) + 0),
                   2 * np * sizeof *dcflx);

            cflx  += 1 * np;
            dcflx += 2 * np;

            pimpl->ratio->coeff[ conn++ ] = dt * (2*(c1 == c) - 1.0);
        }
    }

    assert (conn == nconn + 1);
    assert (cflx == pimpl->ratio->linsolve_buffer + (nconn + 1)*np);

    return nconn;
}


static void
compute_cell_contrib(struct UnstructuredGrid  *G    ,
                     int                       c    ,
                     int                       np   ,
                     double                    pvol ,
                     double                    dt   ,
                     const double             *z    ,
                     const double             *Ac   ,
                     const double             *dAc  ,
                     struct cfs_tpfa_res_impl *pimpl)
{
    int        c1, c2, f, i, off, nconn, p;
    MAT_SIZE_T nrhs;
    double     s, dF1, dF2, *dv, *dv1, *dv2;

    nconn = init_cell_contrib(G, c, np, pvol, dt, z, pimpl);
    nrhs  = 1 + (1 + 2)*nconn;  /* [z, Af*v, Af*dv] */

    factorise_fluid_matrix(np, Ac, pimpl->ratio);
    solve_linear_systems  (np, nrhs, pimpl->ratio,
                           pimpl->ratio->linsolve_buffer);

    /* Sum residual contributions over the connections (+ accumulation):
     *   t1 <- (Ac \ [z, Af*v]) * [-pvol; repmat(dt, [nconn, 1])] */
    matvec(np, nconn + 1, pimpl->ratio->linsolve_buffer,
           pimpl->ratio->coeff, pimpl->ratio->t1);

    /* Compute residual in cell 'c' */
    pimpl->ratio->residual = pvol;
    for (p = 0; p < np; p++) {
        pimpl->ratio->residual += pimpl->ratio->t1[ p ];
    }

    /* Jacobian row */

    vector_zero(1 + (G->cell_facepos[c + 1] - G->cell_facepos[c]),
                pimpl->ratio->mat_row);

    /* t2 <- A \ ((dA/dp) * t1) */
    matvec(np, np, dAc, pimpl->ratio->t1, pimpl->ratio->t2);
    solve_linear_systems(np, 1, pimpl->ratio, pimpl->ratio->t2);

    dF2 = 0.0;
    for (p = 0; p < np; p++) {
        dF2 += pimpl->ratio->t2[ p ];
    }

    pimpl->is_incomp           = pimpl->is_incomp && (! (fabs(dF2) > 0));
    pimpl->ratio->mat_row[ 0 ] = - dF2;

    /* Accumulate inter-cell Jacobian contributions */
    dv  = pimpl->ratio->linsolve_buffer + (1 + nconn)*np;
    off = 1;
    for (i = G->cell_facepos[c]; i < G->cell_facepos[c + 1]; i++, off++) {

        f  = G->cell_faces[i];
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if ((c1 >= 0) && (c2 >= 0)) {
            if (c1 == c) { s =  1.0; dv1 = dv + 0 ; dv2 = dv + np; }
            else         { s = -1.0; dv1 = dv + np; dv2 = dv + 0 ; }

            dF1 = dF2 = 0.0;
            for (p = 0; p < np; p++) {
                dF1 += dv1[ p ];
                dF2 += dv2[ p ];
            }

            pimpl->ratio->mat_row[  0  ] += s * dt * dF1;
            pimpl->ratio->mat_row[ off ] += s * dt * dF2;

            dv += 2 * np;       /* '2' == number of one-sided derivatives. */
        }
    }
}


static void
assemble_sources(double                    dt ,
                 struct compr_src         *src,
                 struct cfs_tpfa_res_data *h  )
{
    int i;

    for (i = 0; i < src->nsrc; i++) {
        assert (src->cell[i]            >= 0      );
        assert (((size_t) src->cell[i]) <  h->J->m);

        h->F[ src->cell[ i ] ] -= dt * src->flux[ i ];
    }
}


/* ---------------------------------------------------------------------- */
static int
assemble_cell_contrib(struct UnstructuredGrid  *G,
                      int                       c,
                      struct cfs_tpfa_res_data *h)
/* ---------------------------------------------------------------------- */
{
    int c1, c2, i, f, j1, j2, off;

    j1 = csrmatrix_elm_index(c, c, h->J);

    h->J->sa[j1] += h->pimpl->ratio->mat_row[ 0 ];

    off = 1;
    for (i = G->cell_facepos[c]; i < G->cell_facepos[c + 1]; i++, off++) {
        f = G->cell_faces[i];

        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        c2 = (c1 == c) ? c2 : c1;

        if (c2 >= 0) {
            j2 = csrmatrix_elm_index(c, c2, h->J);

            h->J->sa[j2] += h->pimpl->ratio->mat_row[ off ];
        }
    }

    h->F[ c ] = h->pimpl->ratio->residual;

    return 0;
}


static void
init_completion_contrib(int                       i    ,
                        int                       np   ,
                        const double             *Ac   ,
                        const double             *dAc  ,
                        struct cfs_tpfa_res_impl *pimpl)
{
    memcpy(pimpl->ratio->linsolve_buffer,
           pimpl->compflux_p + (i * np),
           np * sizeof *pimpl->ratio->linsolve_buffer);

    memcpy(pimpl->ratio->linsolve_buffer + (1 * np),
           pimpl->compflux_deriv_p + (i * 2 * np),
           2 * np * sizeof *pimpl->ratio->linsolve_buffer);

    /* buffer <- Ac \ [A_{wi}q_{wi}, A_{wi} dq_{wi}] */
    factorise_fluid_matrix(np, Ac, pimpl->ratio);
    solve_linear_systems  (np, 1 + 2, pimpl->ratio,
                           pimpl->ratio->linsolve_buffer);

    /* t1 <- Ac \ (A_{wi} q_{wi}) */
    memcpy(pimpl->ratio->t1,
           pimpl->ratio->linsolve_buffer,
           np * sizeof *pimpl->ratio->t1);

    /* t2 <- Ac \ ((dA/dp) * t1) (== -d(Ac^{-1})/dp (A_{wi} q_{wi})) */
    matvec(np, np, dAc, pimpl->ratio->t1, pimpl->ratio->t2);
    solve_linear_systems(np, 1, pimpl->ratio, pimpl->ratio->t2);
}


static void
assemble_completion_to_cell(int c, int wdof, int np, double dt,
                            struct cfs_tpfa_res_data *h)
{
    int    p;
    size_t jc, jw;
    double s1, s2, *d1, *d2;

    /* Accumulate residual contributions and (dA^{-1}/dp) terms as
     * sums of phase contributions. */
    for (p = 0, s1 = 0.0, s2 = 0.0; p < np; p++) {
        s1 += h->pimpl->ratio->t1[ p ];
        s2 += h->pimpl->ratio->t2[ p ];
    }

    /* Assemble residual contributions from well completion.
     *
     * Note negative sign due to perforation flux convention (positive
     * flux into reservoir). */
    h->F[ c ] -= dt * s1;

    /* Assemble Jacobian contributions from well completion. */
    assert (wdof > c);
    jc = csrmatrix_elm_index(c, c   , h->J);
    jw = csrmatrix_elm_index(c, wdof, h->J);

    /* Compressibility-like (diagonal) Jacobian term.  Positive sign
     * since the negative derivative in ->ratio->t2 (see
     * init_completion_contrib()) is cancelled by completion flux
     * convention. */
    h->J->sa[ jc ] += dt * s2;

    /* Linear terms arising from simple differentiation of reservoir
     * volume flux on completion. */
    d1 = h->pimpl->ratio->linsolve_buffer + (1 * np);
    d2 = d1                               + (1 * np);
    for (p = 0, s1 = 0.0, s2 = 0.0; p < np; p++) {
        s1 += d1[ p ];
        s2 += d2[ p ];
    }

    h->J->sa[ jc ] += dt * s1;
    h->J->sa[ jw ] += dt * s2;
}


/* ---------------------------------------------------------------------- */
static void
welleq_coeff_shut(int np, struct cfs_tpfa_res_data *h,
                  double *res, double *w2c, double *w2w)
/* ---------------------------------------------------------------------- */
{
    int           p;
    double        fwi;
    const double *dpflux_w;

    /* Sum reservoir phase flux derivatives set by
     * compute_darcyflux_and_deriv(). */
    dpflux_w = h->pimpl->flux_work + (1 * np);
    for (p = 0, fwi = 0.0; p < np; p++) {
        fwi += dpflux_w[ p ];
    }

    *res = 0.0;
    *w2c = 0.0;
    *w2w = fwi;
}


/* ---------------------------------------------------------------------- */
static void
welleq_coeff_bhp(int np, double dp, struct cfs_tpfa_res_data *h,
                 double *res, double *w2c, double *w2w)
/* ---------------------------------------------------------------------- */
{
    int           p;
    double        fwi;
    const double *dpflux_w;

    /* Sum reservoir phase flux derivatives set by
     * compute_darcyflux_and_deriv(). */
    dpflux_w = h->pimpl->flux_work + (1 * np);
    for (p = 0, fwi = 0.0; p < np; p++) {
        fwi += dpflux_w[ p ];
    }

    *res = fwi * dp;
    *w2c = 0.0;
    *w2w = fwi;
}


/* ---------------------------------------------------------------------- */
static void
welleq_coeff_resv(int np, struct cfs_tpfa_res_data *h,
                  struct WellControls *ctrl,
                  double *res, double *w2c, double *w2w)
/* ---------------------------------------------------------------------- */
{
    int           p;
    double        f;
    const double *pflux, *dpflux_w, *dpflux_c, *distr;

     /* Sum reservoir phase flux and its derivatives set by
      * compute_darcyflux_and_deriv(). */
    pflux    = h->pimpl->flux_work;
    dpflux_w = pflux       + (1             * np);
    dpflux_c = dpflux_w    + (1             * np);
    distr    = ctrl->distr + (ctrl->current * np);

    *res = *w2c = *w2w = 0.0;
    for (p = 0; p < np; p++) {
        f = distr[ p ];

        *res += f * pflux   [ p ];
        *w2w += f * dpflux_w[ p ];
        *w2c += f * dpflux_c[ p ];
    }
}


/* ---------------------------------------------------------------------- */
static void
welleq_coeff_surfrate(int i, int np, struct cfs_tpfa_res_data *h,
                      struct WellControls *ctrl,
                      double *res, double *w2c, double *w2w)
/* ---------------------------------------------------------------------- */
{
    int           p;
    double        f;
    const double *pflux, *dpflux_w, *dpflux_c, *distr;

    pflux    = h->pimpl->compflux_p       + (i             * (1 * np));
    dpflux_w = h->pimpl->compflux_deriv_p + (i             * (2 * np));
    dpflux_c = dpflux_w                   + (1             * (1 * np));
    distr    = ctrl->distr                + (ctrl->current * (1 * np));

    *res = *w2c = *w2w = 0.0;
    for (p = 0; p < np; p++) {
        f = distr[ p ];

        *res += f * pflux   [ p ];
        *w2w += f * dpflux_w[ p ];
        *w2c += f * dpflux_c[ p ];
    }
}


/* ---------------------------------------------------------------------- */
static void
assemble_completion_to_well(int i, int w, int c, int nc, int np,
                            double pw, double dt,
                            struct cfs_tpfa_res_wells *wells,
                            struct cfs_tpfa_res_data  *h    )
/* ---------------------------------------------------------------------- */
{
    int    wdof;
    size_t jc, jw;
    double res = 0.0, w2c = 0.0, w2w = 0.0;

    struct Wells        *W;
    struct WellControls *ctrl;

    assert (wells    != NULL);
    assert (wells->W != NULL);

    W    = wells->W;
    ctrl = W->ctrls[ w ];

    if (ctrl->current < 0) {
        /* Interpreting a negative current control index to mean a shut well */
        welleq_coeff_shut(np, h, &res, &w2c, &w2w);
    }
    else {
        switch (ctrl->type[ ctrl->current ]) {
        case BHP :
            welleq_coeff_bhp(np, pw - ctrl->target[ ctrl->current ],
                             h, &res, &w2c, &w2w);
            break;

        case RESERVOIR_RATE:
            assert (W->number_of_phases == np);
            welleq_coeff_resv(np, h, ctrl, &res, &w2c, &w2w);
            break;

        case SURFACE_RATE:
            welleq_coeff_surfrate(i, np, h, ctrl, &res, &w2c, &w2w);
            break;
        }
    }

    /* Assemble completion contributions */
    wdof = nc + w;
    jc   = csrmatrix_elm_index(wdof, c   , h->J);
    jw   = csrmatrix_elm_index(wdof, wdof, h->J);

    h->F    [ wdof ] += dt * res;
    h->J->sa[ jc   ] += dt * w2c;
    h->J->sa[ jw   ] += dt * w2w;
}


static int
assemble_well_contrib(struct cfs_tpfa_res_wells   *wells ,
                      struct compr_quantities_gen *cq    ,
                      double                       dt    ,
                      const double                *cpress,
                      const double                *wpress,
                      struct cfs_tpfa_res_data    *h     )
{
    int           w, i, c, np, np2, nc;
    int           is_neumann, is_open;
    double        pw, dp, gpot[3] = { 0.0 };
    const double *WI, *wdp, *pmobp;
    const double *Ac, *dAc;

    struct Wells        *W;
    struct WellControls *ctrl;

    nc  = ((int) h->J->m) - wells->W->number_of_wells;
    np  = cq->nphases;
    np2 = np * np;

    W = wells->W;

    WI    = W->WI;
    wdp   = wells->data->wdp;
    pmobp = wells->data->phasemob;

    is_neumann = 1;

    for (w = i = 0; w < W->number_of_wells; w++) {
        pw = wpress[ w ];
        is_open = W->ctrls[w]->current >= 0;

        for (; i < W->well_connpos[w + 1]; i++, pmobp += np) {

            c   = W->well_cells[ i ];
            Ac  = cq->Ac  + (c * np2);
            dAc = cq->dAc + (c * np2);

            dp  = pw + wdp[i] - cpress[ c ];

            init_completion_contrib(i, np, Ac, dAc, h->pimpl);

            if (is_open) {
                assemble_completion_to_cell(c, nc + w, np, dt, h);
            }

            /* Prepare for RESV controls */
            compute_darcyflux_and_deriv(np, WI[i], dp, pmobp, gpot,
                                        h->pimpl->flux_work,
                                        h->pimpl->flux_work + np);

            assemble_completion_to_well(i, w, c, nc, np, pw, dt, wells, h);
        }

        ctrl = W->ctrls[ w ];
        if ((ctrl->current >= 0) && /* OPEN? */
            (ctrl->type[ ctrl->current ] != BHP)) {

            h->F[ nc + w ] -= dt * ctrl->target[ ctrl->current ];
        }
        else {
            is_neumann = 0;
        }
    }

    return is_neumann;
}


/* ---------------------------------------------------------------------- */
static void
compute_fpress(struct UnstructuredGrid *G        ,
               int                      np       ,
               const double            *htrans   ,
               const double            *pmobf    ,
               const double            *gravcap_f,
               const double            *cpress   ,
               const double            *fflux    ,
               double                  *fpress   ,
               double                  *scratch_f)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f; /* , c1, c2; */

    /* Suppress warning about unused parameters. */
    (void) np;  (void) pmobf;  (void) gravcap_f;  (void) fflux;

    /*
     * Define face pressures as weighted average of connecting cell
     * pressures.  Specifically, we define
     *
     *     pf = (t1 p1 + t2 p2) / (t1 + t2)
     *
     * in which t1 and t2 are the one-sided transmissibilities and p1
     * and p2 are the associated cell pressures.
     *
     * NOTE: The formula does not account for effects of gravity or
     * flux boundary conditions.
     */
    for (f = 0; f < G->number_of_faces; f++) {
        scratch_f[f] = fpress[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            scratch_f[f] += htrans[i];
            fpress[f]    += htrans[i] * cpress[c];
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        fpress[f] /= scratch_f[f];
#if 0
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if (((c1 < 0) || (c2 < 0)) &&
            (bc != NULL) && (bc->type[f] == PRESSURE)) {
            fpress[f] = bc->bcval[f];
        }
#endif
    }
}


/* ---------------------------------------------------------------------- */
static void
compute_flux(struct UnstructuredGrid *G        ,
             int                      np       ,
             const double            *trans    ,
             const double            *pmobf    ,
             const double            *gravcap_f,
             const double            *cpress   ,
             double                  *fflux    )
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2, p;
    double t, dp, g;

    for (f = 0; f < G->number_of_faces; f++) {
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

#if 0
        if (((c1 < 0) || (c2 < 0)) &&
            (bc != NULL) && (bc->type[f] == FLUX)) {
            fflux[f] = bc->bcval[f];
            continue;
        }
#endif

        t = g = 0.0;
        for (p = 0; p < np; p++) {
            t += pmobf[f*np + p];
            g += pmobf[f*np + p] * gravcap_f[f*np + p];
        }
        /* t *= trans[f]; */

        if ((c1 >= 0) && (c2 >= 0)) {
            dp  = cpress[c1] - cpress[c2];
#if 0
        } else if ((bc != NULL) && (bc->type[f] == PRESSURE)) {
            if (c1 < 0) {
                dp = bc->bcval[f] - cpress[c2];
                /* g  = -g; */
            } else {
                dp = cpress[c1] - bc->bcval[f];
            }
#endif
        } else {
            /* No BC -> no-flow (== pressure drop offsets gravity) */
            dp = -g / t;
        }

        fflux[f] = trans[f] * (t*dp + g);
    }
}


static void
compute_wflux(int                        np    ,
              struct cfs_tpfa_res_wells *wells ,
              const double              *pmobc ,
              const double              *cpress,
              const double              *wpress,
              double                    *wflux )
{
    int           w, c, i, p;
    double        pw, dp, t;
    const double *pmob, *wdp;

    struct Wells          *W;
    struct CompletionData *cdata;

    assert (wells       != NULL);
    assert (wells->W    != NULL);
    assert (wells->data != NULL);

    W     = wells->W;
    cdata = wells->data;
    wdp   = cdata->wdp;

    for (w = i = 0; w < W->number_of_wells; w++) {
        pw = wpress[w];

        for (; i < W->well_connpos[w + 1]; i++) {
            c  = W->well_cells[ i ];
            dp = pw + wdp[ i ] - cpress[c];

            if (dp > 0) { pmob = cdata->phasemob + (i * np); } /* w->c */
            else        { pmob = pmobc           + (c * np); } /* c->w */

            for (p = 0, t = 0.0; p < np; p++) {
                t += pmob[ p ];
            }

            wflux[ i ] = W->WI[i] * t * dp;
        }
    }
}


static size_t
maxconn(struct UnstructuredGrid *G)
{
    int    c;
    size_t m, n;

    for (c = 0, m = 0; c < G->number_of_cells; c++) {
        n = G->cell_facepos[c + 1] - G->cell_facepos[c];

        if (n > m) { m = n; }
    }

    return m;
}


/* ======================================================================
 * Public interface below separator.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_destroy(struct cfs_tpfa_res_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        csrmatrix_delete(h->J);
        impl_deallocate (h->pimpl);
    }

    free(h);
}


/* ---------------------------------------------------------------------- */
struct cfs_tpfa_res_data *
cfs_tpfa_res_construct(struct UnstructuredGrid   *G      ,
                       struct cfs_tpfa_res_wells *wells  ,
                       int                        nphases)
/* ---------------------------------------------------------------------- */
{
    size_t                    nf, nwperf;
    struct cfs_tpfa_res_data *h;

    h = malloc(1 * sizeof *h);

    if (h != NULL) {
        h->pimpl = impl_allocate(G, wells, maxconn(G), nphases);
        h->J     = construct_matrix(G, wells);

        if ((h->pimpl == NULL) || (h->J == NULL)) {
            cfs_tpfa_res_destroy(h);
            h = NULL;
        }
    }

    if (h != NULL) {
        nf     = G->number_of_faces;
        nwperf = 0;

        if ((wells != NULL) && (wells->W != NULL)) {
            nwperf = wells->W->well_connpos[ wells->W->number_of_wells ];
        }

        /* Allocate linear system components */
        h->F                       = h->pimpl->ddata + 0;

        h->pimpl->compflux_f       = h->F            + h->J->m;
        h->pimpl->compflux_p       =
            h->pimpl->compflux_f                     + (nphases * nf);

        h->pimpl->compflux_deriv_f =
            h->pimpl->compflux_p                     + (nphases * nwperf);
        h->pimpl->compflux_deriv_p =
            h->pimpl->compflux_deriv_f               + (nphases * 2 * nf);

        h->pimpl->flux_work        =
            h->pimpl->compflux_deriv_p               + (nphases * 2 * nwperf);

        h->pimpl->scratch_f        =
            h->pimpl->flux_work                      + (nphases * (1 + 2));
    }

    return h;
}


/* ---------------------------------------------------------------------- */
int
cfs_tpfa_res_assemble(struct UnstructuredGrid     *G        ,
                      double                       dt       ,
                      struct cfs_tpfa_res_forces  *forces   ,
                      const double                *zc       ,
                      struct compr_quantities_gen *cq       ,
                      const double                *trans    ,
                      const double                *gravcap_f,
                      const double                *cpress   ,
                      const double                *wpress   ,
                      const double                *porevol  ,
                      struct cfs_tpfa_res_data    *h        )
/* ---------------------------------------------------------------------- */
{
    int res_is_neumann, well_is_neumann, c, np2, singular;

    csrmatrix_zero(         h->J);
    vector_zero   (h->J->m, h->F);

    h->pimpl->is_incomp = 1;

    compute_compflux_and_deriv(G, cq->nphases, cpress, trans,
                               cq->phasemobf, gravcap_f, cq->Af, h->pimpl);

    res_is_neumann  = 1;
    well_is_neumann = 1;

    np2 = cq->nphases * cq->nphases;
    for (c = 0; c < G->number_of_cells;
         c++, zc += cq->nphases) {

        compute_cell_contrib(G, c, cq->nphases, porevol[c], dt, zc,
                             cq->Ac + (c * np2), cq->dAc + (c * np2),
                             h->pimpl);

        assemble_cell_contrib(G, c, h);
    }

    if ((forces           != NULL) &&
        (forces->wells    != NULL) &&
        (forces->wells->W != NULL)) {
        compute_well_compflux_and_deriv(forces->wells, cq->nphases,
                                        cpress, wpress, h->pimpl);

        well_is_neumann = assemble_well_contrib(forces->wells, cq, dt,
                                                cpress, wpress, h);
    }

    if ((forces != NULL) && (forces->src != NULL)) {
        assert (forces->src->nphases == cq->nphases);
        assemble_sources(dt, forces->src, h);
    }

    singular = res_is_neumann && well_is_neumann && h->pimpl->is_incomp;
    if (singular) {
        h->J->sa[0] *= 2.0;
    }

    return singular;
}


/* ---------------------------------------------------------------------- */
int
cfs_tpfa_res_comprock_assemble(
                      struct UnstructuredGrid     *G        ,
                      double                       dt       ,
                      struct cfs_tpfa_res_forces  *forces   ,
                      const double                *zc       ,
                      struct compr_quantities_gen *cq       ,
                      const double                *trans    ,
                      const double                *gravcap_f,
                      const double                *cpress   ,
                      const double                *wpress   ,
                      const double                *porevol  ,
                      const double                *porevol0 ,
                      const double                *rock_comp,
                      struct cfs_tpfa_res_data    *h        )
/* ---------------------------------------------------------------------- */
{
    /* We want to add this term to the usual residual:
     *
     * (porevol(pressure)-porevol(initial_pressure))/dt.
     *
     * Its derivative (for the diagonal term of the Jacobian) is:
     *
     * porevol(pressure)*rock_comp(pressure)/dt
     */

    int     c, rock_is_incomp, singular;
    size_t  j;
    double  dpv;

    /* Assemble usual system (without rock compressibility). */
    singular = cfs_tpfa_res_assemble(G, dt, forces, zc, cq, trans, gravcap_f,
                                     cpress, wpress, porevol0, h);

    /* If we made a singularity-removing adjustment in the
       regular assembly, we undo it here. */
    if (singular) {
        h->J->sa[0] /= 2.0;
    }

    /* Add new terms to residual and Jacobian. */
    rock_is_incomp = 1;
    for (c = 0; c < G->number_of_cells; c++) {
        j = csrmatrix_elm_index(c, c, h->J);

        dpv = (porevol[c] - porevol0[c]);
        if (dpv != 0.0 || rock_comp[c] != 0.0) {
            rock_is_incomp = 0;
        }

        h->J->sa[j] += porevol[c] * rock_comp[c];
        h->F[c]     += dpv;
    }

    /* Re-do the singularity-removing adjustment if necessary */
    if (rock_is_incomp && singular) {
        h->J->sa[0] *= 2.0;
    }

    return rock_is_incomp && singular;
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_flux(struct UnstructuredGrid    *G        ,
                  struct cfs_tpfa_res_forces *forces   ,
                  int                         np       ,
                  const double               *trans    ,
                  const double               *pmobc    ,
                  const double               *pmobf    ,
                  const double               *gravcap_f,
                  const double               *cpress   ,
                  const double               *wpress   ,
                  double                     *fflux    ,
                  double                     *wflux    )
/* ---------------------------------------------------------------------- */
{
    compute_flux(G, np, trans, pmobf, gravcap_f, cpress, fflux);

    if ((forces           != NULL) &&
        (forces->wells    != NULL) &&
        (forces->wells->W != NULL) && (wpress != NULL) && (wflux != NULL)) {

        compute_wflux(np, forces->wells, pmobc, cpress, wpress, wflux);
    }
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_fpress(struct UnstructuredGrid  *G        ,
                    int                       np       ,
                    const double             *htrans   ,
                    const double             *pmobf    ,
                    const double             *gravcap_f,
                    struct cfs_tpfa_res_data *h        ,
                    const double             *cpress   ,
                    const double             *fflux    ,
                    double                   *fpress   )
/* ---------------------------------------------------------------------- */
{
    compute_fpress(G, np, htrans, pmobf, gravcap_f,
                   cpress, fflux, fpress, h->pimpl->scratch_f);
}


#if 0
/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_retrieve_masstrans(struct UnstructuredGrid               *G,
                            int                   np,
                            struct cfs_tpfa_res_data *h,
                            double               *masstrans_f)
/* ---------------------------------------------------------------------- */
{
    memcpy(masstrans_f, h->pimpl->masstrans_f,
           np * G->number_of_faces * sizeof *masstrans_f);
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_retrieve_gravtrans(struct UnstructuredGrid               *G,
                            int                   np,
                            struct cfs_tpfa_res_data *h,
                            double               *gravtrans_f)
/* ---------------------------------------------------------------------- */
{
    memcpy(gravtrans_f, h->pimpl->gravtrans_f,
           np * G->number_of_faces * sizeof *gravtrans_f);
}


/* ---------------------------------------------------------------------- */
static double
cfs_tpfa_res_impes_maxtime_cell(int                      c,
                            struct UnstructuredGrid                  *G,
                            struct compr_quantities *cq,
                            const double            *trans,
                            const double            *porevol,
                            struct cfs_tpfa_res_data    *h,
                            const double            *dpmobf,
                            const double            *surf_dens,
                            const double            *gravity)
/* ---------------------------------------------------------------------- */
{
    /* Reference:
       K. H. Coats, "IMPES Stability: The Stable Step", SPE 69225

       Capillary pressure parts not included.
    */

    int i, j, k, f, c2;
    double f11, f12, f21, f22;
    double dp, dzg, tr, tmob, detF, eqv_flux;
    const double *pmob;
    const double *A;
    /* This is intended to be compatible with the dune-porsol blackoil
       code. This library is otherwise not depending on particular
       orderings of phases or components, so at some point this
       function should be generalized. */
    enum { Water = 0, Oil = 1, Gas = 2 };
    enum { num_phases = 3 };
    double rho[num_phases];
    double pot[num_phases];
    /* Notation: dpmob[Oil][Water] is d/ds_w(lambda_o) */
    double dpmob[num_phases][num_phases]
        = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };

    assert (cq->nphases == 3);

    f11 = f12 = f21 = f22 = 0.0;

    /* Loop over neighbour faces to accumulate f11, f12 etc. */
    for (i = G->cell_facepos[c]; i < G->cell_facepos[c + 1]; ++i) {
        f  = G->cell_faces[i];
        if ((c2 = G->face_cells[2*f + 0]) == c) {
            c2  = G->face_cells[2*f + 1];
        }

        /* Initially only interior faces */
        if (c2 < 0) {
            continue;
        }

        /* Computing density */
        A = cq->Af + f*(cq->nphases)*(cq->nphases);
        for (j = 0; j < cq->nphases; ++j) {
            rho[j] = 0.0;
            for (k = 0; k < cq->nphases; ++k) {
                rho[j] += A[cq->nphases*j + k]*surf_dens[k];
            }
        }
        /* Computing gravity potentials */
        dp = h->x[c] - h->x[c2];
        dzg = 0.0;
        for (j = 0; j < G->dimensions; ++j) {
            dzg += (G->cell_centroids[G->dimensions*c + j] - G->cell_centroids[G->dimensions*c2 + j])*gravity[j];
        }
        for (j = 0; j < cq->nphases; ++j) {
            pot[j] = fabs(dp - rho[j]*dzg);
        }
        /* Filling the dpmob array from available data.
           Note that we only need the following combinations:
           (Water, Water)
           (Water, Gas)
           (Oil, Water)
           (Oil, Gas)
           (Gas, Gas)

           No derivatives w.r.t. Oil is needed, since there are only two
           independent saturation variables.

           The lack of (Gas, Water) may be due to assumptions on the
           three-phase model used (should be checked to be compatible
           with our choices).
        */
        dpmob[Water][Water] = dpmobf[9*f];
        dpmob[Water][Gas] = dpmobf[9*f + 2];
        dpmob[Oil][Water] = dpmobf[9*f + 3];
        dpmob[Oil][Gas] = dpmobf[9*f + 5];
        dpmob[Gas][Gas] = dpmobf[9*f + 8];
        /* Computing the flux parts f_ij */
        pmob = cq->phasemobf + f*cq->nphases;
        tr = trans[f];
        tmob = pmob[Water] + pmob[Oil] + pmob[Gas];
        f11 += tr*((pmob[Oil] + pmob[Gas])*dpmob[Water][Water]*pot[Water]
                   - pmob[Water]*dpmob[Oil][Water]*pot[Oil])/tmob;
        f12 += -tr*(pmob[Water]*dpmob[Oil][Gas]*pot[Oil]
                    + pmob[Water]*dpmob[Gas][Gas]*pot[Gas]
                    - (pmob[Oil] + pmob[Gas])*dpmob[Water][Gas]*pot[Water])/tmob;
        f21 += -tr*(pmob[Gas]*dpmob[Water][Water]*pot[Water]
                    + pmob[Gas]*dpmob[Oil][Water]*pot[Oil])/tmob;
        f22 += tr*(-pmob[Gas]*dpmob[Oil][Gas]*pot[Oil]
                   + (pmob[Water] + pmob[Oil])*dpmob[Gas][Gas]*pot[Gas]
                   - pmob[Gas]*dpmob[Water][Gas]*pot[Water])/tmob;
    }

    /* (from eq. 3, 4a-e, 5a-c)
       F_i = 1/2 |f11_i + f22_i + \sqrt{G}|
       G = (f11_i + f22_i)^2 - 4 det(F_i)
       fXX_i = \sum_j fXX_ij     (j runs over the neighbours)
       det(F_i) = f11_i f22_i - f12_i f21_i
    */
    detF = f11*f22 - f12*f21;
    eqv_flux = 0.5*fabs(f11 + f22 + sqrt((f11 + f22)*(f11 + f22) - 4*detF));
    /* delta_t < porevol/eqv_flux */
    return porevol[c]/eqv_flux;
}

/* ---------------------------------------------------------------------- */
double
cfs_tpfa_res_impes_maxtime(struct UnstructuredGrid                  *G,
                       struct compr_quantities *cq,
                       const double            *trans,
                       const double            *porevol,
                       struct cfs_tpfa_res_data    *h,
                       const double            *dpmobf,
                       const double            *surf_dens,
                       const double            *gravity)
/* ---------------------------------------------------------------------- */
{
    int c;
    double max_dt, cell_dt;
    max_dt = 1e100;
    for (c = 0; c < G->number_of_cells; ++c) {
        cell_dt = cfs_tpfa_res_impes_maxtime_cell(c, G, cq, trans, porevol, h,
                                              dpmobf, surf_dens, gravity);
        if (cell_dt < max_dt) {
            max_dt = cell_dt;
        }
    }
    return max_dt;
}



/* ---------------------------------------------------------------------- */
void
cfs_tpfa_res_expl_mass_transport(struct UnstructuredGrid                 *G,
                             well_t                 *W,
                             struct completion_data *wdata,
                             int                     np,
                             double                  dt,
                             const double           *porevol,
                             struct cfs_tpfa_res_data   *h,
                             double                 *surf_vol)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f, c2, p, w;
    double dp, dz, gsgn;
    const double *masstrans_f, *gravtrans_f, *masstrans_p, *gravtrans_p;
    const double *cpress, *wpress;

    /* Suppress warning about unused parameter. */
    (void) wdata;

    masstrans_f = h->pimpl->masstrans_f;
    gravtrans_f = h->pimpl->gravtrans_f;
    masstrans_p = h->pimpl->masstrans_p;
    gravtrans_p = h->pimpl->gravtrans_p;
    cpress = h->x;
    wpress = h->x + G->number_of_cells;

    /* Transport through interior faces */
    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f  = G->cell_faces[i];

            if ((c2 = G->face_cells[2*f + 0]) == c) {
              gsgn = 1.0;
                c2  = G->face_cells[2*f + 1];
            } else {
              gsgn = -1.0;
            }

            if (c2 >= 0) {
                dp = cpress[c] - cpress[c2];

                for (p = 0; p < np; p++) {
                    dz  = masstrans_f[f*np + p] * dp;
                    dz += gravtrans_f[f*np + p] * gsgn;

                    /* dz > 0 => flow from c into c2. */
                    surf_vol[c*np + p] -= dz * dt / porevol[c];
                }
            }
        }
    }

    /* Transport through well perforations */
    if (W != NULL) {
        for (w = i = 0; w < W->number_of_wells; w++) {
            for (; i < W->well_connpos[w + 1]; i++) {
                c = W->well_cells[i];
                dp = wpress[w] - cpress[c];

                for (p = 0; p < np; p++) {
                    dz = masstrans_p[i*np + p] * dp;
                    dz += gravtrans_p[i*np + p];

                    /* dz > 0 => flow from perforation into c. */
                    surf_vol[c*np + p] += dz * dt / porevol[c];
                }
            }
        }
    }
}
#endif
