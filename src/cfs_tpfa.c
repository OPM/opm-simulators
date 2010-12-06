#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "flow_bc.h"
#include "well.h"

#include "compr_quant.h"
#include "trans_tpfa.h"
#include "cfs_tpfa.h"
#include "sparse_sys.h"


struct disc_data {
    double *ctrans, *P, *Xf, *Yf, *work;
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
deallocate_disc_data(struct disc_data *data)
/* ---------------------------------------------------------------------- */
{
    if (data != NULL) {
        free(data->ddata);
    }

    free(data);
}


/* ---------------------------------------------------------------------- */
static struct disc_data *
allocate_disc_data(grid_t *g, int np)
/* ---------------------------------------------------------------------- */
{
    size_t nc, nf, ngconn, ddata_sz;
    struct disc_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        nc      = g->number_of_cells;
        nf      = g->number_of_faces;
        ngconn  = g->cell_facepos[nc];

        ddata_sz  = ngconn;      /* ctrans */
        ddata_sz += nc;          /* P */
        ddata_sz += np * nf;     /* Xf */
        ddata_sz += np * ngconn; /* Yf */
        ddata_sz += ngconn;      /* work */

        new->ddata = malloc(ddata_sz * sizeof *new->ddata);

        if (new->ddata == NULL) {
            deallocate_disc_data(new);
            new = NULL;
        } else {
            new->ctrans = new->ddata  + 0      ;
            new->P      = new->ctrans + ngconn ;
            new->Xf     = new->P      + nc     ;
            new->Yf     = new->Xf     + np * nf;
            new->work   = new->Yf     + np * ngconn;
        }
    }

    return new;
}


struct cfs_tpfa_impl {
    double           *fpress;   /* Face pressure */
    double           *accum;

    /* One entry per component per face */
    double           *masstrans_f; /* RB^{-1} [ phase-mobility ] */
    double           *gravtrans_f; /* RB^{-1} [ grav + capillary ] */

    struct disc_data *dd;

    /* Linear storage */
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
impl_deallocate(struct cfs_tpfa_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        free                (pimpl->ddata);
        deallocate_disc_data(pimpl->dd);
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
static struct cfs_tpfa_impl *
impl_allocate(grid_t *G, well_t *W, int np)
/* ---------------------------------------------------------------------- */
{
    size_t                nnu;
    struct cfs_tpfa_impl *new;

    size_t ddata_sz;

    nnu = G->number_of_cells;
    if (W != NULL) {
        nnu += W->number_of_wells;
    }

    ddata_sz  = 2  * nnu;                /* b, x */

    ddata_sz += 1  * G->number_of_faces; /* fpress */
    ddata_sz += 1  * G->number_of_faces; /* accum */
    ddata_sz += np * G->number_of_faces; /* masstrans_f */
    ddata_sz += np * G->number_of_faces; /* gravtrans_f */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->ddata = malloc(ddata_sz * sizeof *new->ddata);
        new->dd    = allocate_disc_data(G, np);

        if (new->ddata == NULL || new->dd == NULL) {
            impl_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct CSRMatrix *
construct_matrix(grid_t *G, well_t *W)
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

        /* Enforce sorted connection structure per row */
        csrmatrix_sortrows(A);
    }

    return A;
}


/* ---------------------------------------------------------------------- */
static void
solve_cellsys_core(grid_t       *G   ,
                   size_t        sz  ,
                   const double *Ac  ,
                   const double *bf  ,
                   double       *xcf ,
                   double       *luAc,
                   MAT_SIZE_T   *ipiv)
/* ---------------------------------------------------------------------- */
{
    int         c, i, f;
    size_t      j, p2;
    double     *v;

    MAT_SIZE_T  nrows, ncols, ldA, ldX, nrhs, info;

    nrows = ncols = ldA = ldX = sz;
    info  = 0;

    v     = xcf;

    for (c = 0, p2 = 0; c < G->number_of_cells; c++) {
        /* Define right-hand sides for local systems */
        for (i = G->cell_facepos[c + 0], nrhs = 0;
             i < G->cell_facepos[c + 1]; i++, nrhs++) {
            f = G->cell_faces[i];

            for (j = 0; j < sz; j++) {
                v[nrhs*sz + j] = bf[f*sz + j];
            }
        }

        /* Factor Ac */
        memcpy(luAc, Ac + p2, sz * sz * sizeof *luAc);
        dgetrf_(&nrows, &ncols, luAc, &ldA, ipiv, &info);

        /* Solve local systems */
        dgetrs_("No Transpose", &nrows, &nrhs,
                luAc, &ldA, ipiv, v, &ldX, &info);

        v  += nrhs * sz;
        p2 += sz   * sz;
    }
}


/* ---------------------------------------------------------------------- */
static void
small_matvec(size_t n, int sz,
             const double *A,
             const double *X,
             double       *Y)
/* ---------------------------------------------------------------------- */
{
    size_t i, p1, p2;

    MAT_SIZE_T nrows, ncols, ld, incx, incy;
    double     a1, a2;

    nrows = ncols = ld = sz;
    incx  = incy  = 1;

    a1 = 1.0;
    a2 = 0.0;
    for (i = p1 = p2 = 0; i < n; i++) {
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, A + p2, &ld, X + p1, &incx,
               &a2,              Y + p1, &incy);

        p1 += sz;
        p2 += sz * sz;
    }
}


/* ---------------------------------------------------------------------- */
static int
solve_cellsys(grid_t       *G ,
              size_t        sz,
              const double *Ac,
              const double *bf,
              double       *xcf)
/* ---------------------------------------------------------------------- */
{
    int     ret;
    double *luAc;

    MAT_SIZE_T *ipiv;

    luAc = malloc(sz * sz * sizeof *luAc);
    ipiv = malloc(sz      * sizeof *ipiv);

    if ((luAc != NULL) && (ipiv != NULL)) {
        solve_cellsys_core(G, sz, Ac, bf, xcf, luAc, ipiv);
        ret = 1;
    } else {
        ret = 0;
    }

    free(ipiv);  free(luAc);

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
set_dynamic_trans(grid_t                  *G    ,
                  const double            *trans,
                  struct compr_quantities *cq   ,
                  struct disc_data        *dd)
/* ---------------------------------------------------------------------- */
{
    int f, p, i;

    for (f = i = 0; f < G->number_of_faces; f++) {
        for (p = 0; p < cq->nphases; p++, i++) {
            dd->Xf[i] = trans[f] * cq->phasemobf[i];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
sum_phase_contrib(grid_t       *G  ,
                  size_t        sz ,
                  const double *xcf,
                  double       *sum)
/* ---------------------------------------------------------------------- */
{
    int    c, i;
    size_t j;

    const double *v;

    for (c = i = 0, v = xcf; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++, v += sz) {

            sum[i] = 0.0;
            for (j = 0; j < sz; j++) {
                sum[i] += v[j];
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
set_dynamic_grav(grid_t                  *G        ,
                 const double            *trans    ,
                 const double            *gravcap_f,
                 struct compr_quantities *cq       ,
                 struct disc_data        *dd)
/* ---------------------------------------------------------------------- */
{
    int f, p, i;

    for (f = i = 0; f < G->number_of_faces; f++) {
        for (p = 0; p < cq->nphases; p++, i++) {
            dd->Xf[i] = trans[f] * gravcap_f[i] * cq->phasemobf[i];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
compute_densrat_update(grid_t                  *G ,
                       struct compr_quantities *cq,
                       struct disc_data        *dd,
                       double                  *q)
/* ---------------------------------------------------------------------- */
{
    small_matvec(G->number_of_faces, cq->nphases, cq->Af, dd->Xf, q);

    solve_cellsys(G, cq->nphases, cq->Ac, q, dd->Yf);

    sum_phase_contrib(G, cq->nphases, dd->Yf, dd->work);
}


/* ---------------------------------------------------------------------- */
static void
compute_psys_contrib(grid_t                  *G,
                     struct compr_quantities *cq,
                     double                   dt,
                     const double            *trans,
                     const double            *gravcap_f,
                     const double            *cpress0,
                     const double            *porevol,
                     struct cfs_tpfa_data    *h)
/* ---------------------------------------------------------------------- */
{
    int    c, nc, i;
    size_t nconn;

    nc    = G->number_of_cells;
    nconn = G->cell_facepos[nc];

    /* Compressible half-trans */
    set_dynamic_trans(G, trans, cq, h->pimpl->dd);
    compute_densrat_update(G, cq, h->pimpl->dd,
                           h->pimpl->masstrans_f);
    memcpy(h->pimpl->dd->ctrans,
           h->pimpl->dd->work,
           nconn * sizeof *h->pimpl->dd->ctrans);

    /* Compressible gravity contributions */
    set_dynamic_grav(G, trans, gravcap_f, cq, h->pimpl->dd);
    compute_densrat_update(G, cq, h->pimpl->dd,
                           h->pimpl->gravtrans_f);

    for (c = 0, i = 0; c < nc; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            h->b[c] -= h->pimpl->dd->work[i];
        }
        h->b[c] += cq->voldiscr[c];
    }

    /* Compressible accumulation term (lhs and rhs) */
    compr_accum_term(nc, dt, porevol, cq->totcompr, h->pimpl->dd->P);
    compr_src_add_press_accum(nc, cpress0, h->pimpl->dd->P, h->b);
}


/* ---------------------------------------------------------------------- */
static void
compute_fpress(grid_t       *G,
               flowbc_t     *bc,
               int           np,
               const double *htrans,
               const double *pmobf,
               const double *cpress,
               const double *fflux,
               double       *fpress)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f, p, c1, c2;
    double t, s;

    for (f = 0; f < G->number_of_faces; f++) {
        fpress[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            t = 0.0;
            for (p = 0; p < np; p++) {
                t += pmobf[f*np + p];
            }
            t *= htrans[i];

            s = 2.0*(G->face_cells[2*f + 0] == c) - 1.0;

            fpress[f] += cpress[c] - (s * fflux[f] / t);
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        fpress[f] /= (c1 >= 0) + (c2 >= 0);

        if (((c1 < 0) || (c2 < 0)) && (bc->type[f] == PRESSURE)) {
            fpress[f] = bc->bcval[f];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
compute_flux(grid_t       *G,
             flowbc_t     *bc,
             int           np,
             const double *trans,
             const double *pmobf,
             const double *cpress,
             double       *fflux)
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2, p;
    double t, dp;

    for (f = 0; f < G->number_of_faces; f++) {
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if (((c1 < 0) || (c2 < 0)) && (bc->type[f] == FLUX)) {
            fflux[f] = bc->bcval[f];
            continue;
        }

        t = 0.0;
        for (p = 0; p < np; p++) {
            t += pmobf[f*np + p];
        }
        t *= trans[f];

        if ((c1 >= 0) && (c2 >= 0)) {
            dp = cpress[c1] - cpress[c2];
        } else if (bc->type[f] == PRESSURE) {
            if (c1 < 0) {
                dp = bc->bcval[f] - cpress[c2];
            } else {
                dp = cpress[c1] - bc->bcval[f];
            }
        } else {
            /* No BC -> no-flow (== zero pressure drop) */
            dp = 0.0;
        }

        fflux[f] = t * dp;
    }
}



/* ======================================================================
 * Public interface below separator.
 * ====================================================================== */

/* ---------------------------------------------------------------------- */
struct cfs_tpfa_data *
cfs_tpfa_construct(grid_t *G, int nphases)
/* ---------------------------------------------------------------------- */
{
    size_t                nf;
    struct cfs_tpfa_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->pimpl = impl_allocate(G, NULL, nphases);
        new->A     = construct_matrix(G, NULL);

        if ((new->pimpl == NULL) || (new->A == NULL)) {
            cfs_tpfa_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        new->b = new->pimpl->ddata;
        new->x = new->b + new->A->m;

        nf = G->number_of_faces;

        new->pimpl->fpress      = new->x                  + new->A->m;
        new->pimpl->accum       = new->pimpl->fpress      + nf;
        new->pimpl->masstrans_f = new->pimpl->accum       + nf;
        new->pimpl->gravtrans_f = new->pimpl->masstrans_f + (nphases * nf);
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_assemble(grid_t                  *G,
                  double                   dt,
                  flowbc_t                *bc,
                  const double            *src,
                  struct compr_quantities *cq,
                  const double            *trans,
                  const double            *gravcap_f,
                  const double            *cpress0,
                  const double            *porevol,
                  struct cfs_tpfa_data    *h)
/* ---------------------------------------------------------------------- */
{
    int c1, c2, c, i, f, j1, j2;
    int is_neumann;

    double *ctrans = h->pimpl->dd->ctrans;

    csrmatrix_zero(         h->A);
    vector_zero   (h->A->m, h->b);

    compute_psys_contrib(G, cq, dt, trans, gravcap_f, cpress0, porevol, h);

    is_neumann = 1;

    for (c = i = 0; c < G->number_of_cells; c++) {
        j1 = csrmatrix_elm_index(c, c, h->A);

        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            c2 = (c1 == c) ? c2 : c1;

            if (c2 >= 0) {
                j2 = csrmatrix_elm_index(c, c2, h->A);

                h->A->sa[j1] += ctrans[i];
                h->A->sa[j2] -= ctrans[i];
            } else if (bc->type[f] == PRESSURE) {
                is_neumann    = 0;

                h->A->sa[j1] += ctrans[i];
                h->b    [c ] += ctrans[i] * bc->bcval[f];
            }
        }

        h->b[c] += src[c];

        /* Compressible accumulation term */
        h->A->sa[j1] += h->pimpl->dd->P[c];
    }

    if (is_neumann) {
        h->A->sa[0] *= 2;
    }
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_press_flux(grid_t               *G,
                    flowbc_t             *bc,
                    int                   np,
                    const double         *trans,
                    const double         *pmobf,
                    struct cfs_tpfa_data *h,
                    double               *cpress,
                    double               *fflux)
/* ---------------------------------------------------------------------- */
{
    /* Assign cell pressure directly from solution vector */
    memcpy(cpress, h->x, G->number_of_cells * sizeof *cpress);

    compute_flux(G, bc, np, trans , pmobf, cpress, fflux);
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_fpress(grid_t       *G,
                flowbc_t     *bc,
                int           np,
                const double *htrans,
                const double *pmobf,
                const double *cpress,
                const double *fflux,
                double       *fpress)
/* ---------------------------------------------------------------------- */
{
    compute_fpress(G, bc, np, htrans, pmobf, cpress, fflux, fpress);
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_retrieve_masstrans(grid_t               *G,
                            int                   np,
                            struct cfs_tpfa_data *h,
                            double               *masstrans_f)
/* ---------------------------------------------------------------------- */
{
    memcpy(masstrans_f, h->pimpl->masstrans_f,
           np * G->number_of_faces * sizeof *masstrans_f);
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_retrieve_gravtrans(grid_t               *G,
                            int                   np,
                            struct cfs_tpfa_data *h,
                            double               *gravtrans_f)
/* ---------------------------------------------------------------------- */
{
    memcpy(gravtrans_f, h->pimpl->gravtrans_f,
           np * G->number_of_faces * sizeof *gravtrans_f);
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_expl_mass_transport(grid_t       *G,
                             int           np,
                             double        dt,
                             const double *porevol,
                             const double *masstrans_f,
                             const double *gravtrans_f,
                             const double *cpress,
                             double       *surf_vol)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f, c2, p;
    double dp, dz;

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f  = G->cell_faces[i];

            if ((c2 = G->face_cells[2*f + 0]) == c) {
                c2  = G->face_cells[2*f + 1];
            }

            if (c2 >= 0) {
                dp = cpress[c] - cpress[c2];

                for (p = 0; p < np; p++) {
                    dz  = masstrans_f[f*np + p] * dp;
                    dz += gravtrans_f[f*np + p]     ;

                    surf_vol[c*np + p] -= dz * dt / porevol[c];
                }
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
void
cfs_tpfa_destroy(struct cfs_tpfa_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        csrmatrix_delete(h->A);
        impl_deallocate (h->pimpl);
    }

    free(h);
}
