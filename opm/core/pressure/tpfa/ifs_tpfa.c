#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>


struct ifs_tpfa_impl {
    double *fgrav;              /* Accumulated grav contrib/face */

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
impl_allocate(struct UnstructuredGrid *G)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_impl *new;

    size_t ddata_sz;

    ddata_sz  = 2 * G->number_of_cells;                /* b, x */
    ddata_sz += 1 * G->number_of_faces;                /* fgrav */

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
ifs_tpfa_construct_matrix(struct UnstructuredGrid *G)
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2;
    size_t nnz;

    struct CSRMatrix *A;

    A = csrmatrix_new_count_nnz(G->number_of_cells);

    if (A != NULL) {
        /* Self connections */
        for (c1 = 0; c1 < G->number_of_cells; c1++) {
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

        nnz = csrmatrix_new_elms_pushback(A);
        if (nnz == 0) {
            csrmatrix_delete(A);
            A = NULL;
        }
    }

    if (A != NULL) {
        /* Fill self connections */
        for (c1 = 0; c1 < G->number_of_cells; c1++) {
            A->ja[ A->ia[ c1 + 1 ] ++ ] = c1;
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

        assert ((size_t) A->ia[ G->number_of_cells ] == nnz);

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



/* ======================================================================
 * Public interface below separator.
 * ====================================================================== */

/* ---------------------------------------------------------------------- */
struct ifs_tpfa_data *
ifs_tpfa_construct(struct UnstructuredGrid *G)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->pimpl = impl_allocate(G);
        new->A     = ifs_tpfa_construct_matrix(G);

        if ((new->pimpl == NULL) || (new->A == NULL)) {
            ifs_tpfa_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        new->b = new->pimpl->ddata;
        new->x = new->b             + new->A->m;

        new->pimpl->fgrav = new->x  + new->A->m;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_assemble(struct UnstructuredGrid               *G,
                  const struct ifs_tpfa_forces *F,
                  const double         *trans,
                  const double         *gpress,
                  struct ifs_tpfa_data *h)
/* ---------------------------------------------------------------------- */
{
    int c1, c2, c, i, f, j1, j2;

    double s;

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

        if ((F != NULL) && (F->src != NULL)) {
            h->b[c] += F->src[c];
        }
    }

    h->A->sa[0] *= 2;
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_press_flux(struct UnstructuredGrid               *G,
                    const double         *trans,
                    struct ifs_tpfa_data *h,
                    double               *cpress,
                    double               *fflux)
/* ---------------------------------------------------------------------- */
{
    int    c1, c2, f;
    double dh;

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
