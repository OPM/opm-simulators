#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dfs.h"
#include "partition.h"


#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/* ---------------------------------------------------------------------- */
static void
partition_coord_idx(int ndims, int idx, const int *size, int *cidx)
/* ---------------------------------------------------------------------- */
{
    int i;

    for (i = 0; i < ndims; i++) {
        cidx[i]  = idx % size[i];
        idx     /=       size[i];
    }

    assert (idx == 0);
}


/* ---------------------------------------------------------------------- */
static int
partition_lin_idx(int ndims, const int *size, const int *cidx)
/* ---------------------------------------------------------------------- */
{
    int i, idx;

    idx = cidx[ndims - 1];
    for (i = ndims - 2; i >= 0; i--) {
        idx = cidx[i] + size[i]*idx;
    }

    return idx;
}


/* ---------------------------------------------------------------------- */
/* Load-balanced linear distribution.
 *
 * See Eric F. Van de Velde, Concurrent Scientific Computing,
 * 1994, Springer Verlag, p. 54 (Sect. 2.3) for details.  */
static void
partition_loadbal_lin_dist(int ndims, const int *size, const int *nbins,
                           int *idx)
/* ---------------------------------------------------------------------- */
{
    int i, L, R, b1, b2;

    for (i = 0; i < ndims; i++) {
        L = size[i] / nbins[i]; /* # entities per bin */
        R = size[i] % nbins[i]; /* # bins containing one extra entity */

        b1 =  idx[i]      / (L + 1);
        b2 = (idx[i] - R) /  L     ;

        idx[i] = MAX(b1, b2);
    }
}


/* ---------------------------------------------------------------------- */
int
partition_unif_idx(int ndims, int nc,
                   const int *fine_d, const int *coarse_d, const int *idx,
                   int *p)
/* ---------------------------------------------------------------------- */
{
    int c, ret, *ix;

    ix = malloc(ndims * sizeof *ix);

    if (ix != NULL) {
        for (c = 0; c < nc; c++) {
            partition_coord_idx(ndims, idx[c], fine_d, ix);

            partition_loadbal_lin_dist(ndims, fine_d, coarse_d, ix);

            p[c] = partition_lin_idx(ndims, coarse_d, ix);
        }

        ret = nc;
    } else {
        ret = -1;
    }

    free(ix);

    return ret;
}


/* ---------------------------------------------------------------------- */
int
partition_compress(int n, int *p)
/* ---------------------------------------------------------------------- */
{
    int ret, i, max, *compr;

    max = -1;
    for (i = 0; i < n; i++) {
        assert (0 <= p[i]);     /* Only non-neg partitions (for now?). */
        max = MAX(max, p[i]);
    }

    compr = calloc(max + 1, sizeof *compr);

    if (compr != NULL) {
        for (i = 0; i < n; i++) { compr[p[i]]++; }

        compr[0] = -1 + (compr[0] > 0);
        for (i = 1; i <= max; i++) {
            compr[i] = compr[i - 1] + (compr[i] > 0);
        }

        for (i = 0; i < n; i++) { p[i] = compr[p[i]]; }

        ret = compr[max];
    } else {
        ret = -1;
    }

    free(compr);

    return ret;
}


/* ---------------------------------------------------------------------- */
void
partition_deallocate_inverse(int *pi, int *inverse)
/* ---------------------------------------------------------------------- */
{
    free(inverse);
    free(pi);
}


/* ---------------------------------------------------------------------- */
int
partition_allocate_inverse(int nc, int max_bin,
                           int **pi, int **inverse)
/* ---------------------------------------------------------------------- */
{
    int nbin, ret, *ptr, *i;

    nbin = max_bin + 1;

    ptr  = malloc((nbin + 1) * sizeof *ptr);
    i    = malloc(nc         * sizeof *i  );

    if ((ptr == NULL) || (i == NULL)) {
        partition_deallocate_inverse(ptr, i);

        *pi      = NULL;
        *inverse = NULL;

        ret = 0;
    } else {
        *pi      = ptr;
        *inverse = i;

        ret = nc;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
reverse_bins(int nbin, const int *pbin, int *elements)
/* ---------------------------------------------------------------------- */
{
    int b, i, j, tmp;

    for (b = 0; b < nbin; b++) {
        i = pbin[b + 0] + 0;
        j = pbin[b + 1] - 1;

        while (i < j) {
            /* Swap reverse (lower <-> upper) */
            tmp         = elements[i];
            elements[i] = elements[j];
            elements[j] = tmp;

            i += 1;             /* Increase lower bound */
            j -= 1;             /* Decrease upper bound */
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
max_block(int nc, const int *p)
/* ---------------------------------------------------------------------- */
{
    int m, i;

    m = -1;

    for (i = 0; i < nc; i++) {
        m = MAX(m, p[i]);
    }

    return m;
}


/* ---------------------------------------------------------------------- */
void
partition_invert(int nc, const int *p, int *pi, int *inverse)
/* ---------------------------------------------------------------------- */
{
    int nbin, b, i;

    nbin = max_block(nc, p) + 1; /* Adjust for bin 0 */

    /* Zero start pointers */
    for (b = 0; b < nbin; b++) { pi[b] = 0; }

    /* Count elements per bin */
    for (i = 0; i < nc  ; i++) { pi[ p[i] ]++; }

    /* Derive start pointers for b=1:nbin (== ubound for b=0:nbin-1) */
    for (b = 1; b < nbin; b++) { pi[b] += pi[b - 1]; }

    /* Set end pointer in last bin */
    assert (pi[nbin - 1] == nc);
    pi[nbin] = nc;

    /* Reverse insert bin elements whilst deriving start pointers */
    for (i = 0; i < nc; i++) {
        inverse[-- pi[ p[i] ]] = i;
    }
    assert (pi[0] == 0);

    /* Reverse the reverse order, creating final inverse mapping */
    reverse_bins(nbin, pi, inverse);
}


/* ---------------------------------------------------------------------- */
void
partition_localidx(int nbin, const int *pi, const int *inverse,
                   int *localidx)
/* ---------------------------------------------------------------------- */
{
    int b, i;

    for (b = 0; b < nbin; b++) {
        for (i = pi[b]; i < pi[b + 1]; i++) {
            localidx[ inverse[i] ] = i - pi[b];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
partition_destroy_c2c(int *pc2c, int *c2c)
/* ---------------------------------------------------------------------- */
{
    free(c2c);  free(pc2c);
}


/* ---------------------------------------------------------------------- */
static int
partition_create_c2c(int nc, int nneigh, const int *neigh,
                     int **pc2c, int **c2c)
/* ---------------------------------------------------------------------- */
{
    int i, ret;

    *pc2c = calloc(nc + 1, sizeof **pc2c);

    if (*pc2c != NULL) {
        for (i = 0; i < nneigh; i++) {
            if ((neigh[2*i + 0] >= 0) && (neigh[2*i + 1] >= 0)) {
                /* Symmetric Laplace matrix (undirected graph) */
                (*pc2c)[neigh[2*i + 0]]++;
                (*pc2c)[neigh[2*i + 1]]++;
            }
        }

        (*pc2c)[0] += 1;        /* Self connection */
        for (i = 1; i < nc; i++) {
            (*pc2c)[i] += (*pc2c)[i - 1];
            (*pc2c)[i] += 1;    /* Self connection */
        }
        (*pc2c)[nc] = (*pc2c)[nc - 1];

        *c2c = malloc((*pc2c)[nc] * sizeof **c2c);

        if (*c2c != NULL) {
            /* Self connections */
            for (i = 0; i < nc; i++) {
                (*c2c)[-- (*pc2c)[i]] = i;
            }

            for (i = 0; i < nneigh; i++) {
                if ((neigh[2*i + 0] >= 0) && (neigh[2*i + 1] >= 0)) {
                    /* Symmetric Laplace matrix (undirected graph) */
                    (*c2c)[-- (*pc2c)[neigh[2*i + 0]]] = neigh[2*i + 1];
                    (*c2c)[-- (*pc2c)[neigh[2*i + 1]]] = neigh[2*i + 0];
                }
            }

            reverse_bins(nc, *pc2c, *c2c);

            ret = nc;
        } else {
            free(*pc2c);
            *pc2c = NULL;

            ret = 0;
        }
    } else {
        *c2c = NULL;

        ret = 0;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
deallocate_dfs_arrays(int *ia, int *ja, int *colour, int *work)
/* ---------------------------------------------------------------------- */
{
    free(work);  free(colour);  free(ja);  free(ia);
}


/* ---------------------------------------------------------------------- */
static int
allocate_dfs_arrays(int n, int nnz,
                    int **ia, int **ja, int **colour, int **work)
/* ---------------------------------------------------------------------- */
{
    int ret;

    *ia     = malloc((n + 1) * sizeof **ia    );
    *ja     = malloc(nnz     * sizeof **ja    );
    *colour = malloc(n       * sizeof **colour);
    *work   = malloc(2 * n   * sizeof **work  );

    if ((*ia     == NULL) || (*ja   == NULL) ||
        (*colour == NULL) || (*work == NULL)) {
        deallocate_dfs_arrays(*ia, *ja, *colour, *work);

        *ia     = NULL;
        *ja     = NULL;
        *colour = NULL;
        *work   = NULL;

        ret = 0;
    } else {
        ret = n;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
count_block_conns(int nblk,
                  const int *pb2c, const int *b2c, const int *pc2c,
                  int *max_blk_cells, int *max_blk_conn)
/* ---------------------------------------------------------------------- */
{
    int b, i, n_blk_conn;

    *max_blk_cells = 0;
    *max_blk_conn  = 0;

    i = 0;                      /* == pb2c[0] */
    for (b = 0; b < nblk; b++) {
        n_blk_conn = 0;

        for (; i < pb2c[b + 1]; i++) {
            n_blk_conn += pc2c[b2c[i] + 1] - pc2c[b2c[i]];
        }

        *max_blk_cells = MAX(*max_blk_cells, pb2c[b + 1] - pb2c[b]);
        *max_blk_conn  = MAX(*max_blk_conn , n_blk_conn);
    }
}


/* ---------------------------------------------------------------------- */
static void
create_block_conns(int        b   ,
                   const int *p   , const int *loc,
                   const int *pb2c, const int *b2c,
                   const int *pc2c, const int *c2c,
                   int       *ia  , int       *ja )
/* ---------------------------------------------------------------------- */
{
    int nc, c, i, j;

    nc = pb2c[b + 1] - pb2c[b];

    /* Clear start pointers */
    memset(ia, 0, (nc + 1) * sizeof *ia);

    for (i = pb2c[b]; i < pb2c[b + 1]; i++) {
        c = b2c[i];   assert (loc[c] == i - pb2c[b]);

        /* Self connections inserted in partition_create_c2c()) */
        for (j = pc2c[c]; j < pc2c[c + 1]; j++) {
            if (p[c2c[j]] == b) {
                /* Connection internal to block 'b'.  Add */
                ia[loc[c]] ++;
            }
        }
    }

    assert (ia[nc] == 0);
    for (i = 1; i <= nc; i++) { ia[i] += ia[i - 1]; }

    for (i = pb2c[b]; i < pb2c[b + 1]; i++) {
        c = b2c[i];

        /* Create connections (self conn automatic) */
        for (j = pc2c[c]; j < pc2c[c + 1]; j++) {
            if (p[c2c[j]] == b) {
                ja[-- ia[loc[c]]] = loc[c2c[j]];
            }
        }
    }
    assert (ia[0] == 0);

    reverse_bins(nc, ia, ja);
}


/* ---------------------------------------------------------------------- */
int
partition_split_disconnected(int nc, int nneigh, const int *neigh, int *p)
/* ---------------------------------------------------------------------- */
{
    int inv_ok, c2c_ok, dfs_ok;
    int i, b, ret, maxblk, ncolour, max_blk_cells, max_blk_conn;

    int *pb2c, *b2c, *loc, *pc2c, *c2c;
    int *ia, *ja, *colour, *work;

    maxblk = max_block(nc, p);

    inv_ok = partition_allocate_inverse(nc, maxblk, &pb2c, &b2c);
    c2c_ok = partition_create_c2c(nc, nneigh, neigh, &pc2c, &c2c);
    loc    = malloc(nc * sizeof *loc);

    if (inv_ok && c2c_ok && (loc != NULL)) {
        partition_invert(nc, p, pb2c, b2c);
        partition_localidx(maxblk + 1, pb2c, b2c, loc);

        count_block_conns(maxblk + 1, pb2c, b2c, pc2c,
                          &max_blk_cells, &max_blk_conn);

        dfs_ok = allocate_dfs_arrays(max_blk_cells, max_blk_conn,
                                     &ia, &ja, &colour, &work);
        if (dfs_ok) {
            /* Target acquired.  Fire. */
            ret = 0;

            for (b = 0; b < maxblk + 1; b++) {
                create_block_conns(b, p, loc, pb2c, b2c, pc2c, c2c, ia, ja);

                dfs(pb2c[b + 1] - pb2c[b], ia, ja, &ncolour, colour, work);

                if (ncolour > 1) {
                    /* Block contains more than one component.  Assign
                     * new block numbers for cells in components
                     * 1:ncomp-1. */
                    for (i = pb2c[b]; i < pb2c[b + 1]; i++) {
                        if (colour[i - pb2c[b]] > 0) {
                            p[b2c[i]] = maxblk + ret + colour[i - pb2c[b]];
                        }
                    }

                    ret += ncolour - 1;
                }
            }
        } else {
            ret = -1;
        }

        deallocate_dfs_arrays(ia, ja, colour, work);
    } else {
        ret = -1;
    }

    free(loc);
    partition_destroy_c2c(pc2c, c2c);
    partition_deallocate_inverse(pb2c, b2c);

    return ret;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
