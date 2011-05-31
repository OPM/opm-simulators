/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dfs.h"
#include "partition.h"


#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/* [cidx{1:ndims}] = ind2sub(size, idx) */
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


/* sub2ind(size, cidx{1:ndims}) */
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


/* Partition 'nc' fine-scale Cartesian indices 'idx' from a box of
 * dimension 'fine_d' into a coarse-scale box of dimension 'coarse_d'.
 *
 * Store partition in vector 'p' (assumed to hold at least 'nc'
 * slots).
 *
 * Allocates a tiny work array to hold 'ndims' ints.  Returns 'nc' if
 * successful and -1 if unable to allocate the work array. */
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


/* Renumber blocks to create contiguous block numbers from 0..n-1
 * (in other words: remove empty coarse blocks).
 *
 * Returns maximum new block number if successful and -1 if not. */
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


/* Free memory resources for block->cell map. */
/* ---------------------------------------------------------------------- */
void
partition_deallocate_inverse(int *pi, int *inverse)
/* ---------------------------------------------------------------------- */
{
    free(inverse);
    free(pi);
}


/* Allocate memory for block->cell map (CSR representation).  Highest
 * block number is 'max_bin'.  Grid contains 'nc' cells.
 *
 * Returns 'nc' (and sets CSR pointer pair (*pi, *inverse)) if
 * successful, -1 and pointers to NULL if not. */
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


/* Invert cell->block mapping 'p' (partition vector) to create
 * block->cell mapping (CSR representation, pointer pair (pi,inverse)). */
/* ---------------------------------------------------------------------- */
void
partition_invert(int nc, const int *p, int *pi, int *inverse)
/* ---------------------------------------------------------------------- */
{
    int nbin, b, i;

    nbin = max_block(nc, p) + 1; /* Adjust for bin 0 */

    /* Zero start pointers */
    memset(pi, 0, (nbin + 1) * sizeof *pi);

    /* Count elements per bin */
    for (i = 0; i < nc; i++) { pi[ p[i] + 1 ]++; }

    for (b = 1; b <= nbin; b++) {
        pi[0] += pi[b];
        pi[b]  = pi[0] - pi[b];
    }

    /* Insert bin elements whilst deriving start pointers */
    for (i = 0; i < nc; i++) {
        inverse[ pi[ p[i] + 1 ] ++ ] = i;
    }

    /* Assert basic sanity */
    assert (pi[nbin] == nc);
    pi[0] = 0;
}


/* Create local cell numbering, within the cell's block, for each
 * global cell. */
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


/* Release memory resources for internal cell-to-cell connectivity
 * (CSR representation). */
/* ---------------------------------------------------------------------- */
static void
partition_destroy_c2c(int *pc2c, int *c2c)
/* ---------------------------------------------------------------------- */
{
    free(c2c);  free(pc2c);
}


/* Create symmetric cell-to-cell (internal) connectivity for domain
 * containing 'nc' cells.  CSR representation (*pc2c,*c2c).
 *
 * Neighbourship 'neigh' is 2*nneigh array such that cell neigh[2*i+0]
 * is connected to cell neigh[2*i+1] for all i=0:nneigh-1.
 *
 * Negative 'neigh' entries represent invalid cells (outside domain).
 *
 * Returns 'nc' (and sets pointer pair) if successful, 0 (and pointer
 * pair to NULL) if not. */
/* ---------------------------------------------------------------------- */
static int
partition_create_c2c(int nc, int nneigh, const int *neigh,
                     int **pc2c, int **c2c)
/* ---------------------------------------------------------------------- */
{
    int i, ret, c1, c2;

    *pc2c = calloc(nc + 1, sizeof **pc2c);

    if (*pc2c != NULL) {
        for (i = 0; i < nneigh; i++) {
            if ((neigh[2*i + 0] >= 0) && (neigh[2*i + 1] >= 0)) {
                /* Symmetric Laplace matrix (undirected graph) */
                (*pc2c)[ neigh[2*i + 0] + 1 ] ++;
                (*pc2c)[ neigh[2*i + 1] + 1 ] ++;
            }
        }

        for (i = 1; i <= nc; i++) {
            (*pc2c)[i] += 1;    /* Self connection */

            (*pc2c)[0] += (*pc2c)[i];
            (*pc2c)[i]  = (*pc2c)[0] - (*pc2c)[i];
        }

        *c2c = malloc((*pc2c)[0] * sizeof **c2c);

        if (*c2c != NULL) {
            /* Self connections */
            for (i = 0; i < nc; i++) {
                (*c2c)[ (*pc2c)[i + 1] ++ ] = i;
            }

            for (i = 0; i < nneigh; i++) {
                c1 = neigh[2*i + 0];
                c2 = neigh[2*i + 1];

                if ((c1 >= 0) && (c2 >= 0)) {
                    /* Symmetric Laplace matrix (undirected graph) */
                    (*c2c)[ (*pc2c)[ c1 + 1 ] ++ ] = c2;
                    (*c2c)[ (*pc2c)[ c2 + 1 ] ++ ] = c1;
                }
            }

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


/* Release dfs() memory resources. */
/* ---------------------------------------------------------------------- */
static void
deallocate_dfs_arrays(int *ia, int *ja, int *colour, int *work)
/* ---------------------------------------------------------------------- */
{
    free(work);  free(colour);  free(ja);  free(ia);
}


/* Allocate dfs() memory resources to support graph containing 'n'
 * nodes and (at most) 'nnz' total connections.  Return 'n' if
 * successful (and set pointers) and 0 (and set pointers to NULL) if
 * not. */
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


/* Compute maximum number of cells (*max_blk_cells) and cell-to-cell
 * connections (*max_blk_conn) over all blocks. */
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


/* Create block-internal (symmetric) connectivity graph (CSR
 * representation ia,ja) for connected component labelling (used in
 * splitting disconnected blocks). */
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
                ia[loc[c] + 1] ++;
            }
        }
    }

    assert (ia[0] == 0);

    for (i = 1; i <= nc; i++) {
        ia[0] += ia[i];
        ia[i]  = ia[0] - ia[i];
    }

    for (i = pb2c[b]; i < pb2c[b + 1]; i++) {
        c = b2c[i];

        /* Create connections (self conn automatic) */
        for (j = pc2c[c]; j < pc2c[c + 1]; j++) {
            if (p[c2c[j]] == b) {
                ja[ ia[loc[c] + 1] ++ ] = loc[c2c[j]];
            }
        }
    }

    ia[0] = 0;
}


/* Split disconnected coarse blocks.  Preserve block numbering where
 * possible.
 *
 * Neighbourship definition 'neigh' is pointer to 2*nneigh array such
 * that cell neigh[2*i+0] is connected to cell neigh[2*i+1] for all
 * i=0:nneigh-1.  Negative entries in 'neigh' represent invalid cells
 * (outside domain).
 *
 * Returns number of new blocks (0 if all blocks internally connected)
 * if successful and -1 otherwise. */
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
