#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

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
void
partition_invert(int nc, const int *p, int *pi, int *inverse)
/* ---------------------------------------------------------------------- */
{
    int nbin, b, i, j, tmp;

    nbin = 0;
    for (i = 0; i < nc; i++) {
        nbin = MAX(nbin, p[i]);
    }
    nbin += 1;            /* Adjust for bin 0 */

    /* Zero start pointers */
    for (b = 0; b < nbin; b++) { pi[b] = 0; }

    /* Count elements per bin */
    for (i = 0; i < nc  ; i++) { pi[ p[i] ]++; }

    /* Derive start pointers for b=1:nbin */
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
    for (b = 0; b < nbin; b++) {
        i = pi[b + 0] + 0;
        j = pi[b + 1] - 1;

        while (i < j) {
            /* Swap reverse (lower <-> upper) */
            tmp        = inverse[i];
            inverse[i] = inverse[j];
            inverse[j] = tmp;

            i += 1;             /* Increase lower bound */
            j -= 1;             /* Decrease upper bound */
        }
    }
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
