/*
Copyright (C) 2012 (c) Jostein R. Natvig <jostein natvig at gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <assert.h>
#include <stddef.h>

#ifdef MATLAB_MEX_FILE
#include "tarjan.h"
#else
#include <opm/core/transport/reorder/tarjan.h>
#endif


static void
clear_vector(size_t n, int *v)
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] = 0; }
}

static int min(int a, int b){ return a < b? a : b;}

/*
  Compute the strong components of directed graph G(edges, vertices),
  return components in reverse topological sorted sequence.
  Complexity O(|vertices|+|edges|). See "http://en.wikipedia.org/wiki/
  Tarjan's_strongly_connected_components_algorithm".

  nv    - number of vertices

  ia,ja - adjacency matrix for directed graph in compressed sparse row
          format: vertex i has directed edges to vertices ja[ia[i]],
          ..., ja[ia[i+1]-1].

  vert  - permutation of vertices into topologically sorted sequence of
          strong components (i.e., loops).

  comp  - pointers to start of each strongly connected component in
          vert, the i'th component has vertices vert[comp[i]], ...,
          vert[comp[i+1]-1].

  ncomp - number of strong components.

  work  - block of memory of size 3*nv*sizeof(int).
 */

/*--------------------------------------------------------------------*/
void
tarjan (int nv, const int *ia, const int *ja, int *vert, int *comp,
        int *ncomp, int *work)
/*--------------------------------------------------------------------*/
{
    /* Hint: end of VERT and COMP are used as stacks. */

    enum {DONE=-2, REMAINING=-1};
    int  c, v, seed, child;
    int  i;

    int *stack  = comp + nv;
    int *bottom = stack;
    int *cstack = vert + nv-1;

#if !defined(NDEBUG)
    int *cbottom = cstack;
#endif

    int  t      = 0;
    int  pos    = 0;

    int *time   = work;
    int *link   = time + nv;
    int *status = link + nv; /* dual usage... */

    clear_vector(3 * ((size_t) nv), work);
    clear_vector(1 * ((size_t) nv), vert);
    clear_vector(1 + ((size_t) nv), comp);

    /* Init status all vertices */
    for (i=0; i<nv; ++i)
    {
        status[i] = REMAINING;
    }

    *ncomp  = 0;
    *comp++ = pos;

    seed = 0;
    while (seed < nv)
    {
        if (status[seed] == DONE)
        {
            ++seed;
            continue;
        }

        /* push seed */
        *stack-- = seed;

        t = 0;

        while ( stack != bottom )
        {
            /* peek c */
            c = *(stack+1);

            assert(status[c] != DONE);
            assert(status[c] >= -2);

            if (status[c] == REMAINING)
            {
                /* number of descendants of c */
                status[c] = ia[c+1]-ia[c];
                time[c]   = link[c] = t++;

                /* push c on strongcomp stack */
                *cstack-- = c;
            }



            /* if all descendants are processed */
            if (status[c] == 0)
            {

                /* if c is root of strong component */
                if (link[c] == time[c])
                {
                    do
                    {
                        assert (cstack != cbottom);

                        /* pop strong component stack */
                        v         = *++cstack;
                        status[v] = DONE;

                        /* store vertex in VERT */
                        vert[pos++]  = v;
                    }
                    while ( v != c );

                    /* store end point of component */
                    *comp++ = pos;
                    ++*ncomp;
                }

                /* pop c */
                ++stack;

                if (stack != bottom)
                {
                    link[*(stack+1)] = min(link[*(stack+1)], link[c]);
                }
            }



            /* if there are more descendants to consider */
            else
            {
                assert(status[c] > 0);

                child = ja[ia[c] + status[c]-1];
                /* decrement descendant count of c*/
                --status[c];

                if (status[child] == REMAINING)
                {
                    /* push child */
                    *stack-- = child;

                }
                else if (status[child] >= 0)
                {
                    link[c] = min(link[c], time[child]);

                }
                else
                {
                    assert(status[child] == DONE);
                }
            }
        }
        assert (cstack == cbottom);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
