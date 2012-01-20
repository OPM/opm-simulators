/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <string.h>
#include <assert.h>

#ifdef MATLAB_MEX_FILE
#include "tarjan.h"
#else
#include <opm/core/transport/reorder/tarjan.h>
#endif




static int min(int a, int b){ return a < b? a : b;}

/* Improved (or obfuscated) version uses less memory
 *
 * Use end of P and Q as stack:
 *   push operation is *s--=elm,
 *   pop  operation is elm=*++s,
 *   peek operation is *(s+1)
 */
void
tarjan (int size, int *ia, int *ja, int *P, int *Q, int *ncomp,
        int *work)
/*--------------------------------------------------------------------*/
{
    enum {DONE=-2, REMAINING=-1};
    int c,v,seed,child;
    int i;

    int *stack    = Q + size,   *bottom  = stack;
    int *cstack   = P + size-1, *cbottom = cstack;

    int t         = 0;
    int pos       = 0;

    int   *time   = work;
    int   *link   = (int *)  time + size;
    int   *status = (int*)   link + size; /* dual usage... */

    (void) cbottom;

    memset(work, 0, 3*size    *  sizeof *work);
    memset(P,    0,   size    *  sizeof *P   );
    memset(Q,    0,  (size+1) *  sizeof *Q   );

    /* Init status all nodes */
    for (i=0; i<size; ++i)
    {
        status[i] = REMAINING;
    }

    *ncomp = 0;
    *Q++ = pos;

    seed = 0;
    while (seed < size)
    {
        if (status[seed] == DONE)
        {
            ++seed;
            continue;
        }

        *stack-- = seed;     /* push seed */

        t = 0;

        while ( stack != bottom )
        {
            c = *(stack+1);        /* peek c */

            assert(status[c] != DONE);
            assert(status[c] >= -2);

            if (status[c] == REMAINING)
            {
                status[c] = ia[c+1]-ia[c]; /* number of descendants of c */
                time[c]   = link[c] = t++;
                *cstack-- = c;             /* push c on strongcomp stack */
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

                        v         = *++cstack; /* pop strong component stack */
                        status[v] = DONE;
                        P[pos++]  = v;          /* store vertex in P */
                    }
                    while ( v != c );

                    *Q++ = pos;       /* store end point of component */
                    ++*ncomp;
                }

                ++stack;  /* pop c */

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
                --status[c];           /* Decrement descendant count of c*/

                if (status[child] == REMAINING)
                {
                    *stack-- = child; /* push child */

                }
                else if (status[child] >= 0)
                {
                    link[c] = min(link[c], time[child]);

                }
                else
                {
                    assert(status[child] = DONE);
                }
            }
        }
        assert (cstack == cbottom);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
