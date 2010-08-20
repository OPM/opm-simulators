/*======================================================================

  File: dfs.c

  Created: Tue May  6 09:27:26 2008

  Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>

  Revision: $Id$
  ====================================================================*/
#include <assert.h>

/* 
 * Assign color (nonnegative number) to each connected component of graph 
 */
void dfs (int size, int *ia, int *ja, int *ncolors, int *color, int* work)
{
    int i, c;
    enum {UNVISITED = -1, VISITED = -2};
    int *stack  = work;
    int *count  = work + size;
    int *bottom = stack;

    *ncolors = 0; /* colors are nonnegative */
  
    for (i=0; i<size; ++i) {
        color [i] = UNVISITED;
        count [i] = ia[i+1]-ia[i];
    }

    /* Push seeds on stack */
    for (i=0; i<size; ++i) {
        if(color[i] >= 0) {  /* FINISHED */
            continue;
        }

        *stack++ = i; /* push i */
        color[i] = VISITED;

        while ( stack != bottom ) {
            c = *(stack-1); /* peek */
            
            if (count[c] > 0){
                int child = ja[ia[c] + count[c]-1];
                count[c]--;
                
                if (color[child] == UNVISITED) {
                    *stack++ = child;
                    color[c] = VISITED;
                }
   
            } else {
                color[c] = *ncolors;
                --stack; /* pop c */
            }
        }
        ++*ncolors;
    }
}




#if TEST
#include <stdlib.h>
#include <stdio.h>


/* Test code.  Reads a sparse matrix from a file specified on the   */
/* command line with file format "m n\n i1 j1 v1\n i2 j2 v2\n ...". */
/* Computes reorder sequence                                        */
int main (int argc, char *argv [])
{
    int      *color, *work;
    int       j, ncolors;
#if 0
    int size = 8;
    int ia[] = {0, 1, 2, 4, 5, 5, 7, 7, 8};
    int ja[] = {1, 2, 0, 3, 4, 5, 6, 6};
#else
    int size = 3;
    int ia[] = {0,2,5,7};
    int ja[] = {0,1,1,0,2,2,1};
#endif


    color = malloc (size * sizeof *color);
    work  = malloc (2*size * sizeof *work);
    dfs(size, ia, ja, &ncolors, color, work);


    fprintf(stderr, "ncolors = %d\n", ncolors);
    for (j=0; j<size; ++j) {
        fprintf(stderr, "%d\n", color[j]);   
    }
 

    free (color);
    free (work);
    return 0;
}
#endif

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
