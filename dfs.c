/*======================================================================

  File: dfs.c

  Created: Tue May  6 09:27:26 2008

  Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>

  Revision: $Id$
  ====================================================================*/
#include <assert.h>

void dfs (int size, int *ia, int *ja, int *ncolors, int *color, int* work)
{
   int i, c;

   int *stack  = work;
   int *bottom = stack;

   *ncolors = 0; /* colors are nonnegative */
  
   for (i=0; i<size; ++i) {
      color [i] = -(ia[i+1]-ia[i]+1);
   }

   /* Push seeds on stack */
   for (i=0; i<size; ++i) {
      if(color[i] >= 0) {
         continue;
      }

      *stack++ = i; /* push i */
     
      while ( stack != bottom ) {
         c = *(stack-1); /* peek */

         if (color[c] < 0){
            *stack++ = ja[ia[c]-color[c]-2];
            ++color[c];
   
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
   int size = 8;
   int ia[] = {0, 1, 2, 4, 5, 5, 7, 7, 8};
   int ja[] = {1, 2, 0, 3, 4, 5, 6, 6};

   color = malloc (size * sizeof *color);
   work  = malloc (size * sizeof *work);
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
