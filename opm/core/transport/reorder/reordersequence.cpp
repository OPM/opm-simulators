/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdio.h>
#include <stdlib.h>
#ifdef MATLAB_MEX_FILE
#include "grid.h"
#include "reordersequence.h"
#include "tarjan.h"
#else
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/transport/reorder/tarjan.h>
#endif

#include <algorithm>
#include <cmath>

struct SortByAbsFlux
{
    SortByAbsFlux(const double* flux)
        : flux_(flux)
    {}
    bool operator() (int f1, int f2)
    {
	return std::fabs(flux_[f1]) < std::fabs(flux_[f2]);
    }
    const double* flux_;
};

/* Construct adjacency matrix of upwind graph wrt flux.  Column
   indices are not sorted. */
static void
make_upwind_graph(int nc, int *cellfaces, int *faceptr, int *face2cell,
                  const double *flux, int *ia, int *ja, int *work)
{
    /* Using topology (conn, cptr), and direction, construct adjacency
       matrix of graph. */

    int i, j, p, f, positive_sign, boundaryface;
    double theflux;

    /* For each face, store upwind cell in work array */
    for (i=0; i<nc; ++i)
    {
        for (j=faceptr[i]; j<faceptr[i+1]; ++j)
        {
            f  = cellfaces[j];
            positive_sign = (i == face2cell[2*f]);
            theflux = positive_sign ? flux[f] : -flux[f];

            if ( theflux > 0  )
            {
                /* i is upwind cell for face f */
                work[f] = i;
            }
        }
    }

    /* Fill ia and ja */
    p = 0;
    ia[0] = p;
    for (i=0; i<nc; ++i)
    {
        for (j=faceptr[i]; j<faceptr[i+1]; ++j)
        {

            f  = cellfaces[j];
            boundaryface = (face2cell[2*f+0] == -1) ||
                           (face2cell[2*f+1] == -1);

            if ( boundaryface )
            {
                continue;
            }

            positive_sign = (i == face2cell[2*f]);
            theflux = positive_sign ? flux[f] : -flux[f];

            if ( theflux < 0)
            {
                // ja[p++] = f;
		ja[p++] = work[f];
            }
        }
        ia[i+1] = p;
	// std::sort(ja + ia[i], ja + ia[i+1], SortByAbsFlux(flux));
	// for (j = ia[i]; j < ia[i+1]; ++j) {
	//     ja[j] = work[ja[j]];
	// }
    }
}

static void
compute_reorder_sequence(int nc, int nf, int *cellfaces, int *faceptr, int *face2cell,
                         const double *flux, int *sequence, int *components, int *ncomponents)
{
    int *ia, *ja;
    int *work;
    int sz = nf;
    if (nf < 3*nc)
    {
        sz = 3*nc;
    }

    work = (int*) malloc( sz    * sizeof *work);
    ia   = (int*) malloc((nc+1) * sizeof *ia);
    ja   = (int*) malloc( nf    * sizeof *ja); /* A bit too much... */


    make_upwind_graph(nc, cellfaces, faceptr, face2cell,
                      flux, ia, ja, work);

    tarjan (nc, ia, ja, sequence, components, ncomponents, work);

    free(ja);
    free(ia);
    free(work);
}


static void
compute_reorder_sequence_graph(int nc, int nf, int *cellfaces, int *faceptr, int *face2cell,
			       const double *flux, int *sequence, int *components, int *ncomponents,
			       int *ia, int *ja)
{
    int *work;
    int sz = nf;
    if (nf < 3*nc)
    {
        sz = 3*nc;
    }

    work = (int*) malloc( sz    * sizeof *work);


    make_upwind_graph(nc, cellfaces, faceptr, face2cell,
                      flux, ia, ja, work);

    tarjan (nc, ia, ja, sequence, components, ncomponents, work);

    free(work);
}


void compute_sequence(const struct UnstructuredGrid *grid, const double *flux,
                      int *sequence,
                      int *components, int *ncomponents)
{

    compute_reorder_sequence(grid->number_of_cells,
                             grid->number_of_faces,
                             grid->cell_faces,
                             grid->cell_facepos,
                             grid->face_cells,
                             flux,
                             sequence,
                             components,
                             ncomponents);
}

void compute_sequence_graph(const struct UnstructuredGrid *grid, const double *flux,
			    int *sequence,
			    int *components, int *ncomponents,
			    int *ia, int *ja)
{

    compute_reorder_sequence_graph(grid->number_of_cells,
				   grid->number_of_faces,
				   grid->cell_faces,
				   grid->cell_facepos,
				   grid->face_cells,
				   flux,
				   sequence,
				   components,
				   ncomponents,
				   ia, ja);
}


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
