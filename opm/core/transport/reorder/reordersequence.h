/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */
#ifndef REORDERSEQUENCE_H_INCLUDED
#define REORDERSEQUENCE_H_INCLUDED

/**
 * \file
 *
 * Facilities for computing a causal permutation of the cells in an
 * UnstructuredGrid such that fluid transport may be subsequently
 * solved by going from sources and up-stream cells to sinks and
 * down-stream cells.
 */

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

struct UnstructuredGrid;


/**
 * Compute causal permutation sequence of grid cells with respect to
 * specific Darcy flux field.
 *
 * \param[in] grid Grid structure for which to compute causal cell
 *                 permutation.
 *
 * \param[in] flux Darcy flux field.  One scalar value for each
 *                 interface/connection in the grid, including the boundary.  We assume that
 *                 <CODE>flux[f]</CODE> is positive if the flow is
 *                 from cell <CODE>grid->face_cells[2*f + 0]</CODE> to
 *                 cell <CODE>grid->face_cells[2*f + 1]</CODE>.
 *
 * \param[out] sequence
 *                 Causal grid cell permutation.  Ordered according to
 *                 topological sorting of the strongly connected
 *                 components of the Darcy flux upwind graph.  Array
 *                 of size <CODE>grid->number_of_cells</CODE>.
 *
 * \param[out] components
 *                 Indirection pointers that describe the strongly
 *                 connected components (i.e., loops) of the Darcy
 *                 flux upwind graph.  Specifically, the \f$i\f$'th
 *                 strongly connected components constitutes cells
 *                 <CODE>sequence[components[i] ... components[i + 1]
 *                 - 1]</CODE>.
 *
 *                 In the ideal case of a perfectly reordered upwind
 *                 graph, this array will hold
 *                 <CODE>grid->number_of_cells + 1</CODE> elements,
 *                 and the relation <CODE>components[i + 1] -
 *                 components[i] == 1</CODE> holds for all \f$i\f$.
 *                 Thus, the <CODE>components</CODE> parameter must
 *                 point to an array of <CODE>grid->number_of_cells +
 *                 1</CODE> elements.
 *
 * \param[out] ncomponents
 *                 Number of strongly connected components.  Pointer
 *                 to a single integer.  The <CODE>components</CODE>
 *                 output are valid for all <CODE>i = 0
 *                 ... *ncomponents - 1</CODE>. Furhtermore, the number of
 *                 components must be in the interval <CODE>[1
 *                 .. grid->number_of_cells]</CODE>.
 */
void
compute_sequence(const struct UnstructuredGrid *grid       ,
                 const double                  *flux       ,
                 int                           *sequence   ,
                 int                           *components ,
                 int                           *ncomponents);


/**
 * Compute causal permutation sequence of grid cells with respect to
 * specific Darcy flux field.  Also return the permuted upwind graph.
 *
 * \param[in] grid Grid structure for which to compute causal cell
 *                 permutation.
 *
 * \param[in] flux Darcy flux field.  One scalar value for each
 *                 interface/connection in the grid, including the boundary.  We assume that
 *                 <CODE>flux[f]</CODE> is positive if the flow is
 *                 from cell <CODE>grid->face_cells[2*f + 0]</CODE> to
 *                 cell <CODE>grid->face_cells[2*f + 1]</CODE>.
 *
 * \param[out] sequence
 *                 Causal grid cell permutation.  Ordered according to
 *                 topological sorting of the strongly connected
 *                 components of the Darcy flux upwind graph.  Array
 *                 of size <CODE>grid->number_of_cells</CODE>.
 *
 * \param[out] components
 *                 Indirection pointers that describe the strongly
 *                 connected components (i.e., loops) of the Darcy
 *                 flux upwind graph.  Specifically, the \f$i\f$'th
 *                 strongly connected components constitutes cells
 *                 <CODE>sequence[components[i] ... components[i + 1]
 *                 - 1]</CODE>.
 *
 *                 In the ideal case of a perfectly reordered upwind
 *                 graph, this array will hold
 *                 <CODE>grid->number_of_cells + 1</CODE> elements,
 *                 and the relation <CODE>components[i + 1] -
 *                 components[i] == 1</CODE> holds for all \f$i\f$.
 *                 Thus, the <CODE>components</CODE> parameter must
 *                 point to an array of <CODE>grid->number_of_cells +
 *                 1</CODE> elements.
 *
 * \param[out] ncomponents
 *                 Number of strongly connected components.  Pointer
 *                 to a single integer.  The <CODE>components</CODE>
 *                 output are valid for all <CODE>i = 0
 *                 ... *ncomponents - 1</CODE>. Furhtermore, the number of
 *                 components must be in the interval <CODE>[1
 *                 .. grid->number_of_cells]</CODE>.
 *
 * \param[out] ia  Indirection pointers into <CODE>ja</CODE> that
 *                 define the "upwind" cells that influence a
 *                 particular grid cell.
 *
 * \param[out] ja  Compressed-sparse representation of the upwind
 *                 graph.  Specifically, the upwind cells that
 *                 influence cell \f$i\f$ are <CODE>ja[ia[i]
 *                 .. ia[i+1]-1]</CODE>.  Array of size at least equal
 *                 to the number of internal faces of
 *                 <CODE>grid</CODE>.  The number
 *                 <CODE>grid->number_of_faces</CODE> is an upper
 *                 bound of the array size.
 */
void
compute_sequence_graph(const struct UnstructuredGrid *grid       ,
                       const double                  *flux       ,
                       int                           *sequence   ,
                       int                           *components ,
                       int                           *ncomponents,
                       int                           *ia         ,
                       int                           *ja         );

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* REORDERSEQUENCE_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
