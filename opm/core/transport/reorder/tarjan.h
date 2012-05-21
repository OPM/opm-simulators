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

/**
 * \file
 *
 * Simple implementation of of Tarjan's algorithm for computing the
 * strongly connected components of a directed graph, \f$G(V,E)\f$.
 * Run-time complexity is \f$O(|V| + |E|)\f$.
 *
 * The implementation is based on
 * "http://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm".
 */

#ifndef TARJAN_H_INCLUDED
#define TARJAN_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */

/**
 * Compute the strongly connected components of a directed graph,
 * \f$G(V,E)\f$.
 *
 * The components are returned in reverse topological sorted sequence.
 *
 * \param[in] nv Number of graph vertices.
 *
 * \param[in] ia
 * \param[in] ja adjacency matrix for directed graph in compressed sparse row
 *               format: vertex i has directed edges to vertices ja[ia[i]],
 *                ..., ja[ia[i+1]-1].
 *
 * \param[out] vert Permutation of vertices into topologically sorted
 *                  sequence of strong components (i.e., loops).
 *                  Array of size <CODE>nv</CODE>.
 *
 * \param[out] comp Pointers to start of each strongly connected
 *                  component in vert, the i'th component has vertices
 *                  vert[comp[i]], ..., vert[comp[i+1] - 1].  Array of
 *                  size <CODE>nv + 1</CODE>.
 *
 * \param[out] ncomp Number of strong components.  Pointer to a single
 *                   <CODE>int</CODE>.
 *
 * \param[out] work Pointer to a scratch array represented as a block
 *                  of memory capable of holding <CODE>3 * nv</CODE>
 *                  elements of type <CODE>int</CODE>.
 */
void
tarjan(int        nv   ,
       const int *ia   ,
       const int *ja   ,
       int       *vert ,
       int       *comp ,
       int       *ncomp,
       int       *work );

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* TARJAN_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
