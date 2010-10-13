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

#ifndef OPM_MIMETIC_HEADER_INCLUDED
#define OPM_MIMETIC_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

void mim_ip_span_nullspace(int nf, int nconn, int d,
                           double *C,
                           double *A,
                           double *X,
                           double *work, int lwork);

void mim_ip_linpress_exact(int nf, int nconn, int d,
                           double vol, double *K,
                           double *N,
                           double *Binv,
                           double *work, int lwork);

void mim_ip_simple(int nf, int nconn, int d,
                   double v, double *K, double *C,
                   double *A, double *N,
                   double *Binv,
                   double *work, int lwork);


/** Compute the mimetic inner products given a grid and cellwise
 *  permeability tensors.
 *
 * @param ncells Number of cells in grid.
 * @param d Number of space dimensions.
 * @param max_ncf Maximum number of faces per cell.
 * @param ncf Number of faces per cell.
 * @param pconn Start indices in conn for each cell, plus end
 *              marker. The size of pconn is (ncells + 1), and for a
 *              cell i, [conn[pconn[i]], conn[pconn[i+1]]) is a
 *              half-open interval containing the indices of faces
 *              adjacent to i.
 * @param conn Cell to face mapping. Size shall be equal to the sum of
 *             ncf. See pconn for explanation.
 * @param fneighbour Face to cell mapping. Its size shall be equal to
 *                   the number of faces times 2. For each face, the
 *                   two entries are either a cell number or -1
 *                   (signifying the outer boundary). The face normal
 *                   points out of the first cell and into the second.
 * @param fcentroid Face centroids. Size shall be equal to the number
 *                  of faces times d.
 * @param fnormal Face normale. Size shall be equal to the number
 *                of faces times d.
 * @param farea Face areas.
 * @param ccentroid Cell centroids. Size shall be ncells*d.
 * @param cvol Cell volumes.
 * @param perm Permeability. Size shall be ncells*d*d, storing a
 *             d-by-d positive definite tensor per cell.
 * @param[out] Binv This is where the inner product will be
 *                  stored. Its size shall be equal to \f$\sum_i
 * n_i^2\f$.
 */
void mim_ip_simple_all(int ncells, int d, int max_ncf, int *ncf,
                       int *pconn, int *conn,
                       int *fneighbour, double *fcentroid, double *fnormal,
                       double *farea, double *ccentroid, double *cvol,
                       double *perm, double *Binv);

#ifdef __cplusplus
}
#endif

#endif /* OPM_MIMETIC_HEADER_INCLUDED */
