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

/**
 * \file
 * Routines to assist mimetic discretisations of the flow equation.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Form linear operator to span the null space of the normal vectors
 * of a grid cell.
 *
 * Specifically,
 * \f[
 * \begin{aligned}
 * X &= \operatorname{diag}(A) (I - QQ^\mathsf{T})
 *      \operatorname{diag}(A), \\
 * Q &= \operatorname{orth}(\operatorname{diag}(A) C)
 * \end{aligned}
 * \f]
 * in which \f$\operatorname{orth}(M)\f$ denotes an orthonormal
 * basis for the colum space (range) of the matrix \f$M\f$,
 * represented as a matrix.
 *
 * @param[in]     nf    Number of faces connected to single grid cell.
 * @param[in]     nconn Total number of grid cell connections.
 *                      Typically equal to @c nf.
 * @param[in]     d     Number of physical dimensions.
 *                      Assumed less than four.
 * @param[in,out] C     Centroid vectors.  Specifically,
 *                      \f$c_{ij} = \Bar{x}_{ij} - \Bar{x}_{cj}\f$.
 *                      Array of size \f$\mathit{nf}\times d\f$
 *                      in column major (Fortran) order.
 *                      Contents destroyed on output.
 * @param[in]     A     Interface areas.
 * @param[out]    X     Null space linear operator.  Array of size
 *                      \f$\mathit{nconn}\times\mathit{nconn}\f$
 *                      in column major (Fortran) order.  On output,
 *                      the upper left \f$\mathit{nf}\times\mathit{nf}\f$
 *                      sub-matrix contains the required null space
 *                      linear operator.
 * @param[out]    work  Scratch array of size at least @c nconn.
 * @param[in]     lwork Actual size of scratch array.
 */
void
mim_ip_span_nullspace(int nf, int nconn, int d,
                      double *C,
                      double *A,
                      double *X,
                      double *work, int lwork);

/**
 * Form (inverse) mimetic inner product that reproduces linear
 * pressure drops (constant velocity) on general polyhedral cells.
 *
 * Specifically
 * \f[
 * B^{-1} = \frac{1}{v} \big(NKN^\mathsf{T} + \frac{6t}{d}\,X\big)
 * \f]
 * in which \f$t = \operatorname{tr}(K)\f$ is the trace of \f$K\f$
 * and \f$X\f$ is the result of function mim_ip_span_nullspace().
 *
 * @param[in]     nf    Number of faces connected to single grid cell.
 * @param[in]     nconn Total number of grid cell connections.
 *                      Typically equal to @c nf.
 * @param[in]     d     Number of physical dimensions.
 *                      Assumed less than four.
 * @param[in]     vol   Cell volume.
 * @param[in]     K     Permeability.  A \f$d\times d\f$ matrix in
 *                      column major (Fortran) order.
 * @param[in]     N     Normal vectors.  An \f$\mathit{nf}\times d\f$
 *                      matrix in column major (Fortran) order.
 * @param[in,out] Binv  Inverse inner product result.  An
 *                      \f$\mathit{nconn}\times\mathit{nconn}\f$
 *                      matrix in column major format.  On input,
 *                      the result of mim_ip_span_nullspace().  On
 *                      output, the upper left
 *                      \f$\mathit{nf}\times\mathit{nf}\f$ sub-matrix
 *                      will be overwritten with \f$B^{-1}\f$.
 * @param[in,out] work  Scratch array of size at least <CODE>nf * d</CODE>.
 * @param[in]     lwork Actual size of scratch array.
 */
void
mim_ip_linpress_exact(int nf, int nconn, int d,
                      double vol, double *K,
                      double *N,
                      double *Binv,
                      double *work, int lwork);


/**
 * Convenience wrapper around the function pair mim_ip_span_nullspace()
 * and mim_ip_linpress_exact().
 *
 * @param[in]     nf    Number of faces connected to single grid cell.
 * @param[in]     nconn Total number of grid cell connections.
 *                      Typically equal to @c nf.
 * @param[in]     d     Number of physical dimensions.
 *                      Assumed less than four.
 * @param[in]     v     Cell volume.
 * @param[in]     K     Permeability.  A \f$d\times d\f$ matrix in
 *                      column major (Fortran) order.
 * @param[in,out] C     Centroid vectors.  Specifically,
 *                      \f$c_{ij} = \Bar{x}_{ij} - \Bar{x}_{cj}\f$.
 *                      Array of size \f$\mathit{nf}\times d\f$
 *                      in column major (Fortran) order.
 *                      Contents destroyed on output.
 * @param[in]     A     Interface areas.
 * @param[in]     N     Outward normal vectors.
 *                      An \f$\mathit{nf}\times d\f$ matrix in
 *                      column major (Fortran) order.
 * @param[out]    Binv  Inverse inner product result.  An
 *                      \f$\mathit{nconn}\times\mathit{nconn}\f$
 *                      matrix in column major format.  On
 *                      output, the upper left
 *                      \f$\mathit{nf}\times\mathit{nf}\f$ sub-matrix
 *                      will be overwritten with \f$B^{-1}\f$
 *                      defined by function mim_ip_linpress_exact().
 * @param[in,out] work  Scratch array of size at least <CODE>nf * d</CODE>.
 * @param[in]     lwork Actual size of scratch array.
 */
void
mim_ip_simple(int nf, int nconn, int d,
              double v, double *K, double *C,
              double *A, double *N,
              double *Binv,
              double *work, int lwork);


/**
 * Compute the mimetic inner products given a grid and cell-wise
 * permeability tensors.
 *
 * This function applies mim_ip_simple() to all specified cells.
 *
 * @param[in]  ncells       Number of cells.
 * @param[in]  d            Number of physical dimensions.
 * @param[in]  max_ncf      Maximum number of connections (faces)
 *                          of any individual cell.
 * @param[in]  pconn        Start pointers of cell-to-face topology
 *                          mapping.
 * @param[in]  conn         Actual cell-to-face topology mapping.
 * @param[in]  fneighbour   Face-to-cell mapping.
 * @param[in]  fcentroid    Face centroids.
 * @param[in]  fnormal      Face normals.
 * @param[in]  farea        Face areas.
 * @param[in]  ccentroid    Cell centroids.
 * @param[in]  cvol         Cell volumes.
 * @param[in]  perm         Cell permeability.
 * @param[out] Binv         Inverse inner product result.  Must point
 *                          to an array of size at least
 *                          \f$\sum_c n_c^2\f$ when \f$n_c\f$ denotes
 *                          the number of connections (faces) of
 *                          cell \f$c\f$.
 */
void
mim_ip_simple_all(int ncells, int d, int max_ncf,
                  int *pconn, int *conn,
                  int *fneighbour, double *fcentroid, double *fnormal,
                  double *farea, double *ccentroid, double *cvol,
                  double *perm, double *Binv);

/**
 * Compute local, static gravity pressure contributions to Darcy
 * flow equation discretised using a mimetic finite-difference method.
 *
 * The pressure contribution of local face \f$i\f$ in cell \f$c\f$ is
 * \f[
 * \mathit{gpress}_{\mathit{pconn}_c + i} =
 * \vec{g}\cdot (\Bar{x}_{\mathit{conn}_{\mathit{pconn}_c + i}}
 * - \Bar{x}_c)
 * \f]
 *
 * @param[in]  nc        Number of cells.
 * @param[in]  d         Number of physcial dimensions.
 * @param[in]  grav      Gravity vector.  Array of size @c d.
 * @param[in]  pconn     Start pointers of cell-to-face topology
 *                       mapping.
 * @param[in]  conn      Actual cell-to-face topology mapping.
 * @param[in]  fcentroid Face centroids.
 * @param[in]  ccentroid Cell centroids.
 * @param[out] gpress    Gravity pressure result.  Array of size
 *                       at least <CODE>pconn[nc]</CODE>.
 */
void
mim_ip_compute_gpress(int nc, int d, const double *grav,
                      const int *pconn, const int *conn,
                      const double *fcentroid, const double *ccentroid,
                      double *gpress);


/**
 * Incorporate effects of multiple phases in mimetic discretisation of
 * flow equations.
 *
 * Specifically, update the (inverse) inner products \f$B^{-1}\f$
 * previously computed using function mim_ip_linpress_exact() according
 * to the rule
 * \f[
 * \Tilde{B}_c^{-1} = \frac{1}{\lambda_{T,c}} B_c^{-1},
 * \quad i=0,\dots,\mathit{nc}-1
 * \f]
 * in which \f$B_c^{-1}\f$ denotes the result of mim_ip_linpress_exact()
 * for cell \f$c\f$ and \f$\lambda_{T,c}\f$ denotes the total mobility
 * of cell \f$c\f$.
 *
 * @param[in]  nc     Number of cells.
 * @param[in]  pconn  Start pointers of cell-to-face topology
 *                    mapping.
 * @param[in]  totmob Total mobility for all cells.  Array of size @c nc.
 * @param[in]  Binv0  Inverse inner product results for all cells.
 * @param[out] Binv   Inverse inner product results incorporating
 *                    effects of multiple fluid phases.
 */
void
mim_ip_mobility_update(int nc, const int *pconn, const double *totmob,
                       const double *Binv0, double *Binv);


/**
 * Incorporate effects of multiple fluid phases into existing, local,
 * static mimetic discretisations of gravity pressure.
 *
 * Specifically, update the result of mim_ip_compute_gpress()
 * according to the rule
 * \f[
 * \Tilde{G}_{\mathit{pconn}_c + i} = \omega_c\cdot
 * G_{\mathit{pconn}_c + i}, \quad i=\mathit{pconn}_c, \dots,
 * \mathit{pconn}_{c+1}-1, \quad c=0,\dots,\mathit{nc}-1
 * \f]
 * in which \f$\omega_c = (\sum_\alpha \lambda_{\alpha,c}
 * \rho_\alpha)/\lambda_{T,c}\f$ and \f$\Tilde{G}\f$ denotes the result
 * of function mim_ip_compute_gpress().
 *
 * @param[in]  nc      Number of cells.
 * @param[in]  pconn   Start pointers of cell-to-face topology
 *                     mapping.
 * @param[in]  omega   Sum of phase densities weighted by
 *                     fractional flow.
 * @param[in]  gpress0 Result of mim_ip_compute_gpress().
 * @param[out] gpress  Gravity pressure incorporating effects
 *                     of multiple fluid phases.
 */
void
mim_ip_density_update(int nc, const int *pconn, const double *omega,
                      const double *gpress0, double *gpress);

#ifdef __cplusplus
}
#endif

#endif /* OPM_MIMETIC_HEADER_INCLUDED */
