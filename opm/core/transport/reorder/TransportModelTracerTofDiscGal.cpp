/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/linalg/blas_lapack.h>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Opm
{


    // ---------------   Helpers for TransportModelTracerTofDiscGal ---------------



    /// A class providing discontinuous Galerkin basis functions.
    struct DGBasis
    {
        static int numBasisFunc(const int dimensions,
                                const int degree)
        {
            switch (dimensions) {
            case 1:
                return degree + 1;
            case 2:
                return (degree + 2)*(degree + 1)/2;
            case 3:
                return (degree + 3)*(degree + 2)*(degree + 1)/6;
            default:
                THROW("Dimensions must be 1, 2 or 3.");
            }
        }

        /// Evaluate all nonzero basis functions at x,
        /// writing to f_x. The array f_x must have
        /// size numBasisFunc(grid.dimensions, degree).
        ///
        /// The basis functions are the following
        ///     Degree 0: 1.
        ///     Degree 1: x - xc, y - yc, z - zc etc.
        /// Further degrees await development.
        static void eval(const UnstructuredGrid& grid,
                         const int cell,
                         const int degree,
                         const double* x,
                         double* f_x)
        {
            const int dim = grid.dimensions;
            const double* cc = grid.cell_centroids + dim*cell;
            // Note intentional fallthrough in this switch statement!
            switch (degree) {
            case 1:
                for (int ix = 0; ix < dim; ++ix) {
                    f_x[1 + ix] = x[ix] - cc[ix];
                }
            case 0:
                f_x[0] = 1;
                break;
            default:
                THROW("Maximum degree is 1 for now.");
            }
        }

        /// Evaluate gradients of all nonzero basis functions at x,
        /// writing to grad_f_x. The array grad_f_x must have size
        /// numBasisFunc(grid.dimensions, degree) * grid.dimensions.
        /// The <grid.dimensions> components of the first basis function
        /// gradient come before the components of the second etc.
        static void evalGrad(const UnstructuredGrid& grid,
                             const int /*cell*/,
                             const int degree,
                             const double* /*x*/,
                             double* grad_f_x)
        {
            const int dim = grid.dimensions;
            const int num_basis = numBasisFunc(dim, degree);
            std::fill(grad_f_x, grad_f_x + num_basis*dim, 0.0);
            if (degree > 1) {
                THROW("Maximum degree is 1 for now.");
            } else if (degree == 1) {
                for (int ix = 0; ix < dim; ++ix) {
                    grad_f_x[dim*(ix + 1) + ix] = 1.0;
                }
            }
        }
    };



    static void cross(const double* a, const double* b, double* res)
    {
        res[0] = a[1]*b[2] - a[2]*b[1];
        res[1] = a[2]*b[0] - a[0]*b[2];
        res[2] = a[0]*b[1] - a[1]*b[0];
    }




    static double triangleArea3d(const double* p0,
                                 const double* p1,
                                 const double* p2)
    {
        double a[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        double b[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        double cr[3];
        cross(a, b, cr);
        return 0.5*std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
    }




    /// Calculates the determinant of a 3 x 3 matrix, represented as
    /// three three-dimensional arrays.
    static double determinantOf(const double* a0,
                                const double* a1,
                                const double* a2)
    {
        return
            a0[0] * (a1[1] * a2[2] - a2[1] * a1[2]) -
            a0[1] * (a1[0] * a2[2] - a2[0] * a1[2]) +
            a0[2] * (a1[0] * a2[1] - a2[0] * a1[1]);
    }




    /// Computes the volume of a tetrahedron consisting of 4 vertices
    /// with 3-dimensional coordinates
    static double tetVolume(const double* p0,
                            const double* p1,
                            const double* p2,
                            const double* p3)
    {
        double a[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        double b[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        double c[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
        return std::fabs(determinantOf(a, b, c) / 6.0);
    }




    /// A class providing numerical quadrature for cells.
    /// In general: \int_{cell} g(x) dx = \sum_{i=0}^{n-1} w_i g(x_i).
    /// Note that this class does multiply weights by cell volume,
    /// so weights always sum to cell volume.
    /// Degree 1 method:
    ///     Midpoint (centroid) method.
    ///         n = 1, w_0 = cell volume, x_0 = cell centroid
    /// Degree 2 method:
    ///    Based on subdivision of each cell face into triangles
    ///    with the face centroid as a common vertex, and then
    ///    subdividing the cell into tetrahedra with the cell
    ///    centroid as a common vertex. Then apply the tetrahedron
    ///    rule with the following 4 nodes (uniform weights):
    ///        a = 0.138196601125010515179541316563436
    ///        x_i has all barycentric coordinates = a, except for
    ///            the i'th coordinate which is = 1 - 3a.
    ///    This rule is from http://nines.cs.kuleuven.be/ecf,
    ///    it is the second degree 2 4-point rule for tets,
    ///    referenced to Stroud(1971).
    ///    The tetrahedra are numbered T_{i,j}, and are given by the
    ///    cell centroid C, the face centroid FC_i, and two nodes
    ///    of face i: FN_{i,j}, FN_{i,j+1}.
    class CellQuadrature
    {
    public:
        CellQuadrature(const UnstructuredGrid& grid,
                       const int cell,
                       const int degree)
            : grid_(grid), cell_(cell), degree_(degree)
        {
            if (degree > 2) {
                THROW("CellQuadrature exact for polynomial degrees > 1 not implemented.");
            }
            if (degree == 2) {
                // Prepare subdivision.
            }
        }

        int numQuadPts() const
        {
            if (degree_ < 2) {
                return 1;
            }
            // Degree 2 case.
            int sumnodes = 0;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                sumnodes += grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
            }
            return 4*sumnodes;
        }

        void quadPtCoord(const int index, double* coord) const
        {
            const int dim = grid_.dimensions;
            const double* cc = grid_.cell_centroids + dim*cell_;
            if (degree_ < 2) {
                std::copy(cc, cc + dim, coord);
                return;
            }
            // Degree 2 case.
            int tetindex = index / 4;
            const int subindex = index % 4;
            const double* nc = grid_.node_coordinates;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
                if (nfn <= tetindex) {
                    // Our tet is not associated with this face.
                    tetindex -= nfn;
                    continue;
                }
                const double* fc = grid_.face_centroids + dim*face;
                const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
                const int node0 = fnodes[tetindex];
                const int node1 = fnodes[(tetindex + 1) % nfn];
                const double* n0c = nc + dim*node0;
                const double* n1c = nc + dim*node1;
                const double a = 0.138196601125010515179541316563436;
                // Barycentric coordinates of our point in the tet.
                double baryc[4] = { a, a, a, a };
                baryc[subindex] = 1.0 - 3.0*a;
                for (int dd = 0; dd < dim; ++dd) {
                    coord[dd] = baryc[0]*cc[dd] + baryc[1]*fc[dd] + baryc[2]*n0c[dd] + baryc[3]*n1c[dd];
                }
                return;
            }
            THROW("Should never reach this point.");
        }

        double quadPtWeight(const int index) const
        {
            if (degree_ < 2) {
                return grid_.cell_volumes[cell_];
            }
            // Degree 2 case.
            const int dim = grid_.dimensions;
            const double* cc = grid_.cell_centroids + dim*cell_;
            int tetindex = index / 4;
            const double* nc = grid_.node_coordinates;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
                if (nfn <= tetindex) {
                    // Our tet is not associated with this face.
                    tetindex -= nfn;
                    continue;
                }
                const double* fc = grid_.face_centroids + dim*face;
                const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
                const int node0 = fnodes[tetindex];
                const int node1 = fnodes[(tetindex + 1) % nfn];
                const double* n0c = nc + dim*node0;
                const double* n1c = nc + dim*node1;
                return 0.25*tetVolume(cc, fc, n0c, n1c);
            }
            THROW("Should never reach this point.");
        }

    private:
        const UnstructuredGrid& grid_;
        const int cell_;
        const int degree_;
    };





    /// A class providing numerical quadrature for faces.
    /// In general: \int_{face} g(x) dx = \sum_{i=0}^{n-1} w_i g(x_i).
    /// Note that this class does multiply weights by face area,
    /// so weights always sum to face area.
    /// Degree 1 method:
    ///     Midpoint (centroid) method.
    ///         n = 1, w_0 = face area, x_0 = face centroid
    /// Degree 2 method:
    ///    Based on subdivision of the face into triangles,
    ///    with the centroid as a common vertex, and the triangle
    ///    edge midpoint rule.
    ///    Triangle i consists of the centroid C, nodes N_i and N_{i+1}.
    ///    Its area is A_i.
    ///        n = 2 * nn  (nn = num nodes in face)
    ///        For i = 0..(nn-1):
    ///        w_i      = 1/3 A_i.
    ///        w_{nn+i} = 1/3 A_{i-1} + 1/3 A_i
    ///        x_i      = (N_i + N_{i+1})/2
    ///        x_{nn+i} = (C + N_i)/2
    ///    All N and A indices are interpreted cyclic, modulus nn.
    class FaceQuadrature
    {
    public:
        FaceQuadrature(const UnstructuredGrid& grid,
                       const int face,
                       const int degree)
            : grid_(grid), face_(face), degree_(degree)
        {
            if (grid_.dimensions != 3) {
                THROW("FaceQuadrature only implemented for 3D case.");
            }
            if (degree_ > 2) {
                THROW("FaceQuadrature exact for polynomial degrees > 2 not implemented.");
            }
        }

        int numQuadPts() const
        {
            if (degree_ < 2) {
                return 1;
            }
            // Degree 2 case.
            return 2 * (grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_]);
        }

        void quadPtCoord(const int index, double* coord) const
        {
            const int dim = grid_.dimensions;
            const double* fc = grid_.face_centroids + dim*face_;
            if (degree_ < 2) {
                std::copy(fc, fc + dim, coord);
                return;
            }
            // Degree 2 case.
            const int nn = grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_];
            const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face_];
            const double* nc = grid_.node_coordinates;
            if (index < nn) {
                // Boundary edge midpoint.
                const int node0 = fnodes[index];
                const int node1 = fnodes[(index + 1)%nn];
                for (int dd = 0; dd < dim; ++dd) {
                    coord[dd] = 0.5*(nc[dim*node0 + dd] + nc[dim*node1 + dd]);
                }
            } else {
                // Interiour edge midpoint.
                // Recall that index is now in [nn, 2*nn).
                const int node = fnodes[index - nn];
                for (int dd = 0; dd < dim; ++dd) {
                    coord[dd] = 0.5*(nc[dim*node + dd] + fc[dd]);
                }
            }
        }

        double quadPtWeight(const int index) const
        {
            if (degree_ < 2) {
                return grid_.face_areas[face_];
            }
            // Degree 2 case.
            const int dim = grid_.dimensions;
            const double* fc = grid_.face_centroids + dim*face_;
            const int nn = grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_];
            const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face_];
            const double* nc = grid_.node_coordinates;
            if (index < nn) {
                // Boundary edge midpoint.
                const int node0 = fnodes[index];
                const int node1 = fnodes[(index + 1)%nn];
                const double area = triangleArea3d(nc + dim*node1, nc + dim*node0, fc);
                return area / 3.0;
            } else {
                // Interiour edge midpoint.
                // Recall that index is now in [nn, 2*nn).
                const int node0 = fnodes[(index - 1) % nn];
                const int node1 = fnodes[index - nn];
                const int node2 = fnodes[(index + 1) % nn];
                const double area0 = triangleArea3d(nc + dim*node1, nc + dim*node0, fc);
                const double area1 = triangleArea3d(nc + dim*node2, nc + dim*node1, fc);
                return (area0 + area1) / 3.0;
            }
        }

    private:
        const UnstructuredGrid& grid_;
        const int face_;
        const int degree_;
    };




    // Initial version: only a constant interpolation.
    static void interpolateVelocity(const UnstructuredGrid& grid,
                                    const int cell,
                                    const double* darcyflux,
                                    const double* /*x*/,
                                    double* v)
    {
        const int dim = grid.dimensions;
        std::fill(v, v + dim, 0.0);
        const double* cc = grid.cell_centroids + cell*dim;
        for (int hface = grid.cell_facepos[cell]; hface < grid.cell_facepos[cell+1]; ++hface) {
            const int face = grid.cell_faces[hface];
            const double* fc = grid.face_centroids + face*dim;
            double flux = 0.0;
            if (cell == grid.face_cells[2*face]) {
                flux = darcyflux[face];
            } else {
                flux = -darcyflux[face];
            }
            for (int dd = 0; dd < dim; ++dd) {
                v[dd] += flux * (fc[dd] - cc[dd]) / grid.cell_volumes[cell];
            }
        }
    }



    // ---------------   Methods of TransportModelTracerTofDiscGal ---------------



    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    TransportModelTracerTofDiscGal::TransportModelTracerTofDiscGal(const UnstructuredGrid& grid)
        : grid_(grid),
          coord_(grid.dimensions),
          velocity_(grid.dimensions)
    {
    }




    /// Solve for time-of-flight.
    /// \param[in]  darcyflux         Array of signed face fluxes.
    /// \param[in]  porevolume        Array of pore volumes.
    /// \param[in]  source            Source term. Sign convention is:
    ///                                 (+) inflow flux,
    ///                                 (-) outflow flux.
    /// \param[in]  degree            Polynomial degree of DG basis functions used.
    /// \param[out] tof_coeff         Array of time-of-flight solution coefficients.
    ///                               The values are ordered by cell, meaning that
    ///                               the K coefficients corresponding to the first
    ///                               cell comes before the K coefficients corresponding
    ///                               to the second cell etc.
    ///                               K depends on degree and grid dimension.
    void TransportModelTracerTofDiscGal::solveTof(const double* darcyflux,
                                                  const double* porevolume,
                                                  const double* source,
                                                  const int degree,
                                                  std::vector<double>& tof_coeff)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        source_ = source;
#ifndef NDEBUG
        // Sanity check for sources.
        const double cum_src = std::accumulate(source, source + grid_.number_of_cells, 0.0);
        if (std::fabs(cum_src) > *std::max_element(source, source + grid_.number_of_cells)*1e-2) {
            THROW("Sources do not sum to zero: " << cum_src);
        }
#endif
        degree_ = degree;
        const int num_basis = DGBasis::numBasisFunc(grid_.dimensions, degree_);
        tof_coeff.resize(num_basis*grid_.number_of_cells);
        std::fill(tof_coeff.begin(), tof_coeff.end(), 0.0);
        tof_coeff_ = &tof_coeff[0];
        rhs_.resize(num_basis);
        jac_.resize(num_basis*num_basis);
        basis_.resize(num_basis);
        basis_nb_.resize(num_basis);
        grad_basis_.resize(num_basis*grid_.dimensions);
        reorderAndTransport(grid_, darcyflux);
    }




    void TransportModelTracerTofDiscGal::solveSingleCell(const int cell)
    {
        // Residual:
        // For each cell K, basis function b_j (spanning V_h),
        // writing the solution u_h|K = \sum_i c_i b_i
        //  Res = - \int_K \sum_i c_i b_i v(x) \cdot \grad b_j dx
        //        + \int_{\partial K} F(u_h, u_h^{ext}, v(x) \cdot n) b_j ds
        //        - \int_K \phi b_j
        // This is linear in c_i, so we do not need any nonlinear iterations.
        // We assemble the jacobian and the right-hand side. The residual is
        // equal to Res = Jac*c - rhs, and we compute rhs directly.

        const int dim = grid_.dimensions;
        const int num_basis = DGBasis::numBasisFunc(dim, degree_);

        std::fill(rhs_.begin(), rhs_.end(), 0.0);
        std::fill(jac_.begin(), jac_.end(), 0.0);

        // Compute cell residual contribution.
        // Note: Assumes that \int_K b_j = 0 for all j > 0
        rhs_[0] += porevolume_[cell];

        // Compute upstream residual contribution.
        for (int hface = grid_.cell_facepos[cell]; hface < grid_.cell_facepos[cell+1]; ++hface) {
            const int face = grid_.cell_faces[hface];
            double flux = 0.0;
            int upstream_cell = -1;
            if (cell == grid_.face_cells[2*face]) {
                flux = darcyflux_[face];
                upstream_cell = grid_.face_cells[2*face+1];
            } else {
                flux = -darcyflux_[face];
                upstream_cell = grid_.face_cells[2*face];
            }
            if (upstream_cell < 0) {
                // This is an outer boundary. Assumed tof = 0 on inflow, so no contribution.
                continue;
            }
            if (flux >= 0.0) {
                // This is an outflow boundary.
                continue;
            }
            // Do quadrature over the face to compute
            // \int_{\partial K} u_h^{ext} (v(x) \cdot n) b_j ds
            // (where u_h^{ext} is the upstream unknown (tof)).
            const double normal_velocity = flux / grid_.face_areas[face];
            FaceQuadrature quad(grid_, face, degree_);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                quad.quadPtCoord(quad_pt, &coord_[0]);
                DGBasis::eval(grid_, cell, degree_, &coord_[0], &basis_[0]);
                DGBasis::eval(grid_, upstream_cell, degree_, &coord_[0], &basis_nb_[0]);
                const double tof_upstream = std::inner_product(basis_nb_.begin(), basis_nb_.end(),
                                                               tof_coeff_ + num_basis*upstream_cell, 0.0);
                const double w = quad.quadPtWeight(quad_pt);
                for (int j = 0; j < num_basis; ++j) {
                    rhs_[j] -= w * tof_upstream * normal_velocity * basis_[j];
                }
            }
        }

        // Compute cell jacobian contribution. We use Fortran ordering
        // for jac_, i.e. rows cycling fastest.
        {
            CellQuadrature quad(grid_, cell, 2*degree_ - 1);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                // b_i (v \cdot \grad b_j)
                quad.quadPtCoord(quad_pt, &coord_[0]);
                DGBasis::eval(grid_, cell, degree_, &coord_[0], &basis_[0]);
                DGBasis::evalGrad(grid_, cell, degree_, &coord_[0], &grad_basis_[0]);
                interpolateVelocity(grid_, cell, darcyflux_, &coord_[0], &velocity_[0]);
                const double w = quad.quadPtWeight(quad_pt);
                for (int j = 0; j < num_basis; ++j) {
                    for (int i = 0; i < num_basis; ++i) {
                        for (int dd = 0; dd < dim; ++dd) {
                            jac_[j*num_basis + i] -= w * basis_[j] * grad_basis_[dim*i + dd] * velocity_[dd];
                        }
                    }
                }
            }
        }

        // Compute downstream jacobian contribution from faces.
        for (int hface = grid_.cell_facepos[cell]; hface < grid_.cell_facepos[cell+1]; ++hface) {
            const int face = grid_.cell_faces[hface];
            double flux = 0.0;
            if (cell == grid_.face_cells[2*face]) {
                flux = darcyflux_[face];
            } else {
                flux = -darcyflux_[face];
            }
            if (flux <= 0.0) {
                // This is an inflow boundary.
                continue;
            }
            // Do quadrature over the face to compute
            // \int_{\partial K} b_i (v(x) \cdot n) b_j ds
            const double normal_velocity = flux / grid_.face_areas[face];
            FaceQuadrature quad(grid_, face, 2*degree_);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                // u^ext flux B   (B = {b_j})
                quad.quadPtCoord(quad_pt, &coord_[0]);
                DGBasis::eval(grid_, cell, degree_, &coord_[0], &basis_[0]);
                const double w = quad.quadPtWeight(quad_pt);
                for (int j = 0; j < num_basis; ++j) {
                    for (int i = 0; i < num_basis; ++i) {
                        jac_[j*num_basis + i] += w * basis_[i] * normal_velocity * basis_[j];
                    }
                }
            }
        }

        // Compute downstream jacobian contribution from sink terms.
        // Contribution from inflow sources would be
        // similar to the contribution from upstream faces, but
        // it is zero since we let all external inflow be associated
        // with a zero tof.
        if (source_[cell] < 0.0) {
            // A sink.
            const double flux = -source_[cell]; // Sign convention for flux: outflux > 0.
            const double flux_density = flux / grid_.cell_volumes[cell];
            // Do quadrature over the cell to compute
            // \int_{K} b_i flux b_j dx
            CellQuadrature quad(grid_, cell, 2*degree_);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                quad.quadPtCoord(quad_pt, &coord_[0]);
                DGBasis::eval(grid_, cell, degree_, &coord_[0], &basis_[0]);
                const double w = quad.quadPtWeight(quad_pt);
                for (int j = 0; j < num_basis; ++j) {
                    for (int i = 0; i < num_basis; ++i) {
                        jac_[j*num_basis + i] += w * basis_[i] * flux_density * basis_[j];
                    }
                }
            }
        }

        // Solve linear equation.
        MAT_SIZE_T n = num_basis;
        MAT_SIZE_T nrhs = 1;
        MAT_SIZE_T lda = num_basis;
        std::vector<MAT_SIZE_T> piv(num_basis);
        MAT_SIZE_T ldb = num_basis;
        MAT_SIZE_T info = 0;
        dgesv_(&n, &nrhs, &jac_[0], &lda, &piv[0], &rhs_[0], &ldb, &info);
        if (info != 0) {
            // Print the local matrix and rhs.
            std::cerr << "Failed solving single-cell system Ax = b in cell " << cell
                      << " with A = \n";
            for (int row = 0; row < n; ++row) {
                for (int col = 0; col < n; ++col) {
                    std::cerr << "    " << jac_[row + n*col];
                }
                std::cerr << '\n';
            }
            std::cerr << "and b = \n";
            for (int row = 0; row < n; ++row) {
                std::cerr << "    " << rhs_[row] << '\n';
            }
            THROW("Lapack error: " << info << " encountered in cell " << cell);
        }
        // The solution ends up in rhs_, so we must copy it.
        std::copy(rhs_.begin(), rhs_.end(), tof_coeff_ + num_basis*cell);
    }




    void TransportModelTracerTofDiscGal::solveMultiCell(const int num_cells, const int* cells)
    {
        std::cout << "Pretending to solve multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        for (int i = 0; i < num_cells; ++i) {
            solveSingleCell(cells[i]);
        }
    }




} // namespace Opm
