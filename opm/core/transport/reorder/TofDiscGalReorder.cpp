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

#include <opm/core/grid/CellQuadrature.hpp>
#include <opm/core/grid/FaceQuadrature.hpp>
#include <opm/core/transport/reorder/TofDiscGalReorder.hpp>
#include <opm/core/transport/reorder/DGBasis.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/VelocityInterpolation.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/linalg/blas_lapack.h>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace Opm
{


    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    /// \param[in] param     Parameters for the solver.
    ///                      The following parameters are accepted (defaults):
    ///   use_cvi (false)                         Use ECVI velocity interpolation.
    ///   use_limiter (false)                     Use a slope limiter. If true, the next three parameters are used.
    ///   limiter_relative_flux_threshold (1e-3)  Ignore upstream fluxes below this threshold, relative to total cell flux.
    ///   limiter_method ("MinUpwindFace")        Limiter method used. Accepted methods are:
    ///                                             MinUpwindFace              Limit cell tof to >= inflow face tofs.
    ///   limiter_usage ("DuringComputations")    Usage pattern for limiter. Accepted choices are:
    ///                                             DuringComputations         Apply limiter to cells as they are computed,
    ///                                                                        so downstream cells' solutions may be affected
    ///                                                                        by limiting in upstream cells.
    ///                                             AsPostProcess              Apply in dependency order, but only after
    ///                                                                        computing (unlimited) solution.
    ///                                             AsSimultaneousPostProcess  Apply to each cell independently, using un-
    ///                                                                        limited solution in neighbouring cells.
    TofDiscGalReorder::TofDiscGalReorder(const UnstructuredGrid& grid,
                                         const parameter::ParameterGroup& param)
        : grid_(grid),
          use_cvi_(false),
          use_limiter_(false),
          limiter_relative_flux_threshold_(1e-3),
          limiter_method_(MinUpwindAverage),
          limiter_usage_(DuringComputations),
          coord_(grid.dimensions),
          velocity_(grid.dimensions)
    {
        const int dg_degree = param.getDefault("dg_degree", 0);
        const bool use_tensorial_basis = param.getDefault("use_tensorial_basis", false);
        if (use_tensorial_basis) {
            basis_func_.reset(new DGBasisMultilin(grid_, dg_degree));
        } else {
            basis_func_.reset(new DGBasisBoundedTotalDegree(grid_, dg_degree));
        }

        use_cvi_ = param.getDefault("use_cvi", use_cvi_);
        use_limiter_ = param.getDefault("use_limiter", use_limiter_);
        if (use_limiter_) {
            limiter_relative_flux_threshold_ = param.getDefault("limiter_relative_flux_threshold",
                                                                limiter_relative_flux_threshold_);
            const std::string limiter_method_str = param.getDefault<std::string>("limiter_method", "MinUpwindAverage");
            if (limiter_method_str == "MinUpwindFace") {
                limiter_method_ = MinUpwindFace;
            } else if (limiter_method_str == "MinUpwindAverage") {
                limiter_method_ = MinUpwindAverage;
            } else {
                THROW("Unknown limiter method: " << limiter_method_str);
            }
            const std::string limiter_usage_str = param.getDefault<std::string>("limiter_usage", "DuringComputations");
            if (limiter_usage_str == "DuringComputations") {
                limiter_usage_ = DuringComputations;
            } else if (limiter_usage_str == "AsPostProcess") {
                limiter_usage_ = AsPostProcess;
            } else if (limiter_usage_str == "AsSimultaneousPostProcess") {
                limiter_usage_ = AsSimultaneousPostProcess;
            } else {
                THROW("Unknown limiter usage spec: " << limiter_usage_str);
            }
        }
        // A note about the use_cvi_ member variable:
        // In principle, we should not need it, since the choice of velocity
        // interpolation is made below, but we may need to use higher order
        // quadrature to exploit CVI, so we store the choice.
        // An alternative would be to add a virtual method isConstant() to
        // the VelocityInterpolationInterface.
        if (use_cvi_) {
            velocity_interpolation_.reset(new VelocityInterpolationECVI(grid_));
        } else {
            velocity_interpolation_.reset(new VelocityInterpolationConstant(grid_));
        }
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
    void TofDiscGalReorder::solveTof(const double* darcyflux,
                                                  const double* porevolume,
                                                  const double* source,
                                                  std::vector<double>& tof_coeff)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        source_ = source;
#ifndef NDEBUG
        // Sanity check for sources.
        const double cum_src = std::accumulate(source, source + grid_.number_of_cells, 0.0);
        if (std::fabs(cum_src) > *std::max_element(source, source + grid_.number_of_cells)*1e-2) {
            // THROW("Sources do not sum to zero: " << cum_src);
            MESSAGE("Warning: sources do not sum to zero: " << cum_src);
        }
#endif
        const int num_basis = basis_func_->numBasisFunc();
        tof_coeff.resize(num_basis*grid_.number_of_cells);
        std::fill(tof_coeff.begin(), tof_coeff.end(), 0.0);
        tof_coeff_ = &tof_coeff[0];
        rhs_.resize(num_basis);
        jac_.resize(num_basis*num_basis);
        orig_jac_.resize(num_basis*num_basis);
        basis_.resize(num_basis);
        basis_nb_.resize(num_basis);
        grad_basis_.resize(num_basis*grid_.dimensions);
        velocity_interpolation_->setupFluxes(darcyflux);
        reorderAndTransport(grid_, darcyflux);
        switch (limiter_usage_) {
        case AsPostProcess:
            applyLimiterAsPostProcess();
            break;
        case AsSimultaneousPostProcess:
            applyLimiterAsSimultaneousPostProcess();
            break;
        case DuringComputations:
            // Do nothing.
            break;
        default:
            THROW("Unknown limiter usage choice: " << limiter_usage_);
        }
    }




    void TofDiscGalReorder::solveSingleCell(const int cell)
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
        const int num_basis = basis_func_->numBasisFunc();

        std::fill(rhs_.begin(), rhs_.end(), 0.0);
        std::fill(jac_.begin(), jac_.end(), 0.0);

        // Compute cell residual contribution.
        {
            const int deg_needed = basis_func_->degree();
            CellQuadrature quad(grid_, cell, deg_needed);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                // Integral of: b_i \phi
                quad.quadPtCoord(quad_pt, &coord_[0]);
                basis_func_->eval(cell, &coord_[0], &basis_[0]);
                const double w = quad.quadPtWeight(quad_pt);
                for (int j = 0; j < num_basis; ++j) {
                    rhs_[j] += w * basis_[j] * porevolume_[cell] / grid_.cell_volumes[cell];
                }
            }
        }

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
            if (flux >= 0.0) {
                // This is an outflow boundary.
                continue;
            }
            if (upstream_cell < 0) {
                // This is an outer boundary. Assumed tof = 0 on inflow, so no contribution.
                continue;
            }
            // Do quadrature over the face to compute
            // \int_{\partial K} u_h^{ext} (v(x) \cdot n) b_j ds
            // (where u_h^{ext} is the upstream unknown (tof)).
            // Quadrature degree set to 2*D, since u_h^{ext} varies
            // with degree D, and b_j too. We assume that the normal
            // velocity is constant (this assumption may have to go
            // for higher order than DG1).
            const double normal_velocity = flux / grid_.face_areas[face];
            const int deg_needed = 2*basis_func_->degree();
            FaceQuadrature quad(grid_, face, deg_needed);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                quad.quadPtCoord(quad_pt, &coord_[0]);
                basis_func_->eval(cell, &coord_[0], &basis_[0]);
                basis_func_->eval(upstream_cell, &coord_[0], &basis_nb_[0]);
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
            // Even with ECVI velocity interpolation, degree of precision 1
            // is sufficient for optimal convergence order for DG1 when we
            // use linear (total degree 1) basis functions.
            // With bi(tri)-linear basis functions, it still seems sufficient
            // for convergence order 2, but the solution looks much better and
            // has significantly lower error with degree of precision 2.
            // For now, we err on the side of caution, and use 2*degree, even
            // though this is wasteful for the pure linear basis functions.
            // const int deg_needed = 2*basis_func_->degree() - 1;
            const int deg_needed = 2*basis_func_->degree();
            CellQuadrature quad(grid_, cell, deg_needed);
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                // b_i (v \cdot \grad b_j)
                quad.quadPtCoord(quad_pt, &coord_[0]);
                basis_func_->eval(cell, &coord_[0], &basis_[0]);
                basis_func_->evalGrad(cell, &coord_[0], &grad_basis_[0]);
                velocity_interpolation_->interpolate(cell, &coord_[0], &velocity_[0]);
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
            FaceQuadrature quad(grid_, face, 2*basis_func_->degree());
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                // u^ext flux B   (B = {b_j})
                quad.quadPtCoord(quad_pt, &coord_[0]);
                basis_func_->eval(cell, &coord_[0], &basis_[0]);
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
            CellQuadrature quad(grid_, cell, 2*basis_func_->degree());
            for (int quad_pt = 0; quad_pt < quad.numQuadPts(); ++quad_pt) {
                quad.quadPtCoord(quad_pt, &coord_[0]);
                basis_func_->eval(cell, &coord_[0], &basis_[0]);
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
        orig_jac_ = jac_;
        orig_rhs_ = rhs_;
        dgesv_(&n, &nrhs, &jac_[0], &lda, &piv[0], &rhs_[0], &ldb, &info);
        if (info != 0) {
            // Print the local matrix and rhs.
            std::cerr << "Failed solving single-cell system Ax = b in cell " << cell
                      << " with A = \n";
            for (int row = 0; row < n; ++row) {
                for (int col = 0; col < n; ++col) {
                    std::cerr << "    " << orig_jac_[row + n*col];
                }
                std::cerr << '\n';
            }
            std::cerr << "and b = \n";
            for (int row = 0; row < n; ++row) {
                std::cerr << "    " << orig_rhs_[row] << '\n';
            }
            THROW("Lapack error: " << info << " encountered in cell " << cell);
        }

        // The solution ends up in rhs_, so we must copy it.
        std::copy(rhs_.begin(), rhs_.end(), tof_coeff_ + num_basis*cell);

        // Apply limiter.
        if (basis_func_->degree() > 0 && use_limiter_ && limiter_usage_ == DuringComputations) {
#ifdef EXTRA_VERBOSE
            std::cout << "Cell: " << cell << "   ";
            std::cout << "v = ";
            for (int dd = 0; dd < dim; ++dd) {
                std::cout << velocity_[dd] << ' ';
            }
            std::cout << "     grad tau = ";
            for (int dd = 0; dd < dim; ++dd) {
                std::cout << tof_coeff_[num_basis*cell + dd + 1] << ' ';
            }
            const double prod = std::inner_product(velocity_.begin(), velocity_.end(),
                                                   tof_coeff_ + num_basis*cell + 1, 0.0);
            const double vv = std::inner_product(velocity_.begin(), velocity_.end(),
                                                 velocity_.begin(), 0.0);
            const double gg = std::inner_product(tof_coeff_ + num_basis*cell + 1,
                                                 tof_coeff_ + num_basis*cell + num_basis,
                                                 tof_coeff_ + num_basis*cell + 1, 0.0);
            std::cout << "     prod = " << std::inner_product(velocity_.begin(), velocity_.end(),
                                                              tof_coeff_ + num_basis*cell + 1, 0.0);
            std::cout << "     normalized = " << prod/std::sqrt(vv*gg);
            std::cout << "     angle = " << std::acos(prod/std::sqrt(vv*gg))*360.0/(2.0*M_PI);
            std::cout << std::endl;
#endif
            applyLimiter(cell, tof_coeff_);
        }
    }




    void TofDiscGalReorder::solveMultiCell(const int num_cells, const int* cells)
    {
        std::cout << "Pretending to solve multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        for (int i = 0; i < num_cells; ++i) {
            solveSingleCell(cells[i]);
        }
    }




    void TransportModelTracerTofDiscGal::applyLimiter(const int cell, double* tof)
    {
        switch (limiter_method_) {
        case MinUpwindFace:
            applyMinUpwindLimiter(cell, true, tof);
            break;
        case MinUpwindAverage:
            applyMinUpwindLimiter(cell, false, tof);
            break;
        default:
            THROW("Limiter type not implemented: " << limiter_method_);
        }
    }




    void TransportModelTracerTofDiscGal::applyMinUpwindLimiter(const int cell, const bool face_min, double* tof)
    {
        if (basis_func_->degree() != 1) {
            THROW("This limiter only makes sense for our DG1 implementation.");
        }

        // Limiter principles:
        // 1. Let M be either:
        //      - the minimum TOF value of all upstream faces,
        //        evaluated in the upstream cells
        //        (chosen if face_min is true).
        //    or:
        //      - the minimum average TOF value of all upstream cells
        //        (chosen if face_min is false).
        //    Then the value at all points in this cell shall be at
        //    least M. Upstream faces whose flux does not exceed the
        //    relative flux threshold are not considered for this
        //    minimum.
        // 2. The TOF shall not be below zero in any point.

        // Find minimum tof on upstream faces/cells and for this cell.
        const int num_basis = basis_func_->numBasisFunc();
        double min_upstream_tof = 1e100;
        double min_here_tof = 1e100;
        int num_upstream_faces = 0;
        const double total_flux = totalFlux(cell);
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
            const bool upstream = (flux < -total_flux*limiter_relative_flux_threshold_);
            const bool interior = (upstream_cell >= 0);

            // Find minimum tof in this cell and upstream.
            // The meaning of minimum upstream tof depends on method.
            min_here_tof = std::min(min_here_tof, minCornerVal(cell, face));
            if (upstream) {
                ++num_upstream_faces;
                double upstream_tof = 0.0;
                if (interior) {
                    if (face_min) {
                        upstream_tof = minCornerVal(upstream_cell, face);
                    } else {
                        upstream_tof = basis_func_->functionAverage(tof_coeff_ + num_basis*upstream_cell);
                    }
                }
                min_upstream_tof = std::min(min_upstream_tof, upstream_tof);
            }
        }

        // Compute slope multiplier (limiter).
        if (num_upstream_faces == 0) {
            min_upstream_tof = 0.0;
            min_here_tof = 0.0;
        }
        if (min_upstream_tof < 0.0) {
            min_upstream_tof = 0.0;
        }
        const double tof_c = basis_func_->functionAverage(tof_coeff_ + num_basis*cell);
        double limiter = (tof_c - min_upstream_tof)/(tof_c - min_here_tof);
        if (tof_c < min_upstream_tof) {
            // Handle by setting a flat solution.
            std::cout << "Trouble in cell " << cell << std::endl;
            limiter = 0.0;
            basis_func_->addConstant(min_upstream_tof - tof_c, tof + num_basis*cell);
        }
        ASSERT(limiter >= 0.0);

        // Actually do the limiting (if applicable).
        if (limiter < 1.0) {
            // std::cout << "Applying limiter in cell " << cell << ", limiter = " << limiter << std::endl;
            basis_func_->multiplyGradient(limiter, tof + num_basis*cell);
        } else {
            // std::cout << "Not applying limiter in cell " << cell << "!" << std::endl;
        }
    }






    void TransportModelTracerTofDiscGal::applyLimiterAsPostProcess()
    {
        // Apply the limiter sequentially to all cells.
        // This means that a cell's limiting behaviour may be affected by
        // any limiting applied to its upstream cells.
        const std::vector<int>& seq = TransportModelInterface::sequence();
        const int nc = seq.size();
        ASSERT(nc == grid_.number_of_cells);
        for (int i = 0; i < nc; ++i) {
            const int cell = seq[i];
            applyLimiter(cell, tof_coeff_);
        }
    }




    void TransportModelTracerTofDiscGal::applyLimiterAsSimultaneousPostProcess()
    {
        // Apply the limiter simultaneously to all cells.
        // This means that each cell is limited independently from all other cells,
        // we write the resulting dofs to a new array instead of writing to tof_coeff_.
        // Afterwards we copy the results back to tof_coeff_.
        const int num_basis = basis_func_->numBasisFunc();
        std::vector<double> tof_coeffs_new(tof_coeff_, tof_coeff_ + num_basis*grid_.number_of_cells);
        for (int c = 0; c < grid_.number_of_cells; ++c) {
            applyLimiter(c, &tof_coeffs_new[0]);
        }
        std::copy(tof_coeffs_new.begin(), tof_coeffs_new.end(), tof_coeff_);
    }




    double TransportModelTracerTofDiscGal::totalFlux(const int cell) const
    {
        // Find total upstream/downstream fluxes.
        double upstream_flux = 0.0;
        double downstream_flux = 0.0;
        for (int hface = grid_.cell_facepos[cell]; hface < grid_.cell_facepos[cell+1]; ++hface) {
            const int face = grid_.cell_faces[hface];
            double flux = 0.0;
            if (cell == grid_.face_cells[2*face]) {
                flux = darcyflux_[face];
            } else {
                flux = -darcyflux_[face];
            }
            if (flux < 0.0) {
                upstream_flux += flux;
            } else {
                downstream_flux += flux;
            }
        }
        // In the presence of sources, significant fluxes may be missing from the computed fluxes,
        // setting the total flux to the (positive) maximum avoids this: since source is either
        // inflow or outflow, not both, either upstream_flux or downstream_flux must be correct.
        return std::max(-upstream_flux, downstream_flux);
    }




    double TransportModelTracerTofDiscGal::minCornerVal(const int cell, const int face) const
    {
        // Evaluate the solution in all corners.
        const int dim = grid_.dimensions;
        const int num_basis = basis_func_->numBasisFunc();
        double min_cornerval = 1e100;
        for (int fnode = grid_.face_nodepos[face]; fnode < grid_.face_nodepos[face+1]; ++fnode) {
            const double* nc = grid_.node_coordinates + dim*grid_.face_nodes[fnode];
            basis_func_->eval(cell, nc, &basis_[0]);
            const double tof_corner = std::inner_product(basis_.begin(), basis_.end(),
                                                         tof_coeff_ + num_basis*cell, 0.0);
            min_cornerval = std::min(min_cornerval, tof_corner);
        }
        return min_cornerval;
    }




} // namespace Opm
