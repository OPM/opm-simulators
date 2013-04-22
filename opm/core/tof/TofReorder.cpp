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

#include "config.h"
#include <opm/core/tof/TofReorder.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Opm
{


    /// Construct solver.
    /// \param[in] gri d      A 2d or 3d grid.
    /// \param[in] linsolver  Linear solver used for multi-cell blocks.
    /// \param[in] use_multidim_upwind  If true, use multidimensional tof upwinding.
    TofReorder::TofReorder(const UnstructuredGrid& grid,
                           const LinearSolverInterface& linsolver,
                           const bool use_multidim_upwind)
        : grid_(grid),
          linsolver_(linsolver),
          darcyflux_(0),
          porevolume_(0),
          source_(0),
          tof_(0),
          tracer_(0),
          num_tracers_(0),
          block_index_(grid.number_of_cells, OutsideBlock),
          use_multidim_upwind_(use_multidim_upwind)
    {
    }




    /// Solve for time-of-flight.
    /// \param[in]  darcyflux         Array of signed face fluxes.
    /// \param[in]  porevolume        Array of pore volumes.
    /// \param[in]  source            Source term. Sign convention is:
    ///                                 (+) inflow flux,
    ///                                 (-) outflow flux.
    /// \param[out] tof               Array of time-of-flight values.
    void TofReorder::solveTof(const double* darcyflux,
                                           const double* porevolume,
                                           const double* source,
                                           std::vector<double>& tof)
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
        tof.resize(grid_.number_of_cells);
        std::fill(tof.begin(), tof.end(), 0.0);
        tof_ = &tof[0];
        if (use_multidim_upwind_) {
            face_tof_.resize(grid_.number_of_faces);
            std::fill(face_tof_.begin(), face_tof_.end(), 0.0);
        }
        num_tracers_ = 0;
        num_multicell_ = 0;
        max_size_multicell_ = 0;
        reorderAndTransport(grid_, darcyflux);
        std::cout << num_multicell_ << " multicell blocks with max size "
                  << max_size_multicell_ << " cells." << std::endl;
    }




    /// Solve for time-of-flight and a number of tracers.
    /// One tracer will be used for each inflow flux specified in
    /// the source parameter.
    /// \param[in]  darcyflux         Array of signed face fluxes.
    /// \param[in]  porevolume        Array of pore volumes.
    /// \param[in]  source            Source term. Sign convention is:
    ///                                 (+) inflow flux,
    ///                                 (-) outflow flux.
    /// \param[out] tof               Array of time-of-flight values (1 per cell).
    /// \param[out] tracer            Array of tracer values (N per cell, where N is
    ///                               the number of cells c for which source[c] > 0.0).
    void TofReorder::solveTofTracer(const double* darcyflux,
                                    const double* porevolume,
                                    const double* source,
                                    std::vector<double>& tof,
                                    std::vector<double>& tracer)
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
        tof.resize(grid_.number_of_cells);
        std::fill(tof.begin(), tof.end(), 0.0);
        tof_ = &tof[0];
        // Find the tracer heads (injectors).
        std::vector<int> tracerheads;
        for (int c = 0; c < grid_.number_of_cells; ++c) {
            if (source[c] > 0.0) {
                tracerheads.push_back(c);
            }
        }
        num_tracers_ = tracerheads.size();
        tracer.resize(grid_.number_of_cells*num_tracers_);
        std::fill(tracer.begin(), tracer.end(), 0.0);
        for (int tr = 0; tr < num_tracers_; ++tr) {
            tracer[tracerheads[tr]*num_tracers_ + tr] = 1.0;
        }
        tracer_ = &tracer[0];
        if (use_multidim_upwind_) {
            face_tof_.resize(grid_.number_of_faces);
            std::fill(face_tof_.begin(), face_tof_.end(), 0.0);
            THROW("Multidimensional upwind not yet implemented for tracer.");
        }
        reorderAndTransport(grid_, darcyflux);
    }




    void TofReorder::solveSingleCell(const int cell)
    {
        if (use_multidim_upwind_) {
            solveSingleCellMultidimUpwind(cell);
            return;
        }
        // Compute flux terms.
        // Sources have zero tof, and therefore do not contribute
        // to upwind_term. Sinks on the other hand, must be added
        // to the downwind_flux (note sign change resulting from
        // different sign conventions: pos. source is injection,
        // pos. flux is outflow).
        double upwind_term = 0.0;
        double downwind_flux = std::max(-source_[cell], 0.0);
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            double flux;
            int other;
            // Compute cell flux
            if (cell == grid_.face_cells[2*f]) {
                flux  = darcyflux_[f];
                other = grid_.face_cells[2*f+1];
            } else {
                flux  =-darcyflux_[f];
                other = grid_.face_cells[2*f];
            }
            // Add flux to upwind_term or downwind_flux
            if (flux < 0.0) {
                // Using tof == 0 on inflow, so we only add a
                // nonzero contribution if we are on an internal
                // face.
                if (other != -1) {
                    upwind_term += flux*tof_[other];
                    for (int tr = 0; tr < num_tracers_; ++tr) {
                        tracer_[num_tracers_*cell + tr] += flux*tracer_[num_tracers_*other + tr];
                    }
                }
            } else {
                downwind_flux += flux;
            }
        }

        // Compute tof.
        tof_[cell] = (porevolume_[cell] - upwind_term)/downwind_flux;

        // Compute tracers (if any).
        // Do not change tracer solution in source cells.
        if (source_[cell] <= 0.0) {
            for (int tr = 0; tr < num_tracers_; ++tr) {
                tracer_[num_tracers_*cell + tr] *= -1.0/downwind_flux;
            }
        }
    }



    void TofReorder::assembleSingleCell(const int cell,
                                        std::vector<int>& local_column,
                                        std::vector<double>& local_coefficient,
                                        double& rhs)
    {
        // Compute flux terms.
        // Sources have zero tof, and therefore do not contribute
        // to upwind_term. Sinks on the other hand, must be added
        // to the downwind_flux (note sign change resulting from
        // different sign conventions: pos. source is injection,
        // pos. flux is outflow).
        double downwind_flux = std::max(-source_[cell], 0.0);
        double upwind_term = 0.0;
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            double flux;
            int other;
            // Compute cell flux
            if (cell == grid_.face_cells[2*f]) {
                flux  = darcyflux_[f];
                other = grid_.face_cells[2*f+1];
            } else {
                flux  =-darcyflux_[f];
                other = grid_.face_cells[2*f];
            }
            // Add flux to upwind_term or downwind_flux
            if (flux < 0.0) {
                // Using tof == 0 on inflow, so we only add a
                // nonzero contribution if we are on an internal
                // face.
                if (other != -1) {
                    if (block_index_[other] == OutsideBlock) {
                        upwind_term += flux*tof_[other];
                    } else {
                        local_column.push_back(block_index_[other]);
                        local_coefficient.push_back(flux);
                    }
                }
            } else {
                downwind_flux += flux;
            }
        }
        local_column.push_back(block_index_[cell]);
        local_coefficient.push_back(downwind_flux);
        rhs = porevolume_[cell] - upwind_term;
    }




    void TofReorder::solveSingleCellMultidimUpwind(const int cell)
    {
        // Compute flux terms.
        // Sources have zero tof, and therefore do not contribute
        // to upwind_term. Sinks on the other hand, must be added
        // to the downwind terms (note sign change resulting from
        // different sign conventions: pos. source is injection,
        // pos. flux is outflow).
        double upwind_term = 0.0;
        double downwind_term_cell_factor = std::max(-source_[cell], 0.0);
        double downwind_term_face = 0.0;
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            double flux;
            // Compute cell flux
            if (cell == grid_.face_cells[2*f]) {
                flux  = darcyflux_[f];
            } else {
                flux  =-darcyflux_[f];
            }
            // Add flux to upwind_term or downwind_term_[face|cell_factor].
            if (flux < 0.0) {
                upwind_term += flux*face_tof_[f];
            } else if (flux > 0.0) {
                double fterm, cterm_factor;
                multidimUpwindTerms(f, cell, fterm, cterm_factor);
                downwind_term_face += fterm*flux;
                downwind_term_cell_factor += cterm_factor*flux;
            }
        }

        // Compute tof for cell.
        tof_[cell] = (porevolume_[cell] - upwind_term - downwind_term_face)/downwind_term_cell_factor;        // }

        // Compute tof for downwind faces.
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            const double outflux_f = (grid_.face_cells[2*f] == cell) ? darcyflux_[f] : -darcyflux_[f];
            if (outflux_f > 0.0) {
                double fterm, cterm_factor;
                multidimUpwindTerms(f, cell, fterm, cterm_factor);
                face_tof_[f] = fterm + cterm_factor*tof_[cell];
            }
        }
    }




    void TofReorder::solveMultiCell(const int num_cells, const int* cells)
    {
#if 0
        std::cout << "Pretending to solve multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        for (int i = 0; i < num_cells; ++i) {
            solveSingleCell(cells[i]);
        }
#else
        // std::cout << "Solving multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        ++num_multicell_;
        max_size_multicell_ = std::max(max_size_multicell_, num_cells);
        // Give each cell a local index in this multi-cell block.
        for (int ci = 0; ci < num_cells; ++ci) {
            block_index_[cells[ci]] = ci;
        }

        // Create sparse matrix and rhs vector for block equation system.
        std::vector<int> ia(num_cells+1);
        std::vector<int> ja;
        std::vector<double> sa;
        ja.reserve(2*num_cells);
        sa.reserve(2*num_cells);
        std::vector<double> rhs(num_cells);
        std::vector< std::pair<int, double> > rowdata;
        for (int ci = 0; ci < num_cells; ++ci) {
            ia[ci] = ja.size();
            assembleSingleCell(cells[ci], ja, sa, rhs[ci]);
            // We standardize sparse row format: must sort row content by column index.
            int num_row_elem = ja.size() - ia[ci];
            rowdata.resize(num_row_elem);
            for (int i = 0; i < num_row_elem; ++i) {
                rowdata[i].first = ja[ia[ci] + i];
                rowdata[i].second = sa[ia[ci] + i];
            }
            std::sort(rowdata.begin(), rowdata.end());
            for (int i = 0; i < num_row_elem; ++i) {
                ja[ia[ci] + i] = rowdata[i].first;
                sa[ia[ci] + i] = rowdata[i].second;
            }
        }
        ia.back() = ja.size();
        ASSERT(ja.size() == sa.size());

        // Solve system.
        std::vector<double> tof_block(num_cells);
        LinearSolverInterface::LinearSolverReport rep
            = linsolver_.solve(num_cells, ja.size(), &ia[0], &ja[0], &sa[0], &rhs[0], &tof_block[0]);
        if (!rep.converged) {
            THROW("Multicell system with " << num_cells << " failed to converge.");
        }

        // Write to global tof vector, and reset block indices to make
        // it usable for next solveMultiCell() call.
        for (int ci = 0; ci < num_cells; ++ci) {
            tof_[cells[ci]] = tof_block[ci];
            block_index_[cells[ci]] = OutsideBlock;
        }
#endif
    }




    // Assumes that face_tof_[f] is known for all upstream faces f of upwind_cell.
    // Assumes that darcyflux_[face] is != 0.0.
    // This function returns factors to compute the tof for 'face':
    //   tof(face) = face_term + cell_term_factor*tof(upwind_cell).
    // It is not computed here, since these factors are needed to
    // compute the tof(upwind_cell) itself.
    void TofReorder::multidimUpwindTerms(const int face,
                                         const int upwind_cell,
                                         double& face_term,
                                         double& cell_term_factor) const
    {
        // Implements multidim upwind according to
        // "Multidimensional upstream weighting for multiphase transport on general grids"
        // by Keilegavlen, Kozdon, Mallison.
        // However, that article does not give a 3d extension other than noting that using
        // multidimensional upwinding in the XY-plane and not in the Z-direction may be
        // a good idea. We have here attempted some generalization, by looking at all
        // face-neighbours across edges as upwind candidates, and giving them all uniform weight.
        // This will over-weight the immediate upstream cell value in an extruded 2d grid with
        // one layer (top and bottom no-flow faces will enter the computation) compared to the
        // original 2d case. Improvements are welcome.
        // Note: Modified algorithm to consider faces that share even a single vertex with
        // the input face. This reduces the problem of non-edge-conformal grids, but does not
        // eliminate it entirely.

        // Identify the adjacent faces of the upwind cell.
        const int* face_nodes_beg = grid_.face_nodes + grid_.face_nodepos[face];
        const int* face_nodes_end = grid_.face_nodes + grid_.face_nodepos[face + 1];
        ASSERT(face_nodes_end - face_nodes_beg == 2 || grid_.dimensions != 2);
        adj_faces_.clear();
        for (int hf = grid_.cell_facepos[upwind_cell]; hf < grid_.cell_facepos[upwind_cell + 1]; ++hf) {
            const int f = grid_.cell_faces[hf];
            if (f != face) {
                const int* f_nodes_beg = grid_.face_nodes + grid_.face_nodepos[f];
                const int* f_nodes_end = grid_.face_nodes + grid_.face_nodepos[f + 1];
                // Find out how many vertices they have in common.
                // Using simple linear searches since sets are small.
                int num_common = 0;
                for (const int* f_iter = f_nodes_beg; f_iter < f_nodes_end; ++f_iter) {
                    num_common += std::count(face_nodes_beg, face_nodes_end, *f_iter);
                }
                // Before: neighbours over an edge (3d) or vertex (2d).
                // Now: neighbours across a vertex.
                // if (num_common == grid_.dimensions - 1) {
                if (num_common > 0) {
                    adj_faces_.push_back(f);
                }
            }
        }

        // Indentify adjacent faces with inflows, compute omega_star, omega,
        // add up contributions.
        const int num_adj = adj_faces_.size();
        // The assertion below only holds if the grid is edge-conformal.
        // No longer testing, since method no longer requires it.
        // ASSERT(num_adj == face_nodes_end - face_nodes_beg);
        const double flux_face = std::fabs(darcyflux_[face]);
        face_term = 0.0;
        cell_term_factor = 0.0;
        for (int ii = 0; ii < num_adj; ++ii) {
            const int f = adj_faces_[ii];
            const double influx_f = (grid_.face_cells[2*f] == upwind_cell) ? -darcyflux_[f] : darcyflux_[f];
            const double omega_star = influx_f/flux_face;
            // SPU
            // const double omega = 0.0;
            // TMU
            // const double omega = omega_star > 0.0 ? std::min(omega_star, 1.0) :  0.0;
            // SMU
            const double omega = omega_star > 0.0 ? omega_star/(1.0 + omega_star) : 0.0;
            face_term += omega*face_tof_[f];
            cell_term_factor += (1.0 - omega);
        }
        face_term /= double(num_adj);
        cell_term_factor /= double(num_adj);
    }



} // namespace Opm
