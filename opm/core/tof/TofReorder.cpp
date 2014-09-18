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
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SparseTable.hpp>

#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

namespace Opm
{


    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    /// \param[in] use_multidim_upwind  If true, use multidimensional tof upwinding.
    TofReorder::TofReorder(const UnstructuredGrid& grid,
                           const bool use_multidim_upwind)
        : grid_(grid),
          darcyflux_(0),
          porevolume_(0),
          source_(0),
          tof_(0),
          gauss_seidel_tol_(1e-3),
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
            // OPM_THROW(std::runtime_error, "Sources do not sum to zero: " << cum_src);
            OPM_MESSAGE("Warning: sources do not sum to zero: " << cum_src);
        }
#endif
        tof.resize(grid_.number_of_cells);
        std::fill(tof.begin(), tof.end(), 0.0);
        tof_ = &tof[0];
        if (use_multidim_upwind_) {
            face_tof_.resize(grid_.number_of_faces);
            std::fill(face_tof_.begin(), face_tof_.end(), 0.0);
            face_part_tof_.resize(grid_.face_nodepos[grid_.number_of_faces]);
            std::fill(face_part_tof_.begin(), face_part_tof_.end(), 0.0);
        }
        compute_tracer_ = false;
        executeSolve();
    }




    /// Solve for time-of-flight and a number of tracers.
    /// \param[in]  darcyflux         Array of signed face fluxes.
    /// \param[in]  porevolume        Array of pore volumes.
    /// \param[in]  source            Source term. Sign convention is:
    ///                                 (+) inflow flux,
    ///                                 (-) outflow flux.
    /// \param[in]  tracerheads       Table containing one row per tracer, and each
    ///                               row contains the source cells for that tracer.
    /// \param[out] tof               Array of time-of-flight values (1 per cell).
    /// \param[out] tracer            Array of tracer values. N per cell, where N is
    ///                               equalt to tracerheads.size().
    void TofReorder::solveTofTracer(const double* darcyflux,
                                    const double* porevolume,
                                    const double* source,
                                    const SparseTable<int>& tracerheads,
                                    std::vector<double>& tof,
                                    std::vector<double>& tracer)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        source_ = source;
        const int num_cells = grid_.number_of_cells;
#ifndef NDEBUG
        // Sanity check for sources.
        const double cum_src = std::accumulate(source, source + num_cells, 0.0);
        if (std::fabs(cum_src) > *std::max_element(source, source + num_cells)*1e-2) {
            OPM_THROW(std::runtime_error, "Sources do not sum to zero: " << cum_src);
        }
#endif
        tof.resize(num_cells);
        std::fill(tof.begin(), tof.end(), 0.0);
        tof_ = &tof[0];

        if (use_multidim_upwind_) {
            face_tof_.resize(grid_.number_of_faces);
            std::fill(face_tof_.begin(), face_tof_.end(), 0.0);
            face_part_tof_.resize(grid_.face_nodepos[grid_.number_of_faces]);
            std::fill(face_part_tof_.begin(), face_part_tof_.end(), 0.0);
        }

        // Execute solve for tof
        compute_tracer_ = false;
        executeSolve();

        // Find the tracer heads (injectors).
        const int num_tracers = tracerheads.size();
        tracer.resize(num_cells*num_tracers);
        std::fill(tracer.begin(), tracer.end(), 0.0);
        if (num_tracers > 0) {
            tracerhead_by_cell_.clear();
            tracerhead_by_cell_.resize(num_cells, NoTracerHead);
        }
        for (int tr = 0; tr < num_tracers; ++tr) {
            for (int i = 0; i < tracerheads[tr].size(); ++i) {
                const int cell = tracerheads[tr][i];
                tracer[num_cells * tr + cell] = 1.0;
                tracerhead_by_cell_[cell] = tr;
            }
        }

        // Execute solve for tracers.
        std::vector<double> fake_pv(num_cells, 0.0);
        porevolume_ = fake_pv.data();
        for (int tr = 0; tr < num_tracers; ++tr) {
            tof_ = tracer.data() + tr * num_cells;
            compute_tracer_ = true;
            executeSolve();
        }

        // Write output tracer data (transposing the computed data).
        std::vector<double> computed = tracer;
        for (int cell = 0; cell < num_cells; ++cell) {
            for (int tr = 0; tr < num_tracers; ++tr) {
                tracer[num_tracers * cell + tr] = computed[num_cells * tr + cell];
            }
        }
    }




    void TofReorder::executeSolve()
    {
        num_multicell_ = 0;
        max_size_multicell_ = 0;
        max_iter_multicell_ = 0;
        reorderAndTransport(grid_, darcyflux_);
        if (num_multicell_ > 0) {
            std::cout << num_multicell_ << " multicell blocks with max size "
                      << max_size_multicell_ << " cells in upto "
                      << max_iter_multicell_ << " iterations." << std::endl;
        }
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
        if (compute_tracer_ && tracerhead_by_cell_[cell] != NoTracerHead) {
            // This is a tracer head cell, already has solution.
            return;
        }
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
                }
            } else {
                downwind_flux += flux;
            }
        }

        // Compute tof.
        tof_[cell] = (porevolume_[cell] - upwind_term)/downwind_flux;
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
        if (compute_tracer_ && tracerhead_by_cell_[cell] != NoTracerHead) {
            // Do nothing to the value in this cell, since we are at a tracer head.
        } else {
            tof_[cell] = (porevolume_[cell] - upwind_term - downwind_term_face)/downwind_term_cell_factor;
        }

        // Compute tof for downwind faces.
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            const double outflux_f = (grid_.face_cells[2*f] == cell) ? darcyflux_[f] : -darcyflux_[f];
            if (outflux_f > 0.0) {
                double fterm, cterm_factor;
                multidimUpwindTerms(f, cell, fterm, cterm_factor);
                face_tof_[f] = fterm + cterm_factor*tof_[cell];

                // Combine locally computed (for each adjacent vertex) terms, with uniform weighting.
                const int* face_nodes_beg = grid_.face_nodes + grid_.face_nodepos[f];
                const int* face_nodes_end = grid_.face_nodes + grid_.face_nodepos[f + 1];
                assert((face_nodes_end - face_nodes_beg) == 2 || grid_.dimensions != 2);
                for (const int* fn_iter = face_nodes_beg; fn_iter < face_nodes_end; ++fn_iter) {
                    double loc_face_term = 0.0;
                    double loc_cell_term_factor = 0.0;
                    const int node_pos = fn_iter - grid_.face_nodes;
                    localMultidimUpwindTerms(f, cell, node_pos,
                                             loc_face_term, loc_cell_term_factor);
                    face_part_tof_[node_pos] = loc_face_term + loc_cell_term_factor * tof_[cell];
                }
            }
        }
    }




    void TofReorder::solveMultiCell(const int num_cells, const int* cells)
    {
        ++num_multicell_;
        max_size_multicell_ = std::max(max_size_multicell_, num_cells);
        // std::cout << "Multiblock solve with " << num_cells << " cells." << std::endl;

        // Using a Gauss-Seidel approach.
        double max_delta = 1e100;
        int num_iter = 0;
        while (max_delta > gauss_seidel_tol_) {
            max_delta = 0.0;
            ++num_iter;
            for (int ci = 0; ci < num_cells; ++ci) {
                const int cell = cells[ci];
                const double tof_before = tof_[cell];
                solveSingleCell(cell);
                max_delta = std::max(max_delta, std::fabs(tof_[cell] - tof_before));
            }
            // std::cout << "Max delta = " << max_delta << std::endl;
        }
        max_iter_multicell_ = std::max(max_iter_multicell_, num_iter);
    }




    // Assumes that face_part_tof_[node_pos] is known for all inflow
    // faces to 'upwind_cell' sharing vertices with 'face'. The index
    // 'node_pos' is the same as the one used for the grid face-node
    // connectivity.
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
        // Implements multidim upwind inspired by
        // "Multidimensional upstream weighting for multiphase transport on general grids"
        // by Keilegavlen, Kozdon, Mallison.
        // However, that article does not give a 3d extension other than noting that using
        // multidimensional upwinding in the XY-plane and not in the Z-direction may be
        // a good idea. We have here attempted some generalization, by treating each face-part
        // (association of a face and a vertex) as possibly influencing all downwind face-parts
        // of the neighbouring cell that share the same vertex.
        // The current implementation aims to reproduce 2d results for extruded 3d grids.

        // Combine locally computed (for each adjacent vertex) terms, with uniform weighting.
        const int* face_nodes_beg = grid_.face_nodes + grid_.face_nodepos[face];
        const int* face_nodes_end = grid_.face_nodes + grid_.face_nodepos[face + 1];
        const int num_terms = face_nodes_end - face_nodes_beg;
        assert(num_terms == 2 || grid_.dimensions != 2);
        face_term = 0.0;
        cell_term_factor = 0.0;
        for (const int* fn_iter = face_nodes_beg; fn_iter < face_nodes_end; ++fn_iter) {
            double loc_face_term = 0.0;
            double loc_cell_term_factor = 0.0;
            localMultidimUpwindTerms(face, upwind_cell, fn_iter - grid_.face_nodes,
                                     loc_face_term, loc_cell_term_factor);
            face_term += loc_face_term;
            cell_term_factor += loc_cell_term_factor;
        }
        face_term /= double(num_terms);
        cell_term_factor /= double(num_terms);

    }




    namespace {
        double weightFunc(const double w)
        {
            // SPU
            // return 0.0;
            // TMU
            return w > 0.0 ? std::min(w, 1.0) :  0.0;
            // SMU
            // return w > 0.0 ? w/(1.0 + w) : 0.0;
        }
    }




    void TofReorder::localMultidimUpwindTerms(const int face,
                                              const int upwind_cell,
                                              const int node_pos,
                                              double& face_term,
                                              double& cell_term_factor) const
    {
        // Loop over all faces adjacent to the given cell and the
        // vertex in position node_pos.
        // If that part's influx is positive, we store it, and also its associated
        // node position.
        std::vector<double> influx;
        std::vector<int> node_pos_influx;
        influx.reserve(5);
        node_pos_influx.reserve(5);
        const int node = grid_.face_nodes[node_pos];
        for (int hf = grid_.cell_facepos[upwind_cell]; hf < grid_.cell_facepos[upwind_cell + 1]; ++hf) {
            const int f = grid_.cell_faces[hf];
            if (f != face) {
                // Find out if the face 'f' is adjacent to vertex 'node'.
                const int* f_nodes_beg = grid_.face_nodes + grid_.face_nodepos[f];
                const int* f_nodes_end = grid_.face_nodes + grid_.face_nodepos[f + 1];
                const int* pos = std::find(f_nodes_beg, f_nodes_end, node);
                const int node_pos2 = pos - grid_.face_nodes;
                const bool is_adj = (pos != f_nodes_end);
                if (is_adj) {
                    const int num_parts = f_nodes_end - f_nodes_beg;
                    const double influx_sign = (grid_.face_cells[2*f] == upwind_cell) ? -1.0 : 1.0;
                    const double part_influx = influx_sign * darcyflux_[f] / double(num_parts);
                    if (part_influx > 0.0) {
                        influx.push_back(part_influx);
                        node_pos_influx.push_back(node_pos2);
                    }
                }
            }
        }

        // Now we may compute the weighting of the upwind terms.
        const int num_parts = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
        const double outflux_sign = (grid_.face_cells[2*face] == upwind_cell) ? 1.0 : -1.0;
        const double part_outflux = outflux_sign * darcyflux_[face] / double(num_parts);
        const double sum_influx = std::accumulate(influx.begin(), influx.end(), 0.0);
        const double w_factor = weightFunc(sum_influx / part_outflux);
        const int num_influx = influx.size();
        std::vector<double> w(num_influx);
        face_term = 0.0;
        for (int ii = 0; ii < num_influx; ++ii) {
            w[ii] = (influx[ii] / sum_influx) * w_factor;
            face_term += w[ii] * face_part_tof_[node_pos_influx[ii]];
        }
        const double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
        cell_term_factor = 1.0 - sum_w;
        const double tol = 1e-5;
        if (cell_term_factor < -tol && cell_term_factor > 1.0 + tol) {
            OPM_THROW(std::logic_error, "cell_term_factor outside [0,1]: " << cell_term_factor);
        }
        cell_term_factor = std::min(std::max(cell_term_factor, 0.0), 1.0);
        assert(cell_term_factor >= 0.0);
        assert(cell_term_factor <= 1.0);
    }

} // namespace Opm
