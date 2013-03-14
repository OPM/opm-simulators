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

#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/utility/RootFinders.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <fstream>
#include <iterator>
#include <numeric>


#define EXPERIMENT_GAUSS_SEIDEL


namespace Opm
{

    // Choose error policy for scalar solves here.
    typedef RegulaFalsi<WarnAndContinueOnError> RootFinder;


    TransportModelTwophase::TransportModelTwophase(const UnstructuredGrid& grid,
                                                   const Opm::IncompPropertiesInterface& props,
                                                   const double tol,
                                                   const int maxit)
        : grid_(grid),
          props_(props),
          tol_(tol),
          maxit_(maxit),
          darcyflux_(0),
          source_(0),
          dt_(0.0),
          saturation_(grid.number_of_cells, -1.0),
          fractionalflow_(grid.number_of_cells, -1.0),
          reorder_iterations_(grid.number_of_cells, 0),
          mob_(2*grid.number_of_cells, -1.0)
#ifdef EXPERIMENT_GAUSS_SEIDEL
        , ia_upw_(grid.number_of_cells + 1, -1),
          ja_upw_(grid.number_of_faces, -1),
          ia_downw_(grid.number_of_cells + 1, -1),
          ja_downw_(grid.number_of_faces, -1)
#endif
    {
        if (props.numPhases() != 2) {
            THROW("Property object must have 2 phases");
        }
        visc_ = props.viscosity();
        int num_cells = props.numCells();
        smin_.resize(props.numPhases()*num_cells);
        smax_.resize(props.numPhases()*num_cells);
        std::vector<int> cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            cells[i] = i;
        }
        props.satRange(props.numCells(), &cells[0], &smin_[0], &smax_[0]);
    }

    void TransportModelTwophase::solve(const double* darcyflux,
                                       const double* porevolume,
                                       const double* source,
                                       const double dt,
                                       std::vector<double>& saturation)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        source_ = source;
        dt_ = dt;
        toWaterSat(saturation, saturation_);

#ifdef EXPERIMENT_GAUSS_SEIDEL
        std::vector<int> seq(grid_.number_of_cells);
        std::vector<int> comp(grid_.number_of_cells + 1);
        int ncomp;
        compute_sequence_graph(&grid_, darcyflux_,
                               &seq[0], &comp[0], &ncomp,
                               &ia_upw_[0], &ja_upw_[0]);
        const int nf = grid_.number_of_faces;
        std::vector<double> neg_darcyflux(nf);
        std::transform(darcyflux, darcyflux + nf, neg_darcyflux.begin(), std::negate<double>());
        compute_sequence_graph(&grid_, &neg_darcyflux[0],
                               &seq[0], &comp[0], &ncomp,
                               &ia_downw_[0], &ja_downw_[0]);
#endif
        std::fill(reorder_iterations_.begin(),reorder_iterations_.end(),0);
        reorderAndTransport(grid_, darcyflux);
        toBothSat(saturation_, saturation);
    }


    const std::vector<int>& TransportModelTwophase::getReorderIterations() const
    {
        return reorder_iterations_;
    }


    // Residual function r(s) for a single-cell implicit Euler transport
    //
    //     r(s) = s - s0 + dt/pv*( influx + outflux*f(s) )
    //
    // where influx is water influx, outflux is total outflux.
    // Influxes are negative, outfluxes positive.
    struct TransportModelTwophase::Residual
    {
        int cell;
        double s0;
        double influx;    // sum_j min(v_ij, 0)*f(s_j) + q_w
        double outflux;   // sum_j max(v_ij, 0) - q
        double comp_term; // q - sum_j v_ij
        double dtpv;    // dt/pv(i)
        const TransportModelTwophase& tm;
        explicit Residual(const TransportModelTwophase& tmodel, int cell_index)
            : tm(tmodel)
        {
            cell    = cell_index;
            s0      = tm.saturation_[cell];
            double src_flux       = -tm.source_[cell];
            bool src_is_inflow = src_flux < 0.0;
            influx  =  src_is_inflow ? src_flux : 0.0;
            outflux = !src_is_inflow ? src_flux : 0.0;
            comp_term = tm.source_[cell];   // Note: this assumes that all source flux is water.
            dtpv    = tm.dt_/tm.porevolume_[cell];

            for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
                int f = tm.grid_.cell_faces[i];
                double flux;
                int other;
                // Compute cell flux
                if (cell == tm.grid_.face_cells[2*f]) {
                    flux  = tm.darcyflux_[f];
                    other = tm.grid_.face_cells[2*f+1];
                } else {
                    flux  =-tm.darcyflux_[f];
                    other = tm.grid_.face_cells[2*f];
                }
                // Add flux to influx or outflux, if interior.
                if (other != -1) {
                    if (flux < 0.0) {
                        influx  += flux*tm.fractionalflow_[other];
                    } else {
                        outflux += flux;
                    }
                    comp_term -= flux;
                }
            }
        }
        double operator()(double s) const
        {
            return s - s0 + dtpv*(outflux*tm.fracFlow(s, cell) + influx + s*comp_term);
        }
    };


    void TransportModelTwophase::solveSingleCell(const int cell)
    {
        Residual res(*this, cell);
        // const double r0 = res(saturation_[cell]);
        // if (std::fabs(r0) < tol_) {
        //     return;
        // }
        int iters_used = 0;
        // saturation_[cell] = modifiedRegulaFalsi(res, smin_[2*cell], smax_[2*cell], maxit_, tol_, iters_used);
        saturation_[cell] = RootFinder::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used);
        // add if it is iteration on an out loop
        reorder_iterations_[cell] = reorder_iterations_[cell] + iters_used;
        fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
    }

    // namespace {
    //  class TofComputer
    //  {
    //  public:
    //      TofComputer(const int num_cells,
    //                  const int* ia,
    //                  const int* ja,
    //                  const int startcell,
    //                  std::vector<int>& tof)
    //          : ia_(ia),
    //            ja_(ja)
    //      {
    //          tof.clear();
    //          tof.resize(num_cells, num_cells);
    //          tof[startcell] = 0;
    //          tof_ = &tof[0];
    //          visitTof(startcell);
    //      }

    //  private:
    //      const int* ia_;
    //      const int* ja_;
    //      int* tof_;

    //      void visitTof(const int cell)
    //      {
    //          for (int j = ia_[cell]; j < ia_[cell+1]; ++j) {
    //              const int nb_cell = ja_[j];
    //              if (tof_[nb_cell] > tof_[cell] + 1) {
    //                  tof_[nb_cell] = tof_[cell] + 1;
    //                  visitTof(nb_cell);
    //              }
    //          }
    //      }

    //  };
    // } // anon namespace


    void TransportModelTwophase::solveMultiCell(const int num_cells, const int* cells)
    {
        // std::ofstream os("dump");
        // std::copy(cells, cells + num_cells, std::ostream_iterator<double>(os, "\n"));

        // Experiment: try a breath-first search to build a more suitable ordering.
        // Verdict: failed to improve #iterations.
        // {
        //     std::vector<int> pos(grid_.number_of_cells, -1);
        //     for (int i = 0; i < num_cells; ++i) {
        //      const int cell = cells[i];
        //      pos[cell] = i;
        //     }
        //     std::vector<int> done_pos(num_cells, 0);
        //     std::vector<int> upstream_pos;
        //     std::vector<int> new_pos;
        //     upstream_pos.push_back(0);
        //     done_pos[0] = 1;
        //     int current = 0;
        //     while (int(new_pos.size()) < num_cells) {
        //      const int i = upstream_pos[current++];
        //      new_pos.push_back(i);
        //      const int cell = cells[i];
        //      for (int j = ia_[cell]; j < ia_[cell+1]; ++j) {
        //          const int opos = pos[ja_[j]];
        //          if (!done_pos[opos]) {
        //              upstream_pos.push_back(opos);
        //              done_pos[opos] = 1;
        //          }
        //      }
        //     }
        //     std::reverse(new_pos.begin(), new_pos.end());
        //     std::copy(new_pos.begin(), new_pos.end(), const_cast<int*>(cells));
        // }

        // Experiment: try a random ordering.
        // Verdict: amazingly, reduced #iterations by approx. 25%!
        // int* c = const_cast<int*>(cells);
        // std::random_shuffle(c, c + num_cells);

        // Experiment: compute topological tof from first cell.
        // Verdict: maybe useful, not tried to exploit it yet.
        // std::vector<int> tof;
        // TofComputer comp(grid_.number_of_cells, &ia_[0], &ja_[0], cells[0], tof);
        // std::ofstream tofdump("tofdump");
        // std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tofdump, "\n"));

        // Experiment: implement a metric measuring badness of ordering
        //             as average distance in (cyclic) ordering from
        //             upstream neighbours.
        // Verdict: does not seem to predict #iterations very well, if at all.
        // std::vector<int> pos(grid_.number_of_cells, -1);
        // for (int i = 0; i < num_cells; ++i) {
        //     const int cell = cells[i];
        //     pos[cell] = i;
        // }
        // double diffsum = 0;
        // for (int i = 0; i < num_cells; ++i) {
        //     const int cell = cells[i];
        //     int num_upstream = 0;
        //     int loc_diffsum = 0;
        //     for (int j = ia_[cell]; j < ia_[cell+1]; ++j) {
        //      const int opos = pos[ja_[j]];
        //      if (opos == -1) {
        //          std::cout << "Hmmm" << std::endl;
        //          continue;
        //      }
        //      ++num_upstream;
        //      const int diff = (i - opos + num_cells) % num_cells;
        //      loc_diffsum += diff;
        //     }
        //     diffsum += double(loc_diffsum)/double(num_upstream);
        // }
        // std::cout << "Average distance from upstream neighbours: " << diffsum/double(num_cells)
        //        << std::endl;

#ifdef EXPERIMENT_GAUSS_SEIDEL
        // Experiment: when a cell changes more than the tolerance,
        //             mark all downwind cells as needing updates. After
        //             computing a single update in each cell, use marks
        //             to guide further updating. Clear mark in cell when
        //             its solution gets updated.
        // Verdict: this is a good one! Approx. halved total time.
        std::vector<int> needs_update(num_cells, 1);
        // This one also needs the mapping from all cells to
        // the strongly connected subset to filter out connections
        std::vector<int> pos(grid_.number_of_cells, -1);
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            pos[cell] = i;
        }

        // Note: partially copied from below.
        const double tol = 1e-9;
        const int max_iters = 300;
        // Must store s0 before we start.
        std::vector<double> s0(num_cells);
        // Must set initial fractional flows before we start.
        // Also, we compute the # of upstream neighbours.
        // std::vector<int> num_upstream(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
            s0[i] = saturation_[cell];
            // num_upstream[i] = ia_upw_[cell + 1] - ia_upw_[cell];
        }
        // Solve once in each cell.
        // std::vector<int> fully_marked_stack;
        // fully_marked_stack.reserve(num_cells);
        int num_iters = 0;
        int update_count = 0; // Change name/meaning to cells_updated?
        do {
            update_count = 0; // Must reset count for every iteration.
            for (int i = 0; i < num_cells; ++i) {
                // while (!fully_marked_stack.empty()) {
                //     // std::cout << "# fully marked cells = " << fully_marked_stack.size() << std::endl;
                //     const int fully_marked_ci = fully_marked_stack.back();
                //     fully_marked_stack.pop_back();
                //     ++update_count;
                //     const int cell = cells[fully_marked_ci];
                //     const double old_s = saturation_[cell];
                //     saturation_[cell] = s0[fully_marked_ci];
                //     solveSingleCell(cell);
                //     const double s_change = std::fabs(saturation_[cell] - old_s);
                //     if (s_change > tol) {
                //      // Mark downwind cells.
                //      for (int j = ia_downw_[cell]; j < ia_downw_[cell+1]; ++j) {
                //          const int downwind_cell = ja_downw_[j];
                //          int ci = pos[downwind_cell];
                //          ++needs_update[ci];
                //          if (needs_update[ci] == num_upstream[ci]) {
                //              fully_marked_stack.push_back(ci);
                //          }
                //      }
                //     }
                //     // Unmark this cell.
                //     needs_update[fully_marked_ci] = 0;
                // }
                if (!needs_update[i]) {
                    continue;
                }
                ++update_count;
                const int cell = cells[i];
                const double old_s = saturation_[cell];
                saturation_[cell] = s0[i];
                solveSingleCell(cell);
                const double s_change = std::fabs(saturation_[cell] - old_s);
                if (s_change > tol) {
                    // Mark downwind cells.
                    for (int j = ia_downw_[cell]; j < ia_downw_[cell+1]; ++j) {
                        const int downwind_cell = ja_downw_[j];
                        int ci = pos[downwind_cell];
                        if (ci != -1) {
                            needs_update[ci] = 1;
                        }
                        // ++needs_update[ci];
                        // if (needs_update[ci] == num_upstream[ci]) {
                        //     fully_marked_stack.push_back(ci);
                        // }
                    }
                }
                // Unmark this cell.
                needs_update[i] = 0;
            }
            // std::cout << "Iter = " << num_iters << "    update_count = " << update_count
            //        << "    # marked cells = "
            //        << std::accumulate(needs_update.begin(), needs_update.end(), 0) << std::endl;
        } while (update_count > 0 && ++num_iters < max_iters);

        // Done with iterations, check if we succeeded.
        if (update_count > 0) {
            THROW("In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Remaining update count = " << update_count);
        }
        std::cout << "Solved " << num_cells << " cell multicell problem in "
                  << num_iters << " iterations." << std::endl;

#else
        double max_s_change = 0.0;
        const double tol = 1e-9;
        const int max_iters = 300;
        int num_iters = 0;
        // Must store s0 before we start.
        std::vector<double> s0(num_cells);
        // Must set initial fractional flows before we start.
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
            s0[i] = saturation_[cell];
        }
        do {
            max_s_change = 0.0;
            for (int i = 0; i < num_cells; ++i) {
                const int cell = cells[i];
                const double old_s = saturation_[cell];
                saturation_[cell] = s0[i];
                solveSingleCell(cell);
                double s_change = std::fabs(saturation_[cell] - old_s);
                // std::cout << "cell = " << cell << "    delta s = " << s_change << std::endl;
                if (max_s_change < s_change) {
                    max_s_change = s_change;
                }
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change
            //        << "    in cell " << max_change_cell << std::endl;
        } while (max_s_change > tol && ++num_iters < max_iters);
        if (max_s_change > tol) {
            THROW("In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_s_change);
        }
        std::cout << "Solved " << num_cells << " cell multicell problem in "
                  << num_iters << " iterations." << std::endl;
#endif // EXPERIMENT_GAUSS_SEIDEL
    }

    double TransportModelTwophase::fracFlow(double s, int cell) const
    {
        double sat[2] = { s, 1.0 - s };
        double mob[2];
        props_.relperm(1, sat, &cell, mob, 0);
        mob[0] /= visc_[0];
        mob[1] /= visc_[1];
        return mob[0]/(mob[0] + mob[1]);
    }





    // Residual function r(s) for a single-cell implicit Euler gravity segregation
    //
    //     r(s) = s - s0 + dt/pv*sum_{j adj i}( gravmod_ij * gf_ij ).
    //
    struct TransportModelTwophase::GravityResidual
    {
        int cell;
        int nbcell[2];
        double s0;
        double dtpv;    // dt/pv(i)
        double gf[2];
        const TransportModelTwophase& tm;
        explicit GravityResidual(const TransportModelTwophase& tmodel,
                                 const std::vector<int>& cells,
                                 const int pos,
                                 const double* gravflux) // Always oriented towards next in column. Size = colsize - 1.
            : tm(tmodel)
        {
            cell = cells[pos];
            nbcell[0] = -1;
            gf[0] = 0.0;
            if (pos > 0) {
                nbcell[0] = cells[pos - 1];
                gf[0] = -gravflux[pos - 1];
            }
            nbcell[1] = -1;
            gf[1] = 0.0;
            if (pos < int(cells.size() - 1)) {
                nbcell[1] = cells[pos + 1];
                gf[1] = gravflux[pos];
            }
            s0      = tm.saturation_[cell];
            dtpv    = tm.dt_/tm.porevolume_[cell];

        }
        double operator()(double s) const
        {
            double res = s - s0;
            double mobcell[2];
            tm.mobility(s, cell, mobcell);
            for (int nb = 0; nb < 2; ++nb) {
                if (nbcell[nb] != -1) {
                    double m[2];
                    if (gf[nb] < 0.0) {
                        m[0] = mobcell[0];
                        m[1] = tm.mob_[2*nbcell[nb] + 1];
                    } else {
                        m[0] = tm.mob_[2*nbcell[nb]];
                        m[1] = mobcell[1];
                    }
                    if (m[0] + m[1] > 0.0) {
                        res += -dtpv*gf[nb]*m[0]*m[1]/(m[0] + m[1]);
                    }
                }
            }
            return res;
        }
    };

    void TransportModelTwophase::mobility(double s, int cell, double* mob) const
    {
        double sat[2] = { s, 1.0 - s };
        props_.relperm(1, sat, &cell, mob, 0);
        mob[0] /= visc_[0];
        mob[1] /= visc_[1];
    }



    void TransportModelTwophase::initGravity(const double* grav)
    {
        // Set up gravflux_ = T_ij g (rho_w - rho_o) (z_i - z_j)
        std::vector<double> htrans(grid_.cell_facepos[grid_.number_of_cells]);
        const int nf = grid_.number_of_faces;
        const int dim = grid_.dimensions;
        gravflux_.resize(nf);
        tpfa_htrans_compute(const_cast<UnstructuredGrid*>(&grid_), props_.permeability(), &htrans[0]);
        tpfa_trans_compute(const_cast<UnstructuredGrid*>(&grid_), &htrans[0], &gravflux_[0]);
        const double delta_rho = props_.density()[0] - props_.density()[1];
        for (int f = 0; f < nf; ++f) {
            const int* c = &grid_.face_cells[2*f];
            double gdz = 0.0;
            if (c[0] != -1 && c[1] != -1) {
                for (int d = 0; d < dim; ++d) {
                    gdz += grav[d]*(grid_.cell_centroids[dim*c[0] + d] - grid_.cell_centroids[dim*c[1] + d]);
                }
            }
            gravflux_[f] *= delta_rho*gdz;
        }
    }



    void TransportModelTwophase::solveSingleCellGravity(const std::vector<int>& cells,
                                                        const int pos,
                                                        const double* gravflux)
    {
        const int cell = cells[pos];
        GravityResidual res(*this, cells, pos, gravflux);
        if (std::fabs(res(saturation_[cell])) > tol_) {
            int iters_used = 0;
            saturation_[cell] = RootFinder::solve(res, smin_[2*cell], smax_[2*cell], maxit_, tol_, iters_used);
            reorder_iterations_[cell] = reorder_iterations_[cell] + iters_used;
        }
        saturation_[cell] = std::min(std::max(saturation_[cell], smin_[2*cell]), smax_[2*cell]);
        mobility(saturation_[cell], cell, &mob_[2*cell]);
    }



    int TransportModelTwophase::solveGravityColumn(const std::vector<int>& cells)
    {
        // Set up column gravflux.
        const int nc = cells.size();
        std::vector<double> col_gravflux(nc - 1);
        for (int ci = 0; ci < nc - 1; ++ci) {
            const int cell = cells[ci];
            const int next_cell = cells[ci + 1];
            for (int j = grid_.cell_facepos[cell]; j < grid_.cell_facepos[cell+1]; ++j) {
                const int face = grid_.cell_faces[j];
                const int c1 = grid_.face_cells[2*face + 0];
                const int c2 = grid_.face_cells[2*face + 1];
                if (c1 == next_cell || c2 == next_cell) {
                    const double gf = gravflux_[face];
                    col_gravflux[ci] = (c1 == cell) ? gf : -gf;
                }
            }
        }

        // Store initial saturation s0
        s0_.resize(nc);
        for (int ci = 0; ci < nc; ++ci) {
            s0_[ci] = saturation_[cells[ci]];
        }

        // Solve single cell problems, repeating if necessary.
        double max_s_change = 0.0;
        int num_iters = 0;
        do {
            max_s_change = 0.0;
            for (int ci = 0; ci < nc; ++ci) {
                const int ci2 = nc - ci - 1;
                double old_s[2] = { saturation_[cells[ci]],
                                    saturation_[cells[ci2]] };
                saturation_[cells[ci]] = s0_[ci];
                solveSingleCellGravity(cells, ci, &col_gravflux[0]);
                saturation_[cells[ci2]] = s0_[ci2];
                solveSingleCellGravity(cells, ci2, &col_gravflux[0]);
                max_s_change = std::max(max_s_change, std::max(std::fabs(saturation_[cells[ci]] - old_s[0]),
                                                               std::fabs(saturation_[cells[ci2]] - old_s[1])));
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change << std::endl;
        } while (max_s_change > tol_ && ++num_iters < maxit_);

        if (max_s_change > tol_) {
            THROW("In solveGravityColumn(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_s_change);
        }
        return num_iters + 1;
    }



    void TransportModelTwophase::solveGravity(const std::vector<std::vector<int> >& columns,
                                              const double* porevolume,
                                              const double dt,
                                              std::vector<double>& saturation)
    {
        // Initialize mobilities.
        const int nc = grid_.number_of_cells;
        std::vector<int> cells(nc);
        for (int c = 0; c < nc; ++c) {
            cells[c] = c;
        }
        mob_.resize(2*nc);
        props_.relperm(cells.size(), &saturation[0], &cells[0], &mob_[0], 0);
        const double* mu = props_.viscosity();
        for (int c = 0; c < nc; ++c) {
            mob_[2*c] /= mu[0];
            mob_[2*c + 1] /= mu[1];
        }

        // Set up other variables.
        porevolume_ = porevolume;
        dt_ = dt;
        toWaterSat(saturation, saturation_);

        // Solve on all columns.
        int num_iters = 0;
        for (std::vector<std::vector<int> >::size_type i = 0; i < columns.size(); i++) {
            // std::cout << "==== new column" << std::endl;
            num_iters += solveGravityColumn(columns[i]);
        }
        std::cout << "Gauss-Seidel column solver average iterations: "
                  << double(num_iters)/double(columns.size()) << std::endl;

        toBothSat(saturation_, saturation);
    }

} // namespace Opm



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
