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
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iterator>

namespace Opm
{


    /// @brief Computes pore volume of all cells in a grid.
    /// @param[in]  grid      a grid
    /// @param[in]  porosity  array of grid.number_of_cells porosity values
    /// @param[out] porevol   the pore volume by cell.
    void computePorevolume(const UnstructuredGrid& grid,
                           const double* porosity,
                           std::vector<double>& porevol)
    {
        int num_cells = grid.number_of_cells;
        porevol.resize(num_cells);
        std::transform(porosity, porosity + num_cells,
                       grid.cell_volumes,
                       porevol.begin(),
                       std::multiplies<double>());
    }


    /// @brief Computes pore volume of all cells in a grid, with rock compressibility effects.
    /// @param[in]  grid      a grid
    /// @param[in]  porosity  array of grid.number_of_cells porosity values
    /// @param[in]  rock_comp rock compressibility properties
    /// @param[in]  pressure  pressure by cell
    /// @param[out] porevol   the pore volume by cell.
    void computePorevolume(const UnstructuredGrid& grid,
                           const double* porosity,
                           const RockCompressibility& rock_comp,
                           const std::vector<double>& pressure,
                           std::vector<double>& porevol)
    {
        int num_cells = grid.number_of_cells;
        porevol.resize(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            porevol[i] = porosity[i]*grid.cell_volumes[i]*rock_comp.poroMult(pressure[i]);
        }
    }

    /// @brief Computes porosity of all cells in a grid, with rock compressibility effects.
    /// @param[in]  grid               a grid
    /// @param[in]  porosity_standard  array of grid.number_of_cells porosity values (at standard conditions)
    /// @param[in]  rock_comp          rock compressibility properties
    /// @param[in]  pressure           pressure by cell
    /// @param[out] porosity           porosity (at reservoir condition)
    void computePorosity(const UnstructuredGrid& grid,
                         const double* porosity_standard,
                         const RockCompressibility& rock_comp,
                         const std::vector<double>& pressure,
                         std::vector<double>& porosity)
    {
        int num_cells = grid.number_of_cells;
        porosity.resize(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            porosity[i] = porosity_standard[i]*rock_comp.poroMult(pressure[i]);
        }
    }


    /// @brief Computes total saturated volumes over all grid cells.
    /// @param[in]  pv        the pore volume by cell.
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[out] sat_vol   must point to a valid array with P elements,
    ///                       where P = s.size()/pv.size().
    ///                       For each phase p, we compute
    ///                       sat_vol_p = sum_i s_p_i pv_i
    void computeSaturatedVol(const std::vector<double>& pv,
                             const std::vector<double>& s,
                             double* sat_vol)
    {
        const int num_cells = pv.size();
        const int np = s.size()/pv.size();
        if (int(s.size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of s and pv vectors do not match.");
        }
        std::fill(sat_vol, sat_vol + np, 0.0);
        for (int c = 0; c < num_cells; ++c) {
            for (int p = 0; p < np; ++p) {
                sat_vol[p] += pv[c]*s[np*c + p];
            }
        }
    }


    /// @brief Computes average saturations over all grid cells.
    /// @param[in]  pv        the pore volume by cell.
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[out] aver_sat  must point to a valid array with P elements,
    ///                       where P = s.size()/pv.size().
    ///                       For each phase p, we compute
    ///                       aver_sat_p = (sum_i s_p_i pv_i) / (sum_i pv_i).
    void computeAverageSat(const std::vector<double>& pv,
                           const std::vector<double>& s,
                           double* aver_sat)
    {
        const int num_cells = pv.size();
        const int np = s.size()/pv.size();
        if (int(s.size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of s and pv vectors do not match.");
        }
        double tot_pv = 0.0;
        // Note that we abuse the output array to accumulate the
        // saturated pore volumes.
        std::fill(aver_sat, aver_sat + np, 0.0);
        for (int c = 0; c < num_cells; ++c) {
            tot_pv += pv[c];
            for (int p = 0; p < np; ++p) {
                aver_sat[p] += pv[c]*s[np*c + p];
            }
        }
        // Must divide by pore volumes to get saturations.
        for (int p = 0; p < np; ++p) {
            aver_sat[p] /= tot_pv;
        }
    }


    /// @brief Computes injected and produced volumes of all phases.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// @param[in]  props     fluid and rock properties.
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[in]  src       if < 0: total outflow, if > 0: first phase inflow.
    /// @param[in]  dt        timestep used
    /// @param[out] injected  must point to a valid array with P elements,
    ///                       where P = s.size()/src.size().
    /// @param[out] produced  must also point to a valid array with P elements.
    void computeInjectedProduced(const IncompPropertiesInterface& props,
                                 const std::vector<double>& s,
                                 const std::vector<double>& src,
                                 const double dt,
                                 double* injected,
                                 double* produced)
    {
        const int num_cells = src.size();
        const int np = s.size()/src.size();
        if (int(s.size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of s and src vectors do not match.");
        }
        std::fill(injected, injected + np, 0.0);
        std::fill(produced, produced + np, 0.0);
        const double* visc = props.viscosity();
        std::vector<double> mob(np);
        for (int c = 0; c < num_cells; ++c) {
            if (src[c] > 0.0) {
                injected[0] += src[c]*dt;
            } else if (src[c] < 0.0) {
                const double flux = -src[c]*dt;
                const double* sat = &s[np*c];
                props.relperm(1, sat, &c, &mob[0], 0);
                double totmob = 0.0;
                for (int p = 0; p < np; ++p) {
                    mob[p] /= visc[p];
                    totmob += mob[p];
                }
                for (int p = 0; p < np; ++p) {
                    produced[p] += (mob[p]/totmob)*flux;
                }
            }
        }
    }



    /// @brief Computes total mobility for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobilities.
    void computeTotalMobility(const Opm::IncompPropertiesInterface& props,
                              const std::vector<int>& cells,
                              const std::vector<double>& s,
                              std::vector<double>& totmob)
    {
        std::vector<double> pmobc;

        computePhaseMobilities(props, cells, s, pmobc);

        const std::size_t                 np = props.numPhases();
        const std::vector<int>::size_type nc = cells.size();

        std::vector<double>(cells.size(), 0.0).swap(totmob);

        for (std::vector<int>::size_type c = 0; c < nc; ++c) {
            for (std::size_t p = 0; p < np; ++p) {
                totmob[ c ] += pmobc[c*np + p];
            }
        }
    }


    /// @brief Computes total mobility and omega for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobility
    /// @param[out] omega     fractional-flow weighted fluid densities.
    void computeTotalMobilityOmega(const Opm::IncompPropertiesInterface& props,
                                   const std::vector<int>& cells,
                                   const std::vector<double>& s,
                                   std::vector<double>& totmob,
                                   std::vector<double>& omega)
    {
        std::vector<double> pmobc;

        computePhaseMobilities(props, cells, s, pmobc);

        const std::size_t                 np = props.numPhases();
        const std::vector<int>::size_type nc = cells.size();

        std::vector<double>(cells.size(), 0.0).swap(totmob);
        std::vector<double>(cells.size(), 0.0).swap(omega );

        const double* rho = props.density();
        for (std::vector<int>::size_type c = 0; c < nc; ++c) {
            for (std::size_t p = 0; p < np; ++p) {
                totmob[ c ] += pmobc[c*np + p];
                omega [ c ] += pmobc[c*np + p] * rho[ p ];
            }

            omega[ c ] /= totmob[ c ];
        }
    }


    /// @brief Computes phase mobilities for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] pmobc     phase mobilities (for all phases).
    void computePhaseMobilities(const Opm::IncompPropertiesInterface& props,
                                const std::vector<int>&               cells,
                                const std::vector<double>&            s    ,
                                std::vector<double>&                  pmobc)
    {
        const std::vector<int>::size_type nc = cells.size();
        const std::size_t                 np = props.numPhases();

        assert(s.size() == nc * np);

        std::vector<double>(nc * np, 0.0).swap(pmobc );
        double*                                dpmobc = 0;
        props.relperm(static_cast<const int>(nc), &s[0], &cells[0],
                      &pmobc[0], dpmobc);

        const double*                 mu  = props.viscosity();
        std::vector<double>::iterator lam = pmobc.begin();
        for (std::vector<int>::size_type c = 0; c < nc; ++c) {
            for (std::size_t p = 0; p < np; ++p, ++lam) {
                *lam /= mu[ p ];
            }
        }
    }

    /// Computes the fractional flow for each cell in the cells argument
    /// @param[in] props                rock and fluid properties
    /// @param[in] cells                cells with which the saturation values are associated
    /// @param[in] saturations          saturation values (for all phases)
    /// @param[out] fractional_flow     the fractional flow for each phase for each cell.

    void computeFractionalFlow(const Opm::IncompPropertiesInterface& props,
                               const std::vector<int>& cells,
                               const std::vector<double>& saturations,
                               std::vector<double>& fractional_flows)
    {
        const int num_phases = props.numPhases();

        computePhaseMobilities(props, cells, saturations, fractional_flows);

        for (std::vector<int>::size_type i = 0; i < cells.size(); ++i) {
            double phase_sum = 0.0;
            for (int phase = 0; phase < num_phases; ++phase) {
                phase_sum += fractional_flows[i * num_phases + phase];
            }
            for (int phase = 0; phase < num_phases; ++phase) {
                fractional_flows[i * num_phases + phase] /= phase_sum;
            }
        }
    }

    /// Compute two-phase transport source terms from face fluxes,
    /// and pressure equation source terms. This puts boundary flows
    /// into the source terms for the transport equation.
    /// \param[in]  grid          The grid used.
    /// \param[in]  src           Pressure eq. source terms. The sign convention is:
    ///                           (+) positive  total inflow (positive velocity divergence)
    ///                           (-) negative  total outflow
    /// \param[in]  faceflux      Signed face fluxes, typically the result from a flow solver.
    /// \param[in]  inflow_frac   Fraction of inflow (boundary and source terms) that consists of first phase.
    ///                           Example: if only water is injected, inflow_frac == 1.0.
    ///                           Note: it is not possible (with this method) to use different fractions
    ///                           for different inflow sources, be they source terms of boundary flows.
    /// \param[in]  wells         Wells data structure.
    /// \param[in]  well_perfrates  Volumetric flow rates per well perforation.
    /// \param[out] transport_src The transport source terms. They are to be interpreted depending on sign:
    ///                           (+) positive  inflow of first phase (water)
    ///                           (-) negative  total outflow of both phases
    void computeTransportSource(const UnstructuredGrid& grid,
                                const std::vector<double>& src,
                                const std::vector<double>& faceflux,
                                const double inflow_frac,
                                const Wells* wells,
                                const std::vector<double>& well_perfrates,
                                std::vector<double>& transport_src)
    {
        int nc = grid.number_of_cells;
        transport_src.resize(nc);
        // Source term and boundary contributions.
        for (int c = 0; c < nc; ++c) {
            transport_src[c] = 0.0;
            transport_src[c] += src[c] > 0.0 ? inflow_frac*src[c] : src[c];
            for (int hf = grid.cell_facepos[c]; hf < grid.cell_facepos[c + 1]; ++hf) {
                int f = grid.cell_faces[hf];
                const int* f2c = &grid.face_cells[2*f];
                double bdy_influx = 0.0;
                if (f2c[0] == c && f2c[1] == -1) {
                    bdy_influx = -faceflux[f];
                } else if (f2c[0] == -1 && f2c[1] == c) {
                    bdy_influx = faceflux[f];
                }
                if (bdy_influx != 0.0) {
                    transport_src[c] += bdy_influx > 0.0 ? inflow_frac*bdy_influx : bdy_influx;
                }
            }
        }

        // Well contributions.
        if (wells) {
            const int nw = wells->number_of_wells;
            const int np = wells->number_of_phases;
            if (np != 2) {
                OPM_THROW(std::runtime_error, "computeTransportSource() requires a 2 phase case.");
            }
            for (int w = 0; w < nw; ++w) {
                const double* comp_frac = wells->comp_frac + np*w;
                for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                    const int perf_cell = wells->well_cells[perf];
                    double perf_rate = well_perfrates[perf];
                    if (perf_rate > 0.0) {
                        // perf_rate is a total inflow rate, we want a water rate.
                        if (wells->type[w] != INJECTOR) {
                            std::cout << "**** Warning: crossflow in well "
                                      << w << " perf " << perf - wells->well_connpos[w]
                                      << " ignored. Rate was "
                                      << perf_rate/Opm::unit::day << " m^3/day." << std::endl;
                            perf_rate = 0.0;
                        } else {
                            assert(std::fabs(comp_frac[0] + comp_frac[1] - 1.0) < 1e-6);
                            perf_rate *= comp_frac[0];
                        }
                    }
                    transport_src[perf_cell] += perf_rate;
                }
            }
        }
    }

    /// @brief Estimates a scalar cell velocity from face fluxes.
    /// @param[in]  grid            a grid
    /// @param[in]  face_flux       signed per-face fluxes
    /// @param[out] cell_velocity   the estimated velocities.
    void estimateCellVelocity(const UnstructuredGrid& grid,
                              const std::vector<double>& face_flux,
                              std::vector<double>& cell_velocity)
    {
        estimateCellVelocity(grid.number_of_cells,
                             grid.number_of_faces,
                             grid.face_centroids,
                             UgGridHelpers::faceCells(grid),
                             grid.cell_centroids,
                             grid.cell_volumes,
                             grid.dimensions,
                             face_flux,
                             cell_velocity);
        
    }

    /// Extract a vector of water saturations from a vector of
    /// interleaved water and oil saturations.
    void toWaterSat(const std::vector<double>& sboth,
                    std::vector<double>& sw)
    {
        int num = sboth.size()/2;
        sw.resize(num);
        for (int i = 0; i < num; ++i) {
            sw[i] = sboth[2*i];
        }
    }

    /// Make a a vector of interleaved water and oil saturations from
    /// a vector of water saturations.
    void toBothSat(const std::vector<double>& sw,
                   std::vector<double>& sboth)
    {
        int num = sw.size();
        sboth.resize(2*num);
        for (int i = 0; i < num; ++i) {
            sboth[2*i] = sw[i];
            sboth[2*i + 1] = 1.0 - sw[i];
        }
    }


    /// Create a src vector equivalent to a wells structure.
    /// For this to be valid, the wells must be all rate-controlled and
    /// single-perforation.
    void wellsToSrc(const Wells& wells, const int num_cells, std::vector<double>& src)
    {
        const int np = wells.number_of_phases;
        if (np != 2) {
            OPM_THROW(std::runtime_error, "wellsToSrc() requires a 2 phase case.");
        }
        src.resize(num_cells);
        for (int w = 0; w < wells.number_of_wells; ++w) {
            const int cur = well_controls_get_current(wells.ctrls[w]);
            if (well_controls_get_num(wells.ctrls[w]) != 1) {
                OPM_MESSAGE("In wellsToSrc(): well has more than one control, all but current control will be ignored.");
            }
            if (well_controls_iget_type(wells.ctrls[w] , cur) != RESERVOIR_RATE) {
                OPM_THROW(std::runtime_error, "In wellsToSrc(): well is something other than RESERVOIR_RATE.");
            }
            if (wells.well_connpos[w+1] - wells.well_connpos[w] != 1) {
                OPM_THROW(std::runtime_error, "In wellsToSrc(): well has multiple perforations.");
            }
            {
                const double * distr = well_controls_iget_distr( wells.ctrls[w] , cur);
                for (int p = 0; p < np; ++p) {
                    if (distr[p] != 1.0) {
                        OPM_THROW(std::runtime_error, "In wellsToSrc(): well not controlled on total rate.");
                    }
                }
            }
            double flow = well_controls_iget_target(wells.ctrls[w] , cur);
            if (wells.type[w] == INJECTOR) {
                flow *= wells.comp_frac[np*w + 0]; // Obtaining water rate for inflow source.
            }
            const int cell = wells.well_cells[wells.well_connpos[w]];
            src[cell] = flow;
        }
    }


    void computeWDP(const Wells& wells, const UnstructuredGrid& grid, const std::vector<double>& saturations,
                    const double* densities, const double gravity, const bool per_grid_cell,
                    std::vector<double>& wdp)
    {
        const int nw = wells.number_of_wells;
        const size_t np = per_grid_cell ?
            saturations.size()/grid.number_of_cells
            : saturations.size()/wells.well_connpos[nw];
        // Simple for now:
        for (int i = 0; i < nw; i++) {
            double depth_ref = wells.depth_ref[i];
            for (int j = wells.well_connpos[i]; j < wells.well_connpos[i + 1]; j++) {
                int cell = wells.well_cells[j];

                // Is this correct wrt. depth_ref?
                double cell_depth = grid.cell_centroids[3 * cell + 2];

                double saturation_sum = 0.0;
                for (size_t p = 0; p < np; p++) {
                    if (!per_grid_cell) {
                        saturation_sum += saturations[j * np + p];
                    } else {
                        saturation_sum += saturations[np * cell + p];
                    }
                }
                if (saturation_sum == 0) {
                    saturation_sum = 1.0;
                }
                double density = 0.0;
                for (size_t p = 0; p < np; p++) {
                    if (!per_grid_cell) {
                        density += saturations[j * np + p] * densities[p] / saturation_sum;
                    } else {
                        // Is this a smart way of doing it?
                        density += saturations[np * cell + p] * densities[p] / saturation_sum;
                    }
                }

                // Is the sign correct?
                wdp.push_back(density * (cell_depth - depth_ref) * gravity);
            }
        }
    }


    void computeFlowRatePerWell(const Wells& wells, const std::vector<double>& flow_rates_per_cell,
            std::vector<double>& flow_rates_per_well)
    {
        int index_in_flow_rates = 0;
        for (int w = 0; w < wells.number_of_wells; w++) {
            int number_of_cells = wells.well_connpos[w + 1] - wells.well_connpos[w];
            double flow_sum = 0.0;
            for (int i = 0; i < number_of_cells; i++) {
                flow_sum += flow_rates_per_cell[index_in_flow_rates++];
            }
            flow_rates_per_well.push_back(flow_sum);
        }
    }


    /// Computes the phase flow rate per well
    /// \param[in] wells                        The wells for which the flow rate should be computed
    /// \param[in] flow_rates_per_well_cell     The total flow rate for each cell (ordered the same
    ///                                         way as the wells struct
    /// \param[in] fractional_flows             the fractional flow for each cell in each well
    /// \param[out] phase_flow_per_well         Will contain the phase flow per well
    void computePhaseFlowRatesPerWell(const Wells& wells,
                                      const std::vector<double>& flow_rates_per_well_cell,
                                      const std::vector<double>& fractional_flows,
                                      std::vector<double>& phase_flow_per_well)
    {
        const int np = wells.number_of_phases;
        const int nw = wells.number_of_wells;
        assert(int(flow_rates_per_well_cell.size()) == wells.well_connpos[nw]);
        phase_flow_per_well.resize(nw * np);
        for (int wix = 0; wix < nw; ++wix) {
            for (int phase = 0; phase < np; ++phase) {
                // Reset vector
                phase_flow_per_well[wix*np + phase] = 0.0;
            }
            for (int i = wells.well_connpos[wix]; i < wells.well_connpos[wix + 1]; ++i) {
                const int cell = wells.well_cells[i];
                for (int phase = 0; phase < np; ++phase) {
                    phase_flow_per_well[wix * np + phase] += flow_rates_per_well_cell[i] * fractional_flows[cell * np + phase];
                }
            }
        }
    }


    void Watercut::push(double time, double fraction, double produced)
    {
        data_.push_back(time);
        data_.push_back(fraction);
        data_.push_back(produced);
    }

    void Watercut::write(std::ostream& os) const
    {
        int sz = data_.size() / 3;
        for (int i = 0; i < sz; ++i) {
            os << data_[3 * i] / Opm::unit::day << "   "
                    << data_[3 * i + 1] << "   "
                    << data_[3 * i + 2] << '\n';
        }
    }


    void WellReport::push(const IncompPropertiesInterface& props,
                          const Wells& wells,
                          const std::vector<double>& saturation,
                          const double time,
                          const std::vector<double>& well_bhp,
                          const std::vector<double>& well_perfrates)
    {
        int nw = well_bhp.size();
        assert(nw == wells.number_of_wells);
        int np = props.numPhases();
        const int max_np = 3;
        if (np > max_np) {
            OPM_THROW(std::runtime_error, "WellReport for now assumes #phases <= " << max_np);
        }
        const double* visc = props.viscosity();
        std::vector<double> data_now;
        data_now.reserve(1 + 3*nw);
        data_now.push_back(time/unit::day);
        for (int w = 0; w < nw; ++w) {
            data_now.push_back(well_bhp[w]/(unit::barsa));
            double well_rate_total = 0.0;
            double well_rate_water = 0.0;
            for (int perf = wells.well_connpos[w]; perf < wells.well_connpos[w + 1]; ++perf) {
                const double perf_rate = unit::convert::to(well_perfrates[perf],
                                                           unit::cubic(unit::meter)/unit::day);
                well_rate_total += perf_rate;
                if (perf_rate > 0.0) {
                    // Injection.
                    well_rate_water += perf_rate*wells.comp_frac[0];
                } else {
                    // Production.
                    const int cell = wells.well_cells[perf];
                    double mob[max_np];
                    props.relperm(1, &saturation[2*cell], &cell, mob, 0);
                    double tmob = 0;
                    for(int i = 0; i < np; ++i) {
                        mob[i] /= visc[i];
                        tmob += mob[i];
                    }
                    const double fracflow = mob[0]/tmob;
                    well_rate_water += perf_rate*fracflow;
                }
            }
            data_now.push_back(well_rate_total);
            if (well_rate_total == 0.0) {
                data_now.push_back(0.0);
            } else {
                data_now.push_back(well_rate_water/well_rate_total);
            }
        }
        data_.push_back(data_now);
    }




    void WellReport::push(const BlackoilPropertiesInterface& props,
                          const Wells& wells,
                          const std::vector<double>& p,
                          const std::vector<double>& z,
                          const std::vector<double>& s,
                          const double time,
                          const std::vector<double>& well_bhp,
                          const std::vector<double>& well_perfrates)
    {
        // TODO: refactor, since this is almost identical to the other push().
        int nw = well_bhp.size();
        assert(nw == wells.number_of_wells);
        int np = props.numPhases();
        const int max_np = 3;
        if (np > max_np) {
            OPM_THROW(std::runtime_error, "WellReport for now assumes #phases <= " << max_np);
        }
        std::vector<double> data_now;
        data_now.reserve(1 + 3*nw);
        data_now.push_back(time/unit::day);
        for (int w = 0; w < nw; ++w) {
            data_now.push_back(well_bhp[w]/(unit::barsa));
            double well_rate_total = 0.0;
            double well_rate_water = 0.0;
            for (int perf = wells.well_connpos[w]; perf < wells.well_connpos[w + 1]; ++perf) {
                const double perf_rate = unit::convert::to(well_perfrates[perf],
                                                           unit::cubic(unit::meter)/unit::day);
                well_rate_total += perf_rate;
                if (perf_rate > 0.0) {
                    // Injection.
                    well_rate_water += perf_rate*wells.comp_frac[0];
                } else {
                    // Production.
                    const int cell = wells.well_cells[perf];
                    double mob[max_np];
                    props.relperm(1, &s[np*cell], &cell, mob, 0);
                    double visc[max_np];
                    props.viscosity(1, &p[cell], &z[np*cell], &cell, visc, 0);
                    double tmob = 0;
                    for(int i = 0; i < np; ++i) {
                        mob[i] /= visc[i];
                        tmob += mob[i];
                    }
                    const double fracflow = mob[0]/(tmob);
                    well_rate_water += perf_rate*fracflow;
                }
            }
            data_now.push_back(well_rate_total);
            if (well_rate_total == 0.0) {
                data_now.push_back(0.0);
            } else {
                data_now.push_back(well_rate_water/well_rate_total);
            }
        }
        data_.push_back(data_now);
    }




    void WellReport::write(std::ostream& os) const
    {
        const int sz = data_.size();
        for (int i = 0; i < sz; ++i) {
            std::copy(data_[i].begin(), data_[i].end(), std::ostream_iterator<double>(os, "\t"));
            os << '\n';
        }
    }



} // namespace Opm
