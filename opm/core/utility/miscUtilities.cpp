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

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid.h>
#include <opm/core/newwells.h>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <functional>
#include <cmath>

namespace Opm
{


    /// @brief Computes pore volume of all cells in a grid.
    /// @param[in]  grid      a grid
    /// @param[in]  props     rock and fluid properties
    /// @param[out] porevol   the pore volume by cell.
    void computePorevolume(const UnstructuredGrid& grid,
			   const Opm::IncompPropertiesInterface& props,
			   std::vector<double>& porevol)
    {
	int num_cells = grid.number_of_cells;
	ASSERT(num_cells == props.numCells());
	porevol.resize(num_cells);
	const double* poro = props.porosity();
	std::transform(poro, poro + num_cells,
		       grid.cell_volumes,
		       porevol.begin(),
		       std::multiplies<double>());
    }


    /// @brief Computes pore volume of all cells in a grid, with rock compressibility effects.
    /// @param[in]  grid      a grid
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  rock_comp rock compressibility properties
    /// @param[in]  pressure  pressure by cell
    /// @param[out] porevol   the pore volume by cell.
    void computePorevolume(const UnstructuredGrid& grid,
                           const IncompPropertiesInterface& props,
                           const RockCompressibility& rock_comp,
                           const std::vector<double>& pressure,
                           std::vector<double>& porevol)
    {
        int num_cells = grid.number_of_cells;
        ASSERT(num_cells == props.numCells());
        porevol.resize(num_cells);
        const double* poro = props.porosity();
        for (int i = 0; i < num_cells; ++i) {
            porevol[i] = poro[i]*grid.cell_volumes[i]*rock_comp.poroMult(pressure[i]);
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
	    THROW("Sizes of s and pv vectors do not match.");
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
	    THROW("Sizes of s and pv vectors do not match.");
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
	    THROW("Sizes of s and src vectors do not match.");
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

        ASSERT (s.size() == nc * np);

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
            for (int w = 0; w < nw; ++w) {
                const double* zfrac = wells->zfrac + 3*w;
                for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                    const int perf_cell = wells->well_cells[perf];
                    double perf_rate = well_perfrates[perf];
                    if (perf_rate > 0.0) {
                        // perf_rate is a total inflow rate, we want a water rate.
                        if (wells->type[w] != INJECTOR) {
                            std::cout << "**** Warning: crossflow in well with index " << w << " ignored." << std::endl;
                            perf_rate = 0.0;
                        } else {
                            ASSERT(std::fabs(zfrac[WATER] + zfrac[OIL] - 1.0) < 1e-6);
                            perf_rate *= zfrac[WATER];
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
	const int dim = grid.dimensions;
	cell_velocity.clear();
	cell_velocity.resize(grid.number_of_cells*dim, 0.0);
	for (int face = 0; face < grid.number_of_faces; ++face) {
	    int c[2] = { grid.face_cells[2*face], grid.face_cells[2*face + 1] };
	    const double* fc = &grid.face_centroids[face*dim];
	    double flux = face_flux[face];
	    for (int i = 0; i < 2; ++i) {
		if (c[i] >= 0) {
		    const double* cc = &grid.cell_centroids[c[i]*dim];
		    for (int d = 0; d < dim; ++d) {
			double v_contrib = fc[d] - cc[d];
			v_contrib *= flux/grid.cell_volumes[c[i]];
			cell_velocity[c[i]*dim + d] += (i == 0) ? v_contrib : -v_contrib;
		    }
		}
	    }
	}
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
	src.resize(num_cells);
	for (int w = 0; w < wells.number_of_wells; ++w) {
	    if (wells.ctrls[w]->num != 1) {
		THROW("In wellsToSrc(): well has more than one control.");
	    }
	    if (wells.ctrls[w]->type[0] != RATE) {
		THROW("In wellsToSrc(): well is BHP, not RATE.");
	    }
	    if (wells.well_connpos[w+1] - wells.well_connpos[w] != 1) {
		THROW("In wellsToSrc(): well has multiple perforations.");
	    }
	    const double flow = wells.ctrls[w]->target[0];
	    const double cell = wells.well_cells[wells.well_connpos[w]];
	    src[cell] = (wells.type[w] == INJECTOR) ? flow : -flow;
	}
    }


    void computeWDP(const Wells& wells, const UnstructuredGrid& grid, const std::vector<double>& saturations,
                    const double* densities, std::vector<double>& wdp, bool per_grid_cell)
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
                wdp.push_back(density * (cell_depth - depth_ref));
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



} // namespace Opm
