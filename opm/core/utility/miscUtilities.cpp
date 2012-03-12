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
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <functional>

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


    /// @brief Computes average saturations over all grid cells.
    /// @param[out] pv        the pore volume by cell.
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
	int num_cells = cells.size();
	int num_phases = props.numPhases();
	totmob.resize(num_cells);
	ASSERT(int(s.size()) == num_cells*num_phases);
	std::vector<double> kr(num_cells*num_phases);
	props.relperm(num_cells, &s[0], &cells[0], &kr[0], 0);
	const double* mu = props.viscosity();
	for (int cell = 0; cell < num_cells; ++cell) {
	    totmob[cell] = 0;
	    for (int phase = 0; phase < num_phases; ++phase) {	
		totmob[cell] += kr[num_phases*cell + phase]/mu[phase];
	    }
	}
    }


    /// @brief Computes total mobility and omega for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobility
    /// @param[out] omega     mobility-weighted (or fractional-flow weighted)
    ///                       fluid densities.
    void computeTotalMobilityOmega(const Opm::IncompPropertiesInterface& props,
				   const std::vector<int>& cells,
				   const std::vector<double>& s,
				   std::vector<double>& totmob,
				   std::vector<double>& omega)
    {
	int num_cells = cells.size();
	int num_phases = props.numPhases();
	totmob.resize(num_cells);
	omega.resize(num_cells);
	ASSERT(int(s.size()) == num_cells*num_phases);
	std::vector<double> kr(num_cells*num_phases);
	props.relperm(num_cells, &s[0], &cells[0], &kr[0], 0);
	const double* mu = props.viscosity();
	for (int cell = 0; cell < num_cells; ++cell) {
	    totmob[cell] = 0.0;
	    for (int phase = 0; phase < num_phases; ++phase) {	
		totmob[cell] += kr[num_phases*cell + phase]/mu[phase];
	    }
	}
	const double* rho = props.density();
	for (int cell = 0; cell < num_cells; ++cell) {
	    omega[cell] = 0.0;
	    for (int phase = 0; phase < num_phases; ++phase) {	
		omega[cell] += rho[phase]*(kr[num_phases*cell + phase]/mu[phase])/totmob[cell];
	    }
	}
    }

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




} // namespace Opm
