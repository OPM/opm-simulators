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

#include <opm/polymer/polymerUtilities.hpp>

namespace Opm
{

    /// @brief Computes total mobility for a set of s/c values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  polyprops polymer properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[in]  c         polymer concentration
    /// @param[out] totmob    total mobilities.
    void computeTotalMobility(const Opm::IncompPropertiesInterface& props,
			      const Opm::PolymerProperties& polyprops,
			      const std::vector<int>& cells,
			      const std::vector<double>& s,
			      const std::vector<double>& c,
			      std::vector<double>& totmob)
    {
	int num_cells = cells.size();
	int num_phases = props.numPhases();
	totmob.resize(num_cells);
	ASSERT(int(s.size()) == num_cells*num_phases);
	std::vector<double> kr(num_cells*num_phases);
	props.relperm(num_cells, &s[0], &cells[0], &kr[0], 0);
	const double* mu = props.viscosity();
	double inv_mu_eff[2] = { 0.0 };
	for (int cell = 0; cell < num_cells; ++cell) {
	    totmob[cell] = 0;
	    polyprops.effectiveInvVisc(c[cell], mu, inv_mu_eff);
	    for (int phase = 0; phase < num_phases; ++phase) {	
		totmob[cell] += kr[num_phases*cell + phase]*inv_mu_eff[phase];
	    }
	}
    }




    /// @brief Computes total mobility and omega for a set of s/c values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  polyprops polymer properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  s         saturation values (for all phases)
    /// @param[in]  c         polymer concentration
    /// @param[out] totmob    total mobility
    /// @param[out] omega     mobility-weighted (or fractional-flow weighted)
    ///                       fluid densities.
    void computeTotalMobilityOmega(const Opm::IncompPropertiesInterface& props,
				   const Opm::PolymerProperties& polyprops,
				   const std::vector<int>& cells,
				   const std::vector<double>& s,
				   const std::vector<double>& c,
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
	double inv_mu_eff[2] = { 0.0 };
	const double* rho = props.density();
	for (int cell = 0; cell < num_cells; ++cell) {
	    totmob[cell] = 0.0;
	    omega[cell] = 0.0;
	    polyprops.effectiveInvVisc(c[cell], mu, inv_mu_eff);
	    for (int phase = 0; phase < num_phases; ++phase) {	
		totmob[cell] += kr[num_phases*cell + phase]*inv_mu_eff[phase];
	    }
	    // Must finish computing totmob before we can use it.
	    for (int phase = 0; phase < num_phases; ++phase) {	
		omega[cell] += rho[phase]*(kr[num_phases*cell + phase]*inv_mu_eff[phase])/totmob[cell];
	    }
	}
    }


} // namespace Opm

