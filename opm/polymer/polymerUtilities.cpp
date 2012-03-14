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


    /// @brief Computes injected and produced volumes of all phases,
    ///        and injeced and produced polymer mass.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// @param[in]  props     fluid and rock properties.
    /// @param[in]  polyprops polymer properties
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[in]  c         polymer concentration
    /// @param[in]  src       if < 0: total outflow, if > 0: first phase inflow.
    /// @param[in]  dt        timestep used
    /// @param[in]  inj_c     injected concentration
    /// @param[out] injected  must point to a valid array with P elements,
    ///                       where P = s.size()/src.size().
    /// @param[out] produced  must also point to a valid array with P elements.
    /// @param[out] polyinj   injected mass of polymer
    /// @param[out] polyprod  produced mass of polymer
    void computeInjectedProduced(const IncompPropertiesInterface& props,
                                 const Opm::PolymerProperties& polyprops,
				 const std::vector<double>& s,
				 const std::vector<double>& c,
				 const std::vector<double>& src,
				 const double dt,
                                 const double inj_c,
				 double* injected,
				 double* produced,
                                 double& polyinj,
                                 double& polyprod)
    {
        const int num_cells = src.size();
        const int np = s.size()/src.size();
        if (int(s.size()) != num_cells*np) {
            THROW("Sizes of s and src vectors do not match.");
        }
        std::fill(injected, injected + np, 0.0);
        std::fill(produced, produced + np, 0.0);
        polyinj = 0.0;
        polyprod = 0.0;
        std::vector<double> inv_eff_visc(np);
        const double* visc = props.viscosity();
        std::vector<double> mob(np);
        for (int cell = 0; cell < num_cells; ++cell) {
            if (src[cell] > 0.0) {
                injected[0] += src[cell]*dt;
                polyinj += src[cell]*dt*inj_c;
            } else if (src[cell] < 0.0) {
                const double flux = -src[cell]*dt;
                const double* sat = &s[np*cell];
                props.relperm(1, sat, &cell, &mob[0], 0);
                polyprops.effectiveInvVisc(c[cell], visc, &inv_eff_visc[0]);
                double totmob = 0.0;
                for (int p = 0; p < np; ++p) {
                    mob[p] *= inv_eff_visc[p];
                    totmob += mob[p];
                }
                for (int p = 0; p < np; ++p) {
                    produced[p] += (mob[p]/totmob)*flux;
                }
                polyprod += (mob[0]/totmob)*flux*c[cell]; // TODO check this term.
            }
        }
    }


    /// @brief Computes total polymer mass over all grid cells.
    /// @param[in]  pv        the pore volume by cell.
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[in]  c         polymer concentration
    /// @param[in]  dps       dead pore space
    /// @return               total polymer mass in grid.
    double computePolymerMass(const std::vector<double>& pv,
                              const std::vector<double>& s,
                              const std::vector<double>& c,
                              const double dps)
    {
	const int num_cells = pv.size();
	const int np = s.size()/pv.size();
	if (int(s.size()) != num_cells*np) {
	    THROW("Sizes of s and pv vectors do not match.");
	}
        double polymass = 0.0;
	for (int cell = 0; cell < num_cells; ++cell) {
            polymass += c[cell]*pv[cell]*(s[np*cell + 0] - dps);
	}
        return polymass;
    }



    /// @brief Computes total absorbed polymer mass over all grid cells.
    /// @param[in]  polyprops polymer properties
    /// @param[in]  pv        the pore volume by cell.
    /// @param[in]  cmax      max polymer concentration for cell
    /// @return               total absorbed polymer mass.
    double computePolymerAdsorbed(const Opm::PolymerProperties& polyprops,
                                  const std::vector<double>& pv,
                                  const std::vector<double>& cmax)
    {
	const int num_cells = pv.size();
        const double rhor = polyprops.rockDensity();
        double abs_mass = 0.0;
	for (int cell = 0; cell < num_cells; ++cell) {
            abs_mass += polyprops.adsorbtion(cmax[cell])*pv[cell]*rhor;
	}
        return abs_mass;
    }



} // namespace Opm

