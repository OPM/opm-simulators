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

#ifndef OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED
#define OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED

#include <opm/core/transport/reorder/TransportModelInterface.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <vector>

class UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;


    /// Containing all the extra information needed to model
    /// polymer-affected flow behaviour. This as an alternative
    /// to changing the IncompPropertiesInterface class.
    /// \TODO Improve encapsulation.
    struct PolymerData
    {
	double c_max_limit;
	double omega;
	double viscMult(double c) const
	{
	    return Opm::linearInterpolation(c_vals_visc, visc_mult_vals, c);
	}
	double rhor;
	double dps;
	double adsorbtion(double c) const
	{
	    return Opm::linearInterpolation(c_vals_ads, ads_vals, c);
	}

	std::vector<double> c_vals_visc;
	std::vector<double> visc_mult_vals;
	std::vector<double> c_vals_ads;
	std::vector<double> ads_vals;
    };



    /// A transport model for two-phase flow with polymer in the
    /// water phase.
    /// \TODO Include permeability reduction effect.
    class TransportModelPolymer : public TransportModelInterface
    {
    public:
	TransportModelPolymer(const UnstructuredGrid& grid,
			      const double* porosity,
			      const double* porevolume,
			      const IncompPropertiesInterface& props,
			      const PolymerData& polyprops);

	/// Solve transport eqn with implicit Euler scheme, reordered.
	/// \TODO Now saturation is expected to be one sw value per cell,
	/// change to [sw so] per cell.
	void solve(const double* darcyflux,
		   const double* source,
		   const double dt,
		   const double inflow_c,
		   double* saturation,
		   double* concentration,
		   double* cmax);

	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);

    private:
	const UnstructuredGrid& grid_;
	const double* porosity_;
	const double* porevolume_;  // one volume per cell
	const IncompPropertiesInterface& props_;
	const PolymerData& polyprops_;

	const double* darcyflux_;   // one flux per grid face
	const double* source_;      // one source per cell
	double dt_;
	double inflow_c_;
	double* saturation_;        // one per cell
	double* concentration_;
	double* cmax_;
	std::vector<double> fractionalflow_;  // one per cell
	std::vector<double> mc_;  // one per cell
	const double* visc_;

	struct ResidualC;
	struct ResidualS;
	double fracFlow(double s, double c, int cell) const;
	double computeMc(double c) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED
