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

#ifndef OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
#define OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED

#include <opm/core/transport/reorder/TransportModelInterface.hpp>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    class TransportModelTwophase : public TransportModelInterface
    {
    public:
	TransportModelTwophase(const UnstructuredGrid& grid,
			       const double* porevolume,
			       const Opm::IncompPropertiesInterface& props,
			       const double tol,
			       const int maxit);

	void solve(const double* darcyflux,
		   const double* source,
		   const double dt,
		   double* saturation);

	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);

    private:
	const UnstructuredGrid& grid_;
	const double* porevolume_;  // one volume per cell
	const IncompPropertiesInterface& props_;
	const double* visc_;
	std::vector<double> smin_;
	std::vector<double> smax_;
	double tol_;
	double maxit_;

	const double* darcyflux_;   // one flux per grid face
	const double* source_;      // one source per cell
	double dt_;
	double* saturation_;        // one per cell
	std::vector<double> fractionalflow_;  // one per cell

	// Storing the upwind and downwind graphs for experiments.
	std::vector<int> ia_upw_;
	std::vector<int> ja_upw_;
	std::vector<int> ia_downw_;
	std::vector<int> ja_downw_;

	struct Residual;
	double fracFlow(double s, int cell) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
