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
#include <map>

struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    class TransportModelTwophase : public TransportModelInterface
    {
    public:
	TransportModelTwophase(const UnstructuredGrid& grid,
			       const Opm::IncompPropertiesInterface& props,
			       const double tol,
			       const int maxit);

	void solve(const double* darcyflux,
                   const double* porevolume,
		   const double* source,
		   const double dt,
		   double* saturation);

	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);

        void initGravity(const double* grav);
        void solveSingleCellGravity(const std::vector<int>& cells,
                                    const int pos,
                                    const double* gravflux);
        int solveGravityColumn(const std::vector<int>& cells);
        void solveGravity(const std::vector<std::vector<int> >& columns,
                          const double* porevolume,
                          const double dt,
                          std::vector<double>& saturation);

    private:
	const UnstructuredGrid& grid_;
	const IncompPropertiesInterface& props_;
	const double* visc_;
	std::vector<double> smin_;
	std::vector<double> smax_;
	double tol_;
	double maxit_;

	const double* darcyflux_;   // one flux per grid face
	const double* porevolume_;  // one volume per cell
	const double* source_;      // one source per cell
	double dt_;
	double* saturation_;        // one per cell
	std::vector<double> fractionalflow_;  // = m[0]/(m[0] + m[1]) per cell
        // For gravity segregation.
        std::vector<double> gravflux_;
        std::vector<double> mob_;
        std::vector<double> s0_;

	// Storing the upwind and downwind graphs for experiments.
	std::vector<int> ia_upw_;
	std::vector<int> ja_upw_;
	std::vector<int> ia_downw_;
	std::vector<int> ja_downw_;

	struct Residual;
	double fracFlow(double s, int cell) const;

        struct GravityResidual;
        void mobility(double s, int cell, double* mob) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
