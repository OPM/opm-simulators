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

class UnstructuredGrid;


namespace Opm
{

    class IncompPropertiesInterface;
    class Parameters;


    class TransportModelTwophase : public TransportModelInterface
    {
    public:
	TransportModelTwophase(const UnstructuredGrid* grid,
			       const Opm::IncompPropertiesInterface* props,
			       const double* darcyflux,
			       const double* porevolume,
			       const double* source,
			       const double dt,
			       double* saturation);

	virtual void solveSingleCell(int cell);

    private:
	const UnstructuredGrid* grid_;
	const IncompPropertiesInterface* props_;
	const double* darcyflux_;   /* one flux per face  in cdata::grid*/
	const double* porevolume_;  /* one volume per cell */
	const double* source_;      /* one source per cell */
	double dt_;
	double* saturation_;        /* one per cell */
	std::vector<double> fractionalflow_;  /* one per cell */

	Parameters getParameters(int cell);
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
