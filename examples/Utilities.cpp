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

#include "Utilities.hpp"

#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/cart_grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/cpgpreprocess/cgridinterface.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <tr1/array>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>




namespace Opm
{






    void
    compute_porevolume(const UnstructuredGrid* g,
		       const Opm::IncompPropertiesInterface& props,
		       std::vector<double>& porevol)
    {
	int num_cells = g->number_of_cells;
	porevol.resize(num_cells);
	const double* poro = props.porosity();
	::std::transform(poro, poro + num_cells,
			 g->cell_volumes,
			 porevol.begin(),
			 ::std::multiplies<double>());
    }


    void
    compute_totmob(const Opm::IncompPropertiesInterface& props,
		   const std::vector<double>& s,
		   std::vector<double>& totmob)
    {
	int num_cells = props.numCells();
	int num_phases = props.numPhases();
	totmob.resize(num_cells);
	ASSERT(int(s.size()) == num_cells*num_phases);
	std::vector<int> cells(num_cells);
	for (int cell = 0; cell < num_cells; ++cell) {
	    cells[cell] = cell;
	}
	std::vector<double> kr(num_cells*num_phases);
	props.relperm(num_cells, &s[0], &cells[0], &kr[0], 0);
	const double* mu = props.viscosity();
	for (int cell = 0; cell < num_cells; ++cell) {
	    totmob[cell] = 0;
	    for (int phase = 0; phase < num_phases; ++phase) {	
		totmob[cell] += kr[2*cell + phase]/mu[phase];
	    }
	}
    }





    void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw)
    {
	int num = sboth.size()/2;
	sw.resize(num);
	for (int i = 0; i < num; ++i) {
	    sw[i] = sboth[2*i];
	}
    }

    void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth)
    {
	int num = sw.size();
	sboth.resize(2*num);
	for (int i = 0; i < num; ++i) {
	    sboth[2*i] = sw[i];
	    sboth[2*i + 1] = 1.0 - sw[i];
	}
    }



} // namespace Opm

