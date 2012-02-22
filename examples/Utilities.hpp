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

#ifndef OPM_UTILITIES_HEADER_INCLUDED
#define OPM_UTILITIES_HEADER_INCLUDED

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
		       std::vector<double>& porevol);


    void
    compute_totmob(const Opm::IncompPropertiesInterface& props,
		   const std::vector<double>& s,
		   std::vector<double>& totmob);


    void
    compute_totmob_omega(const Opm::IncompPropertiesInterface& props,
			 const std::vector<double>& s,
			 std::vector<double>& totmob,
			 std::vector<double>& omega);


    void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw);

    void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth);



} // namespace Opm




#endif // OPM_UTILITIES_HEADER_INCLUDED
