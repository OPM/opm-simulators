/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H


#include <opm/grid/UnstructuredGrid.h>
#include <opm/autodiff/SimulatorSequentialBlackoil.hpp>
#include <opm/autodiff/FlowMainSequential.hpp>
#include <opm/autodiff/BlackoilPressureModel.hpp>
#include <opm/autodiff/BlackoilReorderingTransportModel.hpp>
#include <opm/autodiff/StandardWells.hpp>


// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    typedef UnstructuredGrid Grid;
    typedef Opm::StandardWells WellModel;
    typedef Opm::SimulatorSequentialBlackoil<Grid, WellModel, Opm::BlackoilPressureModel, Opm::BlackoilReorderingTransportModel> Simulator;

    Opm::FlowMainSequential<Grid, Simulator> mainfunc;
    return mainfunc.execute(argc, argv);
}
