/*
  Copyright 2013 University of Stuttgart
  Copyright 2015 Andreas Lauser
  Copyright 2020 NORCE

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

#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/components/co2tables.inc>

TabulatedFunction CO2Tables::tabulatedEnthalpy;
TabulatedFunction CO2Tables::tabulatedDensity;
const double CO2Tables::brineSalinity = 1.000000000000000e-01;

// initialize the static tables once. this is a bit hacky in so far as it uses some
// advanced C++ features (static initializer functions)
int initCO2Tables_()
{
    CO2Tables::tabulatedEnthalpy.resize(TabulatedEnthalpyTraits::xMin,
                                        TabulatedEnthalpyTraits::xMax,
                                        TabulatedEnthalpyTraits::numX,
                                        TabulatedEnthalpyTraits::yMin,
                                        TabulatedEnthalpyTraits::yMax,
                                        TabulatedEnthalpyTraits::numY);

    for (unsigned i = 0; i < TabulatedEnthalpyTraits::numX; ++i)
        for (unsigned j = 0; j < TabulatedEnthalpyTraits::numY; ++j)
            CO2Tables::tabulatedEnthalpy.setSamplePoint(i, j, TabulatedEnthalpyTraits::vals[i][j]);

    CO2Tables::tabulatedDensity.resize(TabulatedDensityTraits::xMin,
                                       TabulatedDensityTraits::xMax,
                                       TabulatedDensityTraits::numX,
                                       TabulatedDensityTraits::yMin,
                                       TabulatedDensityTraits::yMax,
                                       TabulatedDensityTraits::numY);

    for (unsigned i = 0; i < TabulatedDensityTraits::numX; ++i)
        for (unsigned j = 0; j < TabulatedDensityTraits::numY; ++j)
            CO2Tables::tabulatedDensity.setSamplePoint(i, j, TabulatedDensityTraits::vals[i][j]);

    return 0;
}

extern int co2TablesInitDummy__;
int co2TablesInitDummy__ = initCO2Tables_();
