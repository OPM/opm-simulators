// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/models/blackoil/blackoilnewtonmethodparameters.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

template<class Scalar>
void BlackoilNewtonParams<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::DpMaxRel<Scalar>>
        ("Maximum relative change of pressure in a single iteration");
    Parameters::Register<Parameters::DsMax<Scalar>>
        ("Maximum absolute change of any saturation in a single iteration");
    Parameters::Register<Parameters::PriVarOscilationThreshold<Scalar>>
        ("The threshold value for the primary variable switching conditions "
         "after its meaning has switched to hinder oscillations");
    Parameters::Register<Parameters::ProjectSaturations>
        ("Option for doing saturation projection");
    Parameters::Register<Parameters::MaxTemperatureChange<Scalar>>
        ("Maximum absolute change of temperature in a single iteration");
    Parameters::Register<Parameters::TemperatureMax<Scalar>>
        ("Maximum absolute temperature");
    Parameters::Register<Parameters::TemperatureMin<Scalar>>
        ("Minimum absolute temperature");
    Parameters::Register<Parameters::PressureMax<Scalar>>
        ("Maximum absolute pressure");
    Parameters::Register<Parameters::PressureMin<Scalar>>
        ("Minimum absolute pressure");
    Parameters::Register<Parameters::MaximumWaterSaturation<Scalar>>
        ("Maximum water saturation");
    Parameters::Register<Parameters::WaterOnlyThreshold<Scalar>>
        ("Cells with water saturation above or equal is considered one-phase water only");
}

template<class Scalar>
void BlackoilNewtonParams<Scalar>::read()
{
    priVarOscilationThreshold_ = Parameters::Get<Parameters::PriVarOscilationThreshold<Scalar>>();
    dpMaxRel_ = Parameters::Get<Parameters::DpMaxRel<Scalar>>();
    dsMax_ = Parameters::Get<Parameters::DsMax<Scalar>>();
    projectSaturations_ = Parameters::Get<Parameters::ProjectSaturations>();
    maxTempChange_ = Parameters::Get<Parameters::MaxTemperatureChange<Scalar>>();
    tempMax_ = Parameters::Get<Parameters::TemperatureMax<Scalar>>();
    tempMin_ = Parameters::Get<Parameters::TemperatureMin<Scalar>>();
    pressMax_ = Parameters::Get<Parameters::PressureMax<Scalar>>();
    pressMin_ = Parameters::Get<Parameters::PressureMin<Scalar>>();
    waterSaturationMax_ = Parameters::Get<Parameters::MaximumWaterSaturation<Scalar>>();
    waterOnlyThreshold_ = Parameters::Get<Parameters::WaterOnlyThreshold<Scalar>>();
}

template struct BlackoilNewtonParams<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct BlackoilNewtonParams<float>;
#endif

} // namespace Opm
