/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Opm::FluidHeatConductionParams
 */
#ifndef OPM_FLUID_HEAT_CONDUCTION_PARAMS_HPP
#define OPM_FLUID_HEAT_CONDUCTION_PARAMS_HPP

namespace Opm {
/*!
 * \brief Parameters for the heat conduction law which just takes the conductivity of a given fluid phase.
 */
template <class ScalarT>
class FluidHeatConductionParams
{
    // do not copy!
    FluidHeatConductionParams(const FluidHeatConductionParams &)
    {}

public:
    typedef ScalarT Scalar;

    FluidHeatConductionParams()
    { }

};

} // namespace Opm

#endif
