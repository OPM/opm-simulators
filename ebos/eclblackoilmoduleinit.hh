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

#ifndef ECL_BLACKOILMODULE_INIT_HH
#define ECL_BLACKOILMODULE_INIT_HH

namespace Opm {

class EclipseState;
template<class Scalar> struct BlackOilBrineParams;
template<class Scalar> struct BlackOilExtboParams;

//! \brief Setup parameters for brine module from an EclipseState.
template<bool enableSaltPrecipitation, class Scalar>
BlackOilBrineParams<Scalar> setupBrineParams(bool enableBrine,
                                             const EclipseState& eclState);

//! \brief Setup parameters for extbo module from an EclipseState.
template<class Scalar>
BlackOilExtboParams<Scalar> setupExtboParams(bool enableExtbo,
                                             const EclipseState& eclState);
}

#endif
