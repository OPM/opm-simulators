/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015 IRIS AS

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

#ifndef OPM_BLACKOILMODELEBOSTYPETAGS_HEADER_INCLUDED
#define OPM_BLACKOILMODELEBOSTYPETAGS_HEADER_INCLUDED

#include <ebos/eclproblem.hh>

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem));
SET_BOOL_PROP(EclFlowProblem, DisableWells, true);
SET_BOOL_PROP(EclFlowProblem, EnableDebuggingChecks, false);

// SWATINIT is done by the flow part of flow_ebos. this can be removed once the legacy
// code for fluid and satfunc handling gets fully retired.
SET_BOOL_PROP(EclFlowProblem, EnableSwatinit, false);
}}

#endif
