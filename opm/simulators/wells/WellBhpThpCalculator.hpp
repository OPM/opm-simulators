/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED
#define OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Opm
{

class SummaryState;
class WellInterfaceGeneric;

//! \brief Class for computing BHP limits.
class WellBhpThpCalculator {
public:
    //! \brief Constructor sets reference to well.
    WellBhpThpCalculator(const WellInterfaceGeneric& well) : well_(well) {}

    //! \brief Checks if well has THP constraints.
    bool wellHasTHPConstraints(const SummaryState& summaryState) const;

private:
    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_BHP_THP_CALCULATOR_HEADER_INCLUDED
