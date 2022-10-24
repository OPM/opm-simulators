/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BLACKOILWELLMODEL_GUIDE_RATES_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_GUIDE_RATES_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>

#include <string>

namespace Opm {

class BlackoilWellModelGeneric;
namespace data { class GuideRateValue; }
class Group;
class Well;

/// Class for handling the guide rates in the blackoil well model.
class BlackoilWellModelGuideRates
{
public:
    //! \brief Constructor initializes reference to the well model.
    BlackoilWellModelGuideRates(const BlackoilWellModelGeneric& wellModel)
        : wellModel_(wellModel)
    {}

    //! \brief Obtain guide rate values.
    void getGuideRateValues(const GuideRate::RateVector& qs,
                            const bool                   is_inj,
                            const std::string&           wgname,
                            data::GuideRateValue&        grval) const;

    //! \brief Obtain guide rate values for well.
    data::GuideRateValue getGuideRateValues(const Well& well) const;

private:
    const BlackoilWellModelGeneric& wellModel_; //!< Reference to well model
};


} // namespace Opm

#endif
