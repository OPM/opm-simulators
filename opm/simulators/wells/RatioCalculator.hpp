/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#ifndef RATIO_CALCULATOR_HPP
#define RATIO_CALCULATOR_HPP

#include <string>
#include <string_view>
#include <vector>

namespace Opm {

class DeferredLogger;

template<class Value>
class RatioCalculator
{
public:
    RatioCalculator(unsigned gasCompIdx,
                    unsigned oilCompIdx,
                    unsigned waterCompIdx,
                    std::string_view name);

    void gasOilVolumeRatio(Value& volumeRatio,
                           const Value& rv,
                           const Value& rs,
                           const Value& pressure,
                           const std::vector<Value>& cmix_s,
                           const std::vector<Value>& b_perfcells_dense,
                           DeferredLogger& deferred_logger) const;

private:
    unsigned gasComp_;
    unsigned oilComp_;
    unsigned waterComp_;
    std::string name_;
};

} // namespace Opm

#endif // RATIO_CALCULATOR_HPP
