/*
  Copyright 2021 Equinor ASA.

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

#ifndef OPM_GASLIFT_WELL_STATE_HEADER_INCLUDED
#define OPM_GASLIFT_WELL_STATE_HEADER_INCLUDED

#include <optional>
#include <utility>

namespace Opm
{
    class GasLiftWellState
    {
    public:
        //GasLiftWellState() { }
        GasLiftWellState(double oil_rate, bool oil_is_limited,
                     double gas_rate, bool gas_is_limited,
            double alq, bool alq_is_limited, double water_rate, bool water_is_limited, std::optional<bool> increase) :
                oil_rate_{oil_rate},
                oil_is_limited_{oil_is_limited},
                gas_rate_{gas_rate},
                gas_is_limited_{gas_is_limited},
                alq_{alq},
                alq_is_limited_{alq_is_limited},
                water_rate_{water_rate},
                water_is_limited_{water_is_limited},
                increase_{increase}
        {}
        double alq() const { return alq_; }
        bool alqChanged() { return increase_.has_value(); }
        bool alqIsLimited() const { return alq_is_limited_; }
        bool gasIsLimited() const { return gas_is_limited_; }
        double gasRate() const { return gas_rate_; }
        std::pair<double, double> getRates() { return {oil_rate_, gas_rate_}; }
        std::optional<bool> increase() const { return increase_; }
        bool oilIsLimited() const { return oil_is_limited_; }
        double oilRate() const { return oil_rate_; }
        double waterRate() const { return water_rate_; }
        bool waterIsLimited() const { return water_is_limited_; }
        void update(double oil_rate, bool oil_is_limited,
            double gas_rate, bool gas_is_limited,
            double alq, bool alq_is_limited, double water_rate,
            double water_is_limited,
            bool increase)
        {
            oil_rate_ = oil_rate;
            oil_is_limited_ = oil_is_limited;
            gas_rate_ = gas_rate;
            gas_is_limited_ = gas_is_limited;
            alq_ = alq;
            alq_is_limited_ = alq_is_limited;
            water_rate_ = water_rate;
            water_is_limited_ = water_is_limited;
            increase_ = increase;
        }
    private:
        double oil_rate_;
        bool oil_is_limited_;
        double gas_rate_;
        bool gas_is_limited_;
        double alq_;
        bool alq_is_limited_;
        double water_rate_;
        bool water_is_limited_;
        std::optional<bool> increase_;
    };

} // namespace Opm

#endif // OPM_GASLIFT_WELL_STATE_HEADER_INCLUDED
