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

namespace Opm {

template<class Scalar>
class GasLiftWellState
{
public:
    GasLiftWellState(Scalar oil_rate,
                     Scalar oil_pot,
                     bool oil_is_limited,
                     Scalar gas_rate,
                     Scalar gas_pot,
                     bool gas_is_limited,
                     Scalar alq,
                     bool alq_is_limited,
                     Scalar water_rate,
                     Scalar water_pot,
                     bool water_is_limited,
                     Scalar bhp,
                     std::optional<bool> increase)
        : oil_rate_{oil_rate}
        , oil_pot_{oil_pot}
        , oil_is_limited_{oil_is_limited}
        , gas_rate_{gas_rate}
        , gas_pot_{gas_pot}
        , gas_is_limited_{gas_is_limited}
        , alq_{alq}
        , alq_is_limited_{alq_is_limited}
        , water_rate_{water_rate}
        , water_pot_{water_pot}
        , water_is_limited_{water_is_limited}
        , bhp_{bhp}
        , increase_{increase}
    {}

    Scalar alq() const { return alq_; }
    Scalar bhp() const { return bhp_; }
    bool alqChanged() { return increase_.has_value(); }
    bool alqIsLimited() const { return alq_is_limited_; }
    bool gasIsLimited() const { return gas_is_limited_; }
    Scalar gasRate() const { return gas_rate_; }
    Scalar gasPot() const { return gas_pot_; }
    std::pair<Scalar, Scalar> getRates() { return {oil_rate_, gas_rate_}; }
    std::optional<bool> increase() const { return increase_; }
    bool oilIsLimited() const { return oil_is_limited_; }
    Scalar oilRate() const { return oil_rate_; }
    Scalar waterRate() const { return water_rate_; }
    Scalar oilPot() const { return oil_pot_; }
    Scalar waterPot() const { return water_pot_; }
    bool waterIsLimited() const { return water_is_limited_; }
    void update(Scalar oil_rate,
                Scalar oil_pot,
                bool oil_is_limited,
                Scalar gas_rate,
                Scalar gas_pot,
                bool gas_is_limited,
                Scalar alq,
                bool alq_is_limited,
                Scalar water_rate,
                Scalar water_pot,
                Scalar water_is_limited,
                Scalar bhp,
                bool increase)
    {
        oil_rate_ = oil_rate;
        oil_pot_ = oil_pot;
        oil_is_limited_ = oil_is_limited;
        gas_rate_ = gas_rate;
        gas_pot_ = gas_pot;
        gas_is_limited_ = gas_is_limited;
        alq_ = alq;
        alq_is_limited_ = alq_is_limited;
        water_rate_ = water_rate;
        water_pot_ = water_pot;
        water_is_limited_ = water_is_limited;
        bhp_ = bhp;
        increase_ = increase;
    }

private:
    Scalar oil_rate_;
    Scalar oil_pot_;
    bool oil_is_limited_;
    Scalar gas_rate_;
    Scalar gas_pot_;
    bool gas_is_limited_;
    Scalar alq_;
    bool alq_is_limited_;
    Scalar water_rate_;
    Scalar water_pot_;
    bool water_is_limited_;
    Scalar bhp_;
    std::optional<bool> increase_;
};

} // namespace Opm

#endif // OPM_GASLIFT_WELL_STATE_HEADER_INCLUDED
