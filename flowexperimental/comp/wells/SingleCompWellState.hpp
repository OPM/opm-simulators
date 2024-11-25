/*
Copyright 2024, SINTEF Digital

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

#ifndef OPM_SINGLE_COMP_WELL_STATE_HPP
#define OPM_SINGLE_COMP_WELL_STATE_HPP

#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <string>
#include <vector>

namespace Opm {

template <class Scalar>
class ConnectionData {
public:
   Scalar pressure {};
   std::vector<Scalar> phase_rates; // surface rates
   std::vector<Scalar> reservoir_rates; // phase rates
   std::vector<Scalar> total_molar_fractions;
};

template <class Scalar>
class SingleCompWellState {
public:
    std::string name;
    bool producer;

    WellStatus status{WellStatus::OPEN};
    Scalar bhp{0};
    Scalar temperature{0};
    std::vector<Scalar> surface_phase_rates;
    std::vector<Scalar> phase_fractions; // V or L
    std::vector<Scalar> reservoir_phase_rates;
    // WZMF
    std::vector<Scalar> total_molar_fractions;
    // WXMF WYMF and WAMF
    // TODO: potentiall we can use std::vector<std::vector<Sclar>>
    // TODO: to do phase_molar_fractions
    std::vector<std::vector<Scalar> > phase_molar_fractions;

    std::vector<ConnectionData<Scalar>> connection_data;

};

} // namespace Opm

#endif // OPM_SINGLE_COMP_WELL_STATE_HPP
