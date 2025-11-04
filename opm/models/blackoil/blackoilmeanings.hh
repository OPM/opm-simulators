/*
  Copyright 2025 Equinor ASA

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
#ifndef OPM_BLACK_OIL_MEANINGS_HH
#define OPM_BLACK_OIL_MEANINGS_HH
#include <cstdint>

namespace Opm::BlackOil
{
enum class WaterMeaning : std::uint8_t {
    Sw, // water saturation
    Rvw, // vaporized water
    Rsw, // dissolved gas in water
    Disabled, // The primary variable is not used
};

enum class PressureMeaning : std::uint8_t {
    Po, // oil pressure
    Pg, // gas pressure
    Pw, // water pressure
};

enum class GasMeaning : std::uint8_t {
    Sg, // gas saturation
    Rs, // dissolved gas in oil
    Rv, // vapporized oil
    Disabled, // The primary variable is not used
};

enum class BrineMeaning : std::uint8_t {
    Cs, // salt concentration
    Sp, // (precipitated) salt saturation
    Disabled, // The primary variable is not used
};

enum class SolventMeaning : std::uint8_t {
    Ss, // solvent saturation
    Rsolw, // dissolved solvent in water
    Disabled, // The primary variable is not used
};
} // namespace Opm::BlackOil
#endif // OPM_BLACK_OIL_MEANINGS_HH
