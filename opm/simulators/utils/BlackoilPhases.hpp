/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILPHASES_HEADER_INCLUDED
#define OPM_BLACKOILPHASES_HEADER_INCLUDED

#include <array>
#include <vector>

namespace Opm
{

    class BlackoilPhases
    {
    public:
        static const int MaxNumPhases = 3;

        // "Crypto phases" are "phases" (or rather "conservation quantities") in the
        // sense that they can be active or not and canonical indices can be translated
        // to and from active ones. That said, they are not considered by num_phases or
        // MaxNumPhases. The crypto phases which are currently implemented are solvent,
        // polymer, energy, polymer molecular weight, foam and brine.
        static const int NumCryptoPhases = 7;

        // enum ComponentIndex { Water = 0, Oil = 1, Gas = 2 };
        enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2, Solvent = 3, Polymer = 4, Energy = 5, PolymerMW = 6, Foam = 7, Brine = 8, ZFraction = 9 };
    };

    struct PhaseUsage : public BlackoilPhases
    {
        PhaseUsage() = default;
        explicit PhaseUsage(std::vector<BlackoilPhases::PhaseIndex> phases);


        std::array<int, MaxNumPhases + NumCryptoPhases> phase_used;
        std::array<int, MaxNumPhases + NumCryptoPhases> phase_pos;

        int num_phases;
        bool has_solvent{};
        bool has_polymer{};
        bool has_energy{};
        // polymer molecular weight
        bool has_polymermw{};
        bool has_foam{};
        bool has_brine{};
        bool has_zFraction{};

    };

} // namespace Opm

#endif // OPM_BLACKOILPHASES_HEADER_INCLUDED
