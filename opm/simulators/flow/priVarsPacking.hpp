/*
  Copyright 2023 Total SE

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

#ifndef OPM_PRIVARSPACKING_HEADER_INCLUDED
#define OPM_PRIVARSPACKING_HEADER_INCLUDED

#include <cstddef>

namespace Opm {

    namespace PVUtil {
        constexpr int fbits = 4;

        template <class PV>
        std::size_t pack(const PV& privar) {
            std::size_t m1 = static_cast<std::size_t>(privar.primaryVarsMeaningWater());
            std::size_t m2 = static_cast<std::size_t>(privar.primaryVarsMeaningPressure());
            std::size_t m3 = static_cast<std::size_t>(privar.primaryVarsMeaningGas());
            std::size_t m4 = static_cast<std::size_t>(privar.primaryVarsMeaningBrine());
            std::size_t m5 = static_cast<std::size_t>(privar.primaryVarsMeaningSolvent());
            return m1 + (m2 << fbits*1) + (m3 << fbits*2) + (m4 << fbits*3) + (m5 << fbits*4);
        }

        template <class PV>
        void unPack(PV& privar, const std::size_t meanings) {
            const std::size_t filter = ((1 << fbits) - 1);
            std::size_t m1 = (meanings >> fbits*0) & filter;
            std::size_t m2 = (meanings >> fbits*1) & filter;
            std::size_t m3 = (meanings >> fbits*2) & filter;
            std::size_t m4 = (meanings >> fbits*3) & filter;
            std::size_t m5 = (meanings >> fbits*4) & filter;
            privar.setPrimaryVarsMeaningWater(typename PV::WaterMeaning(m1));
            privar.setPrimaryVarsMeaningPressure(typename PV::PressureMeaning(m2));
            privar.setPrimaryVarsMeaningGas(typename PV::GasMeaning(m3));
            privar.setPrimaryVarsMeaningBrine(typename PV::BrineMeaning(m4));
            privar.setPrimaryVarsMeaningSolvent(typename PV::SolventMeaning(m5));
        }
    } // namespace PVMeanings
} // namespace Opm

#endif // OPM_PRIVARSPACKING_HEADER_INCLUDED
