/*
  Copyright 2023 SINTEF.

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

#include <config.h>

#include <opm/simulators/flow/priVarsPacking.hpp>

#define BOOST_TEST_MODULE priVarsPacking
#include <boost/test/unit_test.hpp>


// Must define a class for testing, using extracts from BlackoilPrimaryVariables,
// but without the typetags.
class PriVarMeaning
{
public:

    enum class WaterMeaning {
        Sw,  // water saturation
        Rvw, // vaporized water
        Rsw, // dissolved gas in water
        Disabled, // The primary variable is not used
    };

    enum class PressureMeaning {
        Po, // oil pressure
        Pg, // gas pressure
        Pw, // water pressure
    };
    enum class GasMeaning {
        Sg, // gas saturation
        Rs, // dissolved gas in oil
        Rv, // vapporized oil
        Disabled, // The primary variable is not used
    };

    enum class BrineMeaning {
        Cs, // salt concentration
        Sp, // (precipitated) salt saturation
        Disabled, // The primary variable is not used
    };

    enum class SolventMeaning {
        Ss, // solvent saturation
        Rsolw, // dissolved solvent in water
        Disabled, // The primary variable is not used
    };

    WaterMeaning primaryVarsMeaningWater() const
    { return primaryVarsMeaningWater_; }
    void setPrimaryVarsMeaningWater(WaterMeaning newMeaning)
    { primaryVarsMeaningWater_ = newMeaning; }

    PressureMeaning primaryVarsMeaningPressure() const
    { return primaryVarsMeaningPressure_; }
    void setPrimaryVarsMeaningPressure(PressureMeaning newMeaning)
    { primaryVarsMeaningPressure_ = newMeaning; }

    GasMeaning primaryVarsMeaningGas() const
    { return primaryVarsMeaningGas_; }
    void setPrimaryVarsMeaningGas(GasMeaning newMeaning)
    { primaryVarsMeaningGas_ = newMeaning; }

    BrineMeaning primaryVarsMeaningBrine() const
    { return primaryVarsMeaningBrine_; }
    void setPrimaryVarsMeaningBrine(BrineMeaning newMeaning)
    { primaryVarsMeaningBrine_ = newMeaning; }

    SolventMeaning primaryVarsMeaningSolvent() const
    { return primaryVarsMeaningSolvent_; }
    void setPrimaryVarsMeaningSolvent(SolventMeaning newMeaning)
    { primaryVarsMeaningSolvent_ = newMeaning; }

    bool operator==(const PriVarMeaning& other) const
    {
        return primaryVarsMeaningWater_ == other.primaryVarsMeaningWater_
            && primaryVarsMeaningPressure_ == other.primaryVarsMeaningPressure_
            && primaryVarsMeaningGas_ == other.primaryVarsMeaningGas_
            && primaryVarsMeaningBrine_ == other.primaryVarsMeaningBrine_
            && primaryVarsMeaningSolvent_ == other.primaryVarsMeaningSolvent_;
    }
private:
    WaterMeaning primaryVarsMeaningWater_ = WaterMeaning::Disabled;
    PressureMeaning primaryVarsMeaningPressure_ = PressureMeaning::Pw;
    GasMeaning primaryVarsMeaningGas_ = GasMeaning::Disabled;
    BrineMeaning primaryVarsMeaningBrine_ = BrineMeaning::Disabled;
    SolventMeaning primaryVarsMeaningSolvent_ = SolventMeaning::Disabled;
};



BOOST_AUTO_TEST_CASE(meanings)
{
    // Test default meanings.
    {
        PriVarMeaning pv1, pv2;
        std::size_t p = Opm::PVUtil::pack(pv1);
        Opm::PVUtil::unPack(pv2, p);
        BOOST_CHECK(pv1 == pv2);
    }

    // Test explicitly set meanings.
    {
        PriVarMeaning pv1, pv2, pv3;
        pv1.setPrimaryVarsMeaningPressure(PriVarMeaning::PressureMeaning::Pw);
        pv1.setPrimaryVarsMeaningWater(PriVarMeaning::WaterMeaning::Rvw);
        pv1.setPrimaryVarsMeaningGas(PriVarMeaning::GasMeaning::Disabled);
        pv1.setPrimaryVarsMeaningBrine(PriVarMeaning::BrineMeaning::Cs);
        pv1.setPrimaryVarsMeaningSolvent(PriVarMeaning::SolventMeaning::Disabled);
        std::size_t p = Opm::PVUtil::pack(pv1);
        Opm::PVUtil::unPack(pv2, p);
        BOOST_CHECK(pv1 == pv2);
        BOOST_CHECK(!(pv1 == pv3));
        BOOST_CHECK(!(pv2 == pv3));
    }
}
