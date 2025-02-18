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

#define BOOST_TEST_MODULE TestBlackOilFluidSystemGPU
#include <config.h>

#include <boost/test/unit_test.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystemDynamic.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

BOOST_AUTO_TEST_CASE(TestThrowOnUnInitialized) {
    BOOST_CHECK(!Opm::BlackOilFluidSystem<double>::isInitialized());
    BOOST_CHECK_THROW(Opm::BlackOilFluidSystem<double>::getDynamicInstance(), std::logic_error);
}

BOOST_AUTO_TEST_CASE(TestDynamicCreationCPU) {
    Opm::BlackOilFluidSystem<double>::initBegin(5);
    Opm::BlackOilFluidSystem<double>::initEnd();
    BOOST_CHECK(Opm::BlackOilFluidSystem<double>::isInitialized());
    auto staticDummyInstance = Opm::BlackOilFluidSystem<double>{};
    BOOST_CHECK(staticDummyInstance.isInitialized());

    Opm::BlackOilFluidSystemDynamic<double>& fluidSystem = Opm::BlackOilFluidSystem<double>::getDynamicInstance();
}
BOOST_AUTO_TEST_CASE(TestDynamicGettersCPU) {
    Opm::BlackOilFluidSystemDynamic<double>& fluidSystem = Opm::BlackOilFluidSystem<double>::getDynamicInstance();
    BOOST_CHECK_EQUAL(fluidSystem.numActivePhases(), Opm::BlackOilFluidSystem<double>::numActivePhases());
    for (int phase = 0; phase < fluidSystem.numActivePhases(); ++phase) {
      BOOST_CHECK_EQUAL(fluidSystem.phaseIsActive(phase), Opm::BlackOilFluidSystem<double>::phaseIsActive(phase));
    }
    //BOOST_CHECK_EQUAL(fluidSystem.surfaceTemperature(), Opm::BlackOilFluidSystem<double>::surfaceTemperature());
    //BOOST_CHECK_EQUAL(fluidSystem.surfacePressure(), Opm::BlackOilFluidSystem<double>::surfacePressure());
    BOOST_CHECK_EQUAL(fluidSystem.enableDissolvedGas(), Opm::BlackOilFluidSystem<double>::enableDissolvedGas());
    BOOST_CHECK_EQUAL(fluidSystem.enableDissolvedGasInWater(), Opm::BlackOilFluidSystem<double>::enableDissolvedGasInWater());
    BOOST_CHECK_EQUAL(fluidSystem.enableVaporizedOil(), Opm::BlackOilFluidSystem<double>::enableVaporizedOil());
    BOOST_CHECK_EQUAL(fluidSystem.enableVaporizedWater(), Opm::BlackOilFluidSystem<double>::enableVaporizedWater());
    BOOST_CHECK_EQUAL(fluidSystem.enableDiffusion(), Opm::BlackOilFluidSystem<double>::enableDiffusion());
    BOOST_CHECK_EQUAL(fluidSystem.isInitialized(), Opm::BlackOilFluidSystem<double>::isInitialized());
    BOOST_CHECK_EQUAL(fluidSystem.useSaturatedTables(), Opm::BlackOilFluidSystem<double>::useSaturatedTables());
    //BOOST_CHECK_EQUAL(fluidSystem.enthalpyEqEnergy(), Opm::BlackOilFluidSystem<double>::enthalpyEqEnergy());
}

