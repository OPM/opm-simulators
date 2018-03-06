// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This test makes sure that the programming interface is
 *        observed by all fluid systems
 */
#include "config.h"

#include <opm/material/checkFluidSystem.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

// include all fluid systems in opm-material
#include <opm/material/fluidsystems/SinglePhaseFluidSystem.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/H2ON2LiquidPhaseFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirXyleneFluidSystem.hpp>

#include <opm/material/thermal/FluidThermalConductionLaw.hpp>

// include all fluid states
#include <opm/material/fluidstates/PressureOverlayFluidState.hpp>
#include <opm/material/fluidstates/SaturationOverlayFluidState.hpp>
#include <opm/material/fluidstates/TemperatureOverlayFluidState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/NonEquilibriumFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

// include the tables for CO2 which are delivered with opm-material by default
#include <opm/material/common/UniformTabulated2DFunction.hpp>

namespace Opm {
namespace FluidSystemsTest {
#include <opm/material/components/co2tables.inc>
} }

#include <dune/common/parallel/mpihelper.hh>

// check that the blackoil fluid system implements all non-standard functions
template <class Evaluation, class FluidSystem>
void ensureBlackoilApi()
{
    // here we don't want to call these methods at runtime, we just want to make sure
    // that they compile
    while (false) {
#if HAVE_ECL_INPUT
        Opm::Deck deck;
        Opm::EclipseState eclState(deck);
        FluidSystem::initFromDeck(deck, eclState);
#endif

        typedef typename FluidSystem::Scalar Scalar;
        typedef Opm::CompositionalFluidState<Evaluation, FluidSystem> FluidState;
        FluidState fluidState;
        Evaluation XoG = 0.0;
        Evaluation XgO = 0.0;
        Evaluation Rs = 0.0;
        Evaluation Rv = 0.0;
        Evaluation dummy;

        // some additional typedefs
        typedef typename FluidSystem::OilPvt OilPvt;
        typedef typename FluidSystem::GasPvt GasPvt;
        typedef typename FluidSystem::WaterPvt WaterPvt;

        // check the black-oil specific enums
        static_assert(FluidSystem::numPhases == 3, "");
        static_assert(FluidSystem::numComponents == 3, "");

        static_assert(0 <= FluidSystem::oilPhaseIdx && FluidSystem::oilPhaseIdx < 3, "");
        static_assert(0 <= FluidSystem::gasPhaseIdx && FluidSystem::gasPhaseIdx < 3, "");
        static_assert(0 <= FluidSystem::waterPhaseIdx && FluidSystem::waterPhaseIdx < 3, "");

        static_assert(0 <= FluidSystem::oilCompIdx && FluidSystem::oilCompIdx < 3, "");
        static_assert(0 <= FluidSystem::gasCompIdx && FluidSystem::gasCompIdx < 3, "");
        static_assert(0 <= FluidSystem::waterCompIdx && FluidSystem::waterCompIdx < 3, "");

        // check the non-parser initialization
        std::shared_ptr<OilPvt> oilPvt;
        std::shared_ptr<GasPvt> gasPvt;
        std::shared_ptr<WaterPvt> waterPvt;

        unsigned numPvtRegions = 2;
        FluidSystem::initBegin(numPvtRegions);
        FluidSystem::setEnableDissolvedGas(true);
        FluidSystem::setEnableVaporizedOil(true);
        FluidSystem::setGasPvt(gasPvt);
        FluidSystem::setOilPvt(oilPvt);
        FluidSystem::setWaterPvt(waterPvt);
        FluidSystem::setReferenceDensities(/*oil=*/600.0,
                                           /*water=*/1000.0,
                                           /*gas=*/1.0,
                                           /*regionIdx=*/0);
        FluidSystem::initEnd();

        // the molarMass() method has an optional argument for the PVT region
        unsigned numRegions OPM_UNUSED = FluidSystem::numRegions();
        Scalar Mg OPM_UNUSED = FluidSystem::molarMass(FluidSystem::gasCompIdx,
                                                      /*regionIdx=*/0);
        bool b1 OPM_UNUSED = FluidSystem::enableDissolvedGas();
        bool b2 OPM_UNUSED = FluidSystem::enableVaporizedOil();
        Scalar rhoRefOil OPM_UNUSED = FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx,
                                                                    /*regionIdx=*/0);
        dummy = FluidSystem::convertXoGToRs(XoG, /*regionIdx=*/0);
        dummy = FluidSystem::convertXgOToRv(XgO, /*regionIdx=*/0);
        dummy = FluidSystem::convertXoGToxoG(XoG, /*regionIdx=*/0);
        dummy = FluidSystem::convertXgOToxgO(XgO, /*regionIdx=*/0);
        dummy = FluidSystem::convertRsToXoG(Rs, /*regionIdx=*/0);
        dummy = FluidSystem::convertRvToXgO(Rv, /*regionIdx=*/0);

        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
            dummy = FluidSystem::density(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::saturatedDensity(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::saturatedInverseFormationVolumeFactor(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::viscosity(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::saturatedDissolutionFactor(fluidState, phaseIdx, /*regionIdx=*/0);
            dummy = FluidSystem::saturatedDissolutionFactor(fluidState, phaseIdx, /*regionIdx=*/0, /*maxSo=*/1.0);
            dummy = FluidSystem::saturationPressure(fluidState, phaseIdx, /*regionIdx=*/0);
            for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++ compIdx)
                dummy = FluidSystem::fugacityCoefficient(fluidState, phaseIdx, compIdx,  /*regionIdx=*/0);
        }

        // prevent GCC from producing a "variable assigned but unused" warning
        dummy = 2.0*dummy;


        // the "not considered safe to use directly" API
        const OilPvt& oilPvt2 OPM_UNUSED = FluidSystem::oilPvt();
        const GasPvt& gasPvt2 OPM_UNUSED = FluidSystem::gasPvt();
        const WaterPvt& waterPvt2 OPM_UNUSED = FluidSystem::waterPvt();
    }
}

// check the API of all fluid states
template <class Scalar>
void testAllFluidStates()
{
    typedef Opm::FluidSystems::H2ON2<Scalar> FluidSystem;

    // SimpleModularFluidState
    {   Opm::SimpleModularFluidState<Scalar,
                                     /*numPhases=*/2,
                                     /*numComponents=*/0,
                                     /*FluidSystem=*/void,
                                     /*storePressure=*/false,
                                     /*storeTemperature=*/false,
                                     /*storeComposition=*/false,
                                     /*storeFugacity=*/false,
                                     /*storeSaturation=*/false,
                                     /*storeDensity=*/false,
                                     /*storeViscosity=*/false,
                                     /*storeEnthalpy=*/false> fs;

        checkFluidState<Scalar>(fs); }

    {   Opm::SimpleModularFluidState<Scalar,
                                     /*numPhases=*/2,
                                     /*numComponents=*/2,
                                     FluidSystem,
                                     /*storePressure=*/true,
                                     /*storeTemperature=*/true,
                                     /*storeComposition=*/true,
                                     /*storeFugacity=*/true,
                                     /*storeSaturation=*/true,
                                     /*storeDensity=*/true,
                                     /*storeViscosity=*/true,
                                     /*storeEnthalpy=*/true> fs;

        checkFluidState<Scalar>(fs); }

    // CompositionalFluidState
    {   Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        checkFluidState<Scalar>(fs); }

    // NonEquilibriumFluidState
    {   Opm::NonEquilibriumFluidState<Scalar, FluidSystem> fs;
        checkFluidState<Scalar>(fs); }

    // ImmiscibleFluidState
    {   Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        checkFluidState<Scalar>(fs); }

    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
    BaseFluidState baseFs;

    // TemperatureOverlayFluidState
    {   Opm::TemperatureOverlayFluidState<BaseFluidState> fs(baseFs);
        checkFluidState<Scalar>(fs); }

    // PressureOverlayFluidState
    {   Opm::PressureOverlayFluidState<BaseFluidState> fs(baseFs);
        checkFluidState<Scalar>(fs); }

    // SaturationOverlayFluidState
    {   Opm::SaturationOverlayFluidState<BaseFluidState> fs(baseFs);
        checkFluidState<Scalar>(fs); }
}

template <class Scalar, class FluidStateEval, class LhsEval>
void testAllFluidSystems()
{
    typedef Opm::LiquidPhase<Scalar, Opm::H2O<Scalar>> Liquid;
    typedef Opm::GasPhase<Scalar, Opm::N2<Scalar>> Gas;

    // black-oil
    {
        typedef Opm::FluidSystems::BlackOil<Scalar> FluidSystem;
        if (false) checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>();

        typedef Opm::DenseAd::Evaluation<Scalar, 1> BlackoilDummyEval;
        ensureBlackoilApi<Scalar, FluidSystem>();
        ensureBlackoilApi<BlackoilDummyEval, FluidSystem>();
    }

    // Brine -- CO2
    {   typedef Opm::FluidSystems::BrineCO2<Scalar, Opm::FluidSystemsTest::CO2Tables> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // H2O -- N2
    {   typedef Opm::FluidSystems::H2ON2<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // H2O -- N2 -- liquid phase
    {   typedef Opm::FluidSystems::H2ON2LiquidPhase<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // H2O -- Air
    {   typedef Opm::SimpleH2O<Scalar> H2O;
        typedef Opm::FluidSystems::H2OAir<Scalar, H2O> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Opm::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // H2O -- Air -- Xylene
    {   typedef Opm::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // 2p-immiscible
    {   typedef Opm::FluidSystems::TwoPhaseImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    {   typedef Opm::FluidSystems::TwoPhaseImmiscible<Scalar, Liquid, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    {  typedef Opm::FluidSystems::TwoPhaseImmiscible<Scalar, Gas, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    // 1p
    {   typedef Opm::FluidSystems::SinglePhase<Scalar, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }

    {   typedef Opm::FluidSystems::SinglePhase<Scalar, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem, FluidStateEval, LhsEval>(); }
}

template <class Scalar>
inline void testAll()
{
    typedef Opm::DenseAd::Evaluation<Scalar, 3> Evaluation;

    // ensure that all fluid states are API-compliant
    testAllFluidStates<Scalar>();
    testAllFluidStates<Evaluation>();

    // ensure that all fluid systems are API-compliant: Each fluid system must be usable
    // for both, scalars and function evaluations. The fluid systems for function
    // evaluations must also be usable for scalars.
    testAllFluidSystems<Scalar, /*FluidStateEval=*/Scalar, /*LhsEval=*/Scalar>();
    testAllFluidSystems<Scalar, /*FluidStateEval=*/Evaluation, /*LhsEval=*/Evaluation>();
    testAllFluidSystems<Scalar, /*FluidStateEval=*/Evaluation, /*LhsEval=*/Scalar>();
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
