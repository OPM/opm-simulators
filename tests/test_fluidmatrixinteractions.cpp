/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \brief This test makes sure that the API for fluid-matrix
 *        interactions is observed by all capillary pressure / relperm
 *        laws.
 */
#include "config.h"

// include all capillary pressure laws
#include <opm/material/fluidmatrixinteractions/BrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/ParkerLenhard.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/VanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/SplineTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/ThreePhaseParkerVanGenuchten.hpp>

// include the helper classes to construct traits
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

// include some fluid states
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

// include some fluid systems
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

// include some components
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/N2.hpp>

#include <opm/core/utility/Unused.hpp>

// include the MPI header if available
#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

// class to call MPI_Init() on construction and MPI_Finalize() in the
// destructor
class MyMpiHelper
{
public:
    MyMpiHelper(int &argc, char **&argv)
    {
#if HAVE_MPI
        MPI_Init(&argc, &argv);
#endif // HAVE_MPI
    };

    ~MyMpiHelper()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif // HAVE_MPI
    };
};

// this function makes sure that a capillary pressure law adheres to
// the generic programming interface for such laws. This API _must_ be
// implemented by all capillary pressure laws. If there are no _very_
// strong reasons to do otherwise, numerical models should only use on
// this API.
template <class MaterialLaw, class FluidState>
void testGenericApi()
{
    while (0) {
        // ensure the presence of the required values
        static const int numPhases = MaterialLaw::numPhases;

        // check for the presence of the is*Dependent values
        OPM_UNUSED static const bool isSaturationDependent = MaterialLaw::isSaturationDependent;
        OPM_UNUSED static const bool isPressureDependent = MaterialLaw::isPressureDependent;
        OPM_UNUSED static const bool isTemperatureDependent = MaterialLaw::isTemperatureDependent;
        OPM_UNUSED static const bool isCompositionDependent = MaterialLaw::isCompositionDependent;

        // Make sure that the Traits, Params and Scalar typedefs are
        // exported by the material law
        typedef typename MaterialLaw::Params Params;
        typedef typename MaterialLaw::Traits Traits;
        typedef typename MaterialLaw::Scalar Scalar;
        typedef typename MaterialLaw::Traits::Scalar TraitsScalar;

        static_assert(std::is_same<Scalar, TraitsScalar>::value,
                      "The traits and the material law must use the same type as scalar value");
        static_assert(numPhases == Traits::numPhases,
                      "The traits and the material law must use the number of fluid phases");

        // check the API of the parameter class. setting the actual
        // parameter values is implementation specific. But all
        // parameters must be default and copy constructible as well
        // as exhibit the finalize() method!
        Params params;
        params.finalize();
        const Params paramsConst(params);

        // test the generic methods which need to be implemented by
        // all material laws
        const FluidState fs;
        double destValues[numPhases];
        MaterialLaw::capillaryPressures(destValues, paramsConst, fs);
        MaterialLaw::saturations(destValues, paramsConst, fs);
        MaterialLaw::relativePermeabilities(destValues, paramsConst, fs);

        std::array<Scalar, numPhases> dpc;
        MaterialLaw::dCapillaryPressures_dSaturation(dpc, paramsConst, fs, /*phaseIdx=*/0);
        MaterialLaw::dCapillaryPressures_dPressure(dpc, paramsConst, fs, /*phaseIdx=*/0);
        MaterialLaw::dCapillaryPressures_dTemperature(dpc, paramsConst, fs);
        MaterialLaw::dCapillaryPressures_dMoleFraction(dpc, paramsConst, fs, /*phaseIdx=*/0, /*compIdx=*/0);

        std::array<Scalar, numPhases> dkr;
        MaterialLaw::dRelativePermeabilities_dSaturation(dkr, paramsConst, fs, /*phaseIdx=*/0);
        MaterialLaw::dRelativePermeabilities_dPressure(dkr, paramsConst, fs, /*phaseIdx=*/0);
        MaterialLaw::dRelativePermeabilities_dTemperature(dkr, paramsConst, fs);
        MaterialLaw::dRelativePermeabilities_dMoleFraction(dkr, paramsConst, fs, /*phaseIdx=*/0, /*compIdx=*/0);
    }
}

// this function makes ensures that a pressure law adheres to the
// covenience programming interface for two-phase material laws. The
// main purpose of this interface is to simplify the implementation of
// nested material laws.
template <class MaterialLaw, class FluidState>
void testTwoPhaseApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static const int numPhases = MaterialLaw::numPhases;
        static_assert(numPhases == 2,
                      "The number of fluid phases for a twophase "
                      "capillary pressure law must be 2");
        static_assert(MaterialLaw::implementsTwoPhaseApi,
                      "This material law is expected to implement "
                      "the two-phase API!");

        OPM_UNUSED static const int wettingPhaseIdx = MaterialLaw::wettingPhaseIdx;
        OPM_UNUSED static const int nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx;

        // make sure the two-phase specific methods are present
        const FluidState fs;
        const typename MaterialLaw::Params params;

        OPM_UNUSED Scalar v;
        v = MaterialLaw::pcnw(params, fs);
        v = MaterialLaw::Sw(params, fs);
        v = MaterialLaw::Sn(params, fs);
        v = MaterialLaw::krw(params, fs);
        v = MaterialLaw::krn(params, fs);
    }
}

template <class MaterialLaw, class FluidState>
void testTwoPhaseSatApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static_assert(MaterialLaw::implementsTwoPhaseSatApi,
                      "This material law is expected to implement "
                      "the two-phase saturation only API!");
        static_assert(!MaterialLaw::isPressureDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on the absolute phase pressures!");
        static_assert(!MaterialLaw::isTemperatureDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on temperature!");
        static_assert(!MaterialLaw::isCompositionDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on the phase compositions!");

        OPM_UNUSED static const int numPhases = MaterialLaw::numPhases;

        // make sure the two-phase specific methods are present
        const typename MaterialLaw::Params params;

        Scalar Sw = 0;
        OPM_UNUSED Scalar v;
        v = MaterialLaw::twoPhaseSatPcnw(params, Sw);
        v = MaterialLaw::twoPhaseSatSw(params, Sw);
        v = MaterialLaw::twoPhaseSatSn(params, Sw);
        v = MaterialLaw::twoPhaseSatKrw(params, Sw);
        v = MaterialLaw::twoPhaseSatKrn(params, Sw);
    }
}

template <class MaterialLaw, class FluidState>
void testThreePhaseApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static const int numPhases = MaterialLaw::numPhases;
        static_assert(numPhases == 3,
                      "The number of fluid phases for a threephase "
                      "capillary pressure law must be 3");

        OPM_UNUSED static const int wettingPhaseIdx = MaterialLaw::wettingPhaseIdx;
        OPM_UNUSED static const int nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx;
        OPM_UNUSED static const int gasPhaseIdx = MaterialLaw::gasPhaseIdx;

        // make sure the two-phase specific methods are present
        const FluidState fs;
        const typename MaterialLaw::Params params;

        OPM_UNUSED Scalar v;
        v = MaterialLaw::pcnw(params, fs);
        v = MaterialLaw::Sw(params, fs);
        v = MaterialLaw::Sn(params, fs);
        v = MaterialLaw::Sg(params, fs);
        v = MaterialLaw::krw(params, fs);
        v = MaterialLaw::krn(params, fs);
        v = MaterialLaw::krg(params, fs);
    }
}

template <class MaterialLaw>
void testThreePhaseSatApi()
{
}

int main(int argc, char **argv)
{
    typedef double Scalar;
    typedef Opm::SimpleH2O<Scalar> H2O;
    typedef Opm::N2<Scalar> N2;

    typedef Opm::LiquidPhase<Scalar, H2O> Liquid;
    typedef Opm::GasPhase<Scalar, N2> Gas;

    typedef Opm::FluidSystems::TwoPhaseImmiscible<Scalar, Liquid, Gas> TwoPFluidSystem;
    typedef Opm::FluidSystems::BlackOil<Scalar> ThreePFluidSystem;

    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        TwoPFluidSystem::wettingPhaseIdx,
                                        TwoPFluidSystem::nonWettingPhaseIdx> TwoPhaseTraits;

    typedef Opm::ThreePhaseMaterialTraits<Scalar,
                                          ThreePFluidSystem::waterPhaseIdx,
                                          ThreePFluidSystem::oilPhaseIdx,
                                          ThreePFluidSystem::gasPhaseIdx> ThreePhaseTraits;

    typedef Opm::ImmiscibleFluidState<Scalar, TwoPFluidSystem> TwoPhaseFluidState;
    typedef Opm::ImmiscibleFluidState<Scalar, ThreePFluidSystem> ThreePhaseFluidState;

    MyMpiHelper mpiHelper(argc, argv);

    // test conformance to the capillary pressure APIs
    {
        typedef Opm::BrooksCorey<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::LinearMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();

        typedef Opm::EffToAbsLaw<MaterialLaw> TwoPAbsLaw;
        testGenericApi<TwoPAbsLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<TwoPAbsLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<TwoPAbsLaw, TwoPhaseFluidState>();

        typedef Opm::LinearMaterial<ThreePhaseTraits> ThreePMaterialLaw;
        testGenericApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePMaterialLaw, ThreePhaseFluidState>();

        typedef Opm::EffToAbsLaw<ThreePMaterialLaw> ThreePAbsLaw;
        testGenericApi<ThreePAbsLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePAbsLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePAbsLaw, ThreePhaseFluidState>();
    }
    {
        typedef Opm::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Opm::EclDefaultMaterial<ThreePhaseTraits,
                                        /*GasOilMaterial=*/TwoPhaseMaterial,
                                        /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Opm::ThreePhaseParkerVanGenuchten<ThreePhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Opm::NullMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();

        typedef Opm::NullMaterial<ThreePhaseTraits> ThreePMaterialLaw;
        testGenericApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePMaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Opm::ParkerLenhard<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::SplineTwoPhaseMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::VanGenuchten<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::RegularizedBrooksCorey<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Opm::RegularizedVanGenuchten<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }

    return 0;
}
