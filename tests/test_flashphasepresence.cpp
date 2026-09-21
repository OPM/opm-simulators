/*
  Copyright 2026 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#define BOOST_TEST_MODULE FlashPhasePresence
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/ptflash/flashintensivequantities.hh>

#include <array>

namespace {

using Evaluation = Opm::DenseAd::Evaluation<double, 3>;

// Fixed phase split and molar volumes isolate saturation bookkeeping from the
// EOS. The production intensive-quantity update still computes its AD values,
// applies the hydrocarbon floor, and provides the reporting accessors.
template<bool waterEnabled>
struct TestFluidSystem
{
    using Scalar = double;
    static constexpr int numComponents = 2;
    static constexpr int numPhases = waterEnabled ? 3 : 2;
    static constexpr unsigned oilPhaseIdx = 0;
    static constexpr unsigned gasPhaseIdx = 1;
    static constexpr unsigned waterPhaseIdx = 2;

    static bool phaseIsActive(unsigned phase) { return phase < numPhases; }
    static double molarMass(unsigned) { return 0.02; }
    static double acentricFactor(unsigned) { return 0.0; }
    static double criticalTemperature(unsigned) { return 300.0; }
    static double criticalPressure(unsigned) { return 1.0e6; }

    template<class Eval>
    struct ParameterCache {
        explicit ParameterCache(int) {}
        template<class FluidState>
        void updatePhase(const FluidState&, unsigned) {}
        Eval molarVolume(unsigned) const { return Eval{1.0}; }
        Eval correctedMolarVolume(unsigned) const { return Eval{1.0}; }
    };

    template<class FluidState, class Cache>
    static auto density(const FluidState&, const Cache&, unsigned)
    { return typename FluidState::ValueType{1000.0}; }

    template<class FluidState, class Cache>
    static auto viscosity(const FluidState&, const Cache&, unsigned)
    { return typename FluidState::ValueType{0.001}; }
};

struct FixedFlash
{
    static inline double liquidFraction = 0.4;

    template<class FluidState, class Method>
    static void solve(FluidState& fs, Method, double, int, int)
    {
        fs.setLvalue(Evaluation{liquidFraction});
        for (unsigned phase = 0; phase < 2; ++phase) {
            for (unsigned comp = 0; comp < 2; ++comp) {
                fs.setMoleFraction(phase, comp, fs.moleFraction(comp));
            }
        }
    }
};

struct TestMaterialLaw
{
    struct Params {};
    template<class Values, class FluidState>
    static void relativePermeabilities(Values& values, const Params&, const FluidState& fs)
    {
        for (unsigned phase = 0; phase < values.size(); ++phase) {
            values[phase] = fs.saturation(phase);
        }
    }
};

struct TestDiscretization
{
    template<class Context>
    void update(const Context&, unsigned, unsigned) {}
};

struct TestFluxModule
{
    struct FluxIntensiveQuantities {
        template<class Context>
        void update_(const Context&, unsigned, unsigned) {}
    };
};

struct TestGridView { static constexpr int dimensionworld = 3; };

template<bool waterEnabled>
struct TestContext;

} // namespace

namespace Opm::Properties {
namespace TTag {
template<bool waterEnabled>
struct FlashPresenceTest {};
} // namespace TTag

#define TEST_TYPE_PROPERTY(Name, ...) \
    template<class TypeTag, bool W> \
    struct Name<TypeTag, TTag::FlashPresenceTest<W>> { using type = __VA_ARGS__; }
TEST_TYPE_PROPERTY(Scalar, double);
TEST_TYPE_PROPERTY(Evaluation, ::Evaluation);
TEST_TYPE_PROPERTY(FluidSystem, TestFluidSystem<W>);
TEST_TYPE_PROPERTY(ElementContext, TestContext<W>);
TEST_TYPE_PROPERTY(DiscIntensiveQuantities, TestDiscretization);
TEST_TYPE_PROPERTY(FluxModule, TestFluxModule);
TEST_TYPE_PROPERTY(MaterialLaw, TestMaterialLaw);
TEST_TYPE_PROPERTY(MaterialLawParams, TestMaterialLaw::Params);
TEST_TYPE_PROPERTY(GridView, TestGridView);
TEST_TYPE_PROPERTY(ThreadManager, void);
TEST_TYPE_PROPERTY(FlashSolver, FixedFlash);
TEST_TYPE_PROPERTY(Indices, Opm::FlashIndices<TypeTag, 0>);
#undef TEST_TYPE_PROPERTY

#define TEST_VALUE_PROPERTY(Name, Value) \
    template<class TypeTag, bool W> \
    struct Name<TypeTag, TTag::FlashPresenceTest<W>> { static constexpr auto value = Value; }
TEST_VALUE_PROPERTY(NumComponents, 2);
TEST_VALUE_PROPERTY(NumPhases, W ? 3 : 2);
TEST_VALUE_PROPERTY(EnableEnergy, false);
TEST_VALUE_PROPERTY(EnableDiffusion, false);
TEST_VALUE_PROPERTY(EnableWater, W);
#undef TEST_VALUE_PROPERTY
} // namespace Opm::Properties

namespace {

template<bool waterEnabled>
using IntensiveQuantities = Opm::FlashIntensiveQuantities<
    Opm::Properties::TTag::FlashPresenceTest<waterEnabled>>;

struct TestPrimaryVariables
{
    std::array<double, 3> values{1.0e7, 0.3, 1.0};
    Evaluation makeEvaluation(unsigned index, unsigned) const
    { return Evaluation{values[index], static_cast<int>(index)}; }
};

struct TestProblem
{
    template<class Context>
    double temperature(const Context&, unsigned, unsigned) const { return 300.0; }
    int getEosType() const { return 0; }
    template<class Context>
    TestMaterialLaw::Params materialLawParams(const Context&, unsigned, unsigned) const { return {}; }
    template<class Context>
    double porosity(const Context&, unsigned, unsigned) const { return 0.2; }
    template<class Context>
    Dune::FieldMatrix<double, 3, 3> intrinsicPermeability(const Context&, unsigned, unsigned) const
    { return Dune::FieldMatrix<double, 3, 3>{1.0}; }
};

template<bool waterEnabled>
struct TestContext
{
    TestPrimaryVariables priVars;
    TestProblem problem_;
    const TestPrimaryVariables& primaryVars(unsigned, unsigned) const { return priVars; }
    const TestProblem& problem() const { return problem_; }
    const IntensiveQuantities<waterEnabled>* thermodynamicHint(unsigned, unsigned) const
    { return nullptr; }
    unsigned globalSpaceIndex(unsigned, unsigned) const { return 0; }
};

struct ParametersFixture
{
    ParametersFixture()
    {
        Opm::Parameters::Register<Opm::Parameters::FlashTolerance<double>>("Flash tolerance");
        Opm::Parameters::Register<Opm::Parameters::FlashVerbosity>("Flash verbosity");
        Opm::Parameters::Register<Opm::Parameters::FlashTwoPhaseMethod>("Flash method");
        Opm::Parameters::endRegistration();
    }
};

} // namespace

BOOST_GLOBAL_FIXTURE(ParametersFixture);

BOOST_AUTO_TEST_CASE(PureWaterKeepsSolverDerivativesButReportsNoHydrocarbon)
{
    TestContext<true> context;
    IntensiveQuantities<true> iq;
    FixedFlash::liquidFraction = 0.4;
    iq.update(context, 0, 0);

    BOOST_CHECK(!iq.hasHydrocarbon());
    for (unsigned phase = 0; phase < 2; ++phase) {
        BOOST_CHECK_GT(iq.fluidState().saturation(phase).value(), 0.0);
        BOOST_CHECK_LT(iq.fluidState().saturation(phase).derivative(2), 0.0);
        BOOST_CHECK(!iq.phaseIsPresent(phase));
        BOOST_CHECK_EQUAL(iq.saturationForOutput(phase), 0.0);
    }
    BOOST_CHECK(iq.phaseIsPresent(2));
    BOOST_CHECK_EQUAL(iq.saturationForOutput(2), 1.0);

    // A component's storage retains composition sensitivity in the solver.
    const auto storage = iq.fluidState().saturation(0) * iq.fluidState().moleFraction(0, 0);
    BOOST_CHECK_GT(storage.derivative(1), 0.0);
}

BOOST_AUTO_TEST_CASE(PositiveHydrocarbonBelowTheFloorIsStillPresent)
{
    TestContext<true> context;
    context.priVars.values[2] = 1.0 - 1.0e-10;
    IntensiveQuantities<true> iq;
    FixedFlash::liquidFraction = 0.4;
    iq.update(context, 0, 0);

    BOOST_CHECK(iq.hasHydrocarbon());
    for (unsigned phase = 0; phase < 3; ++phase) {
        BOOST_CHECK(iq.phaseIsPresent(phase));
        BOOST_CHECK_EQUAL(iq.saturationForOutput(phase), iq.fluidState().saturation(phase).value());
    }
}

BOOST_AUTO_TEST_CASE(WaterFreeTraceGasIsNotSuppressed)
{
    TestContext<false> context;
    IntensiveQuantities<false> iq;
    FixedFlash::liquidFraction = 1.0 - 5.0e-8;
    iq.update(context, 0, 0);

    BOOST_CHECK(iq.hasHydrocarbon());
    BOOST_CHECK(iq.phaseIsPresent(0));
    BOOST_CHECK(iq.phaseIsPresent(1));
    BOOST_CHECK_CLOSE(iq.saturationForOutput(1), 5.0e-8, 1e-6);
    // Querying the inactive water phase must not index its missing phase slot.
    BOOST_CHECK(!iq.phaseIsPresent(2));
    BOOST_CHECK_EQUAL(iq.saturationForOutput(2), 0.0);
}

BOOST_AUTO_TEST_CASE(PresenceIsRefreshedWhenHydrocarbonEntersAndLeaves)
{
    TestContext<true> context;
    IntensiveQuantities<true> iq;
    FixedFlash::liquidFraction = 0.0; // All gas when hydrocarbon is present.
    iq.update(context, 0, 0);
    BOOST_CHECK(!iq.hasHydrocarbon());

    context.priVars.values[2] = 0.8;
    iq.update(context, 0, 0);
    BOOST_CHECK(iq.hasHydrocarbon());
    BOOST_CHECK(!iq.phaseIsPresent(0));
    BOOST_CHECK(iq.phaseIsPresent(1));
    BOOST_CHECK_CLOSE(iq.saturationForOutput(1), 0.2, 1e-10);
    const auto copied = iq;
    BOOST_CHECK(copied.phaseIsPresent(1));

    context.priVars.values[2] = 1.0;
    iq.update(context, 0, 0);
    BOOST_CHECK(!iq.hasHydrocarbon());
    BOOST_CHECK(!iq.phaseIsPresent(1));
    BOOST_CHECK_EQUAL(iq.saturationForOutput(2), 1.0);
}
