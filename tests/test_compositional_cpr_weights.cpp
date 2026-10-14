/*
  Copyright 2026 SINTEF Digital

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

#include "config.h"

#define BOOST_TEST_MODULE CompositionalCprWeights
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/material/components/C1.hpp>
#include <opm/material/components/C10.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>
#include <opm/models/ptflash/flashintensivequantities.hh>
#include <opm/models/ptflash/flashlocalresidual.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>

namespace {

constexpr int numComponents = 3;
template<bool waterEnabled>
constexpr int numEq = numComponents + (waterEnabled ? 1 : 0);

template<bool waterEnabled>
using Evaluation = Opm::DenseAd::Evaluation<double, numEq<waterEnabled>>;
template<bool waterEnabled>
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<double, numComponents, waterEnabled>;

template<bool waterEnabled> struct TestContext;
template<bool waterEnabled> struct TestPrimaryVariables;
template<bool waterEnabled> struct TestModel;
struct TestExtensiveQuantities {};

// Stub the surrounding cell/linearizer interfaces. Flash, storage derivatives,
// and true-IMPES weight construction all use their production implementations.
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
    double extrusionFactor() const { return 1.0; }
};

template<bool waterEnabled>
struct TestDiscretizationResidual
{
    Dune::FieldVector<Evaluation<waterEnabled>, numEq<waterEnabled>> residual_{};
    const auto& residual(int) const { return residual_; }
};

struct TestFluxModule
{
    struct FluxIntensiveQuantities {
        template<class Context>
        void update_(const Context&, unsigned, unsigned) {}
    };
};

struct TestGridView { static constexpr int dimensionworld = 3; };

} // namespace

namespace Opm::Properties {
namespace TTag {
template<bool waterEnabled>
struct CompositionalCprTest {};
} // namespace TTag

#define TEST_TYPE_PROPERTY(Name, ...) \
    template<class TypeTag, bool W> \
    struct Name<TypeTag, TTag::CompositionalCprTest<W>> { using type = __VA_ARGS__; }
TEST_TYPE_PROPERTY(Scalar, double);
TEST_TYPE_PROPERTY(Evaluation, ::Evaluation<W>);
TEST_TYPE_PROPERTY(FluidSystem, ::FluidSystem<W>);
TEST_TYPE_PROPERTY(ElementContext, TestContext<W>);
TEST_TYPE_PROPERTY(DiscIntensiveQuantities, TestDiscretization);
TEST_TYPE_PROPERTY(DiscLocalResidual, TestDiscretizationResidual<W>);
TEST_TYPE_PROPERTY(FluxModule, TestFluxModule);
TEST_TYPE_PROPERTY(MaterialLaw, TestMaterialLaw);
TEST_TYPE_PROPERTY(MaterialLawParams, TestMaterialLaw::Params);
TEST_TYPE_PROPERTY(GridView, TestGridView);
TEST_TYPE_PROPERTY(ThreadManager, void);
TEST_TYPE_PROPERTY(FlashSolver, Opm::PTFlash<double, ::FluidSystem<W>>);
TEST_TYPE_PROPERTY(Indices, Opm::FlashIndices<TypeTag, 0>);
TEST_TYPE_PROPERTY(IntensiveQuantities, Opm::FlashIntensiveQuantities<TypeTag>);
TEST_TYPE_PROPERTY(EqVector, Dune::FieldVector<double, numEq<W>>);
TEST_TYPE_PROPERTY(RateVector, Dune::FieldVector<::Evaluation<W>, numEq<W>>);
TEST_TYPE_PROPERTY(PrimaryVariables, TestPrimaryVariables<W>);
TEST_TYPE_PROPERTY(ExtensiveQuantities, TestExtensiveQuantities);
TEST_TYPE_PROPERTY(Model, TestModel<W>);
#undef TEST_TYPE_PROPERTY

#define TEST_VALUE_PROPERTY(Name, Value) \
    template<class TypeTag, bool W> \
    struct Name<TypeTag, TTag::CompositionalCprTest<W>> { static constexpr auto value = Value; }
TEST_VALUE_PROPERTY(NumComponents, numComponents);
TEST_VALUE_PROPERTY(NumEq, numEq<W>);
TEST_VALUE_PROPERTY(NumPhases, W ? 3 : 2);
TEST_VALUE_PROPERTY(EnableEnergy, false);
TEST_VALUE_PROPERTY(EnableDiffusion, false);
TEST_VALUE_PROPERTY(EnableWater, W);
#undef TEST_VALUE_PROPERTY
} // namespace Opm::Properties

namespace {

template<bool waterEnabled>
using TypeTag = Opm::Properties::TTag::CompositionalCprTest<waterEnabled>;
template<bool waterEnabled>
using IntensiveQuantities = Opm::FlashIntensiveQuantities<TypeTag<waterEnabled>>;
template<bool waterEnabled>
using LocalResidual = Opm::FlashLocalResidual<TypeTag<waterEnabled>>;

template<bool waterEnabled>
struct TestPrimaryVariables
{
    std::array<double, numEq<waterEnabled>> values{};
    Evaluation<waterEnabled> makeEvaluation(unsigned index, unsigned) const
    { return Evaluation<waterEnabled>{values[index], static_cast<int>(index)}; }
};

struct TestProblem
{
    double temperature_ = 300.0;
    template<class Context>
    double temperature(const Context&, unsigned, unsigned) const { return temperature_; }
    Opm::CompositionalConfig::EOSType getEosType() const
    { return Opm::CompositionalConfig::EOSType::PR; }
    template<class Context>
    TestMaterialLaw::Params materialLawParams(const Context&, unsigned, unsigned) const { return {}; }
    template<class Context>
    auto porosity(const Context& context, unsigned, unsigned) const
    {
        const auto pressure = context.primaryVars(0, 0).makeEvaluation(0, 0);
        return 0.2 * (1.0 + 1.0e-10 * (pressure - 1.0e7));
    }
    template<class Context>
    Dune::FieldMatrix<double, 3, 3> intrinsicPermeability(const Context&, unsigned, unsigned) const
    { return Dune::FieldMatrix<double, 3, 3>{1.0}; }
};

struct TestGrid
{
    Opm::Parallel::Communication comm() const
    { return Opm::Parallel::Communication{Dune::MPIHelper::getCommunicator()}; }
};

struct TestVanguard
{
    TestGrid grid_;
    const TestGrid& grid() const { return grid_; }
    int cartesianIndex(unsigned index) const { return index; }
};

template<bool waterEnabled>
struct TestSimulator
{
    TestPrimaryVariables<waterEnabled> primary;
    TestProblem problem;
    TestVanguard vanguard_;
    double timeStepSize() const { return 86400.0; }
    const TestVanguard& vanguard() const { return vanguard_; }
};

template<bool waterEnabled>
struct TestContext
{
    const TestSimulator<waterEnabled>& simulator_;
    IntensiveQuantities<waterEnabled> iq;
    explicit TestContext(const TestSimulator<waterEnabled>& simulator) : simulator_(simulator) {}
    const auto& simulator() const { return simulator_; }
    void updatePrimaryStencil(int) {}
    void updatePrimaryIntensiveQuantities(unsigned timeIdx) { iq.update(*this, 0, timeIdx); }
    const auto& primaryVars(unsigned, unsigned) const { return simulator_.primary; }
    const TestProblem& problem() const { return simulator_.problem; }
    const auto& intensiveQuantities(unsigned, unsigned) const { return iq; }
    const IntensiveQuantities<waterEnabled>* thermodynamicHint(unsigned, unsigned) const
    { return nullptr; }
    unsigned globalSpaceIndex(unsigned, unsigned) const { return 0; }
    struct SubControlVolume { double volume() const { return 100.0; } };
    struct Stencil {
        SubControlVolume subControlVolume(unsigned) const { return {}; }
    };
    Stencil stencil(unsigned) const { return {}; }
};

template<bool waterEnabled>
struct TestLocalLinearizer
{
    LocalResidual<waterEnabled> residual_;
    const auto& localResidual() const { return residual_; }
};

template<bool waterEnabled>
struct TestLinearizer
{
    struct MatrixAdapter {
        using MatrixBlock = Opm::MatrixBlock<double, numEq<waterEnabled>, numEq<waterEnabled>>;
    };
    MatrixAdapter jacobian() const { return {}; }
};

template<bool waterEnabled>
struct TestModel
{
    TestLocalLinearizer<waterEnabled> local;
    TestLinearizer<waterEnabled> linearizer_;
    const auto& localLinearizer(std::size_t) const { return local; }
    const auto& linearizer() const { return linearizer_; }
};

template<bool waterEnabled>
void initFluidSystem()
{
    using FS = FluidSystem<waterEnabled>;
    FS::init();
    const auto addComponent = []<class Component>() {
        FS::addComponent(typename FS::ComponentParam{
            Component::name(), Component::molarMass(), Component::criticalTemperature(),
            Component::criticalPressure(), Component::criticalVolume(),
            Component::acentricFactor(), 0.0});
    };
    addComponent.template operator()<Opm::SimpleCO2<double>>();
    addComponent.template operator()<Opm::C1<double>>();
    addComponent.template operator()<Opm::C10<double>>();

    if constexpr (waterEnabled) {
        auto water = std::make_shared<typename FS::WaterPvt>();
        constexpr auto approach = Opm::WaterPvtApproach::ConstantCompressibilityWater;
        water->setApproach(approach);
        auto& pvt = water->template getRealPvt<approach>();
        pvt.setNumRegions(1);
        pvt.setReferenceDensities(0, 800.0, 1.0, 1000.0);
        pvt.setReferencePressure(0, 1.0e5);
        pvt.setCompressibility(0, 4.5e-10);
        pvt.setViscosity(0, 1.0e-3);
        pvt.setReferenceFormationVolumeFactor(0, 1.0);
        FS::setWaterPvt(water);
    }
}

enum class HydrocarbonPhase { Gas, Oil, TwoPhase };

template<bool waterEnabled>
void checkWeights(double temperature,
                  const std::array<double, numEq<waterEnabled>>& primaryVariables,
                  HydrocarbonPhase expectedPhase)
{
    using Indices = Opm::FlashIndices<TypeTag<waterEnabled>, 0>;
    constexpr int pressureIdx = Indices::pressureSwitchIdx;
    TestSimulator<waterEnabled> simulator;
    simulator.problem.temperature_ = temperature;
    simulator.primary.values = primaryVariables;
    TestContext<waterEnabled> context{simulator};
    TestModel<waterEnabled> model;
    const std::array<std::array<int, 1>, 1> chunks{{{{0}}}};
    Dune::BlockVector<Dune::FieldVector<double, numEq<waterEnabled>>> weights(1);

    // Check the phase split explicitly so each case exercises its intended EOS
    // branch. This also keeps any flash failure outside the OpenMP weight loop.
    context.updatePrimaryIntensiveQuantities(0);
    const double liquidFraction = context.iq.fluidState().L().value();
    switch (expectedPhase) {
    case HydrocarbonPhase::Gas:
        BOOST_REQUIRE_SMALL(liquidFraction, 1.0e-10);
        break;
    case HydrocarbonPhase::Oil:
        BOOST_REQUIRE_SMALL(1.0 - liquidFraction, 1.0e-10);
        break;
    case HydrocarbonPhase::TwoPhase:
        BOOST_REQUIRE_GT(liquidFraction, 1.0e-10);
        BOOST_REQUIRE_LT(liquidFraction, 1.0 - 1.0e-10);
        break;
    }
    if constexpr (waterEnabled) {
        BOOST_CHECK_EQUAL(context.iq.hasHydrocarbon(), primaryVariables[Indices::water0Idx] < 1.0);
    }

    Opm::Amg::getTrueImpesWeights(pressureIdx, weights, context, model, chunks, false);
    double maxAbsWeight = 0.0;
    for (const double weight : weights[0]) {
        BOOST_REQUIRE(std::isfinite(weight));
        maxAbsWeight = std::max(maxAbsWeight, std::abs(weight));
    }
    BOOST_CHECK_CLOSE(maxAbsWeight, 1.0, 1.0e-10);

    Dune::FieldVector<Evaluation<waterEnabled>, numEq<waterEnabled>> storage;
    model.local.localResidual().computeStorage(storage, context, 0, 0);
    // A CPR pressure equation must retain pressure dependence while eliminating
    // the storage sensitivities to the independent mole fractions (and Sw).
    // Scale cancellation by the contributing terms to cover both normal and
    // nearly water-only cells without a unit-dependent absolute tolerance.
    for (int variable = 0; variable < numEq<waterEnabled>; ++variable) {
        BOOST_TEST_CONTEXT("water enabled = " << waterEnabled << ", primary variable = " << variable) {
            double weightedDerivative = 0.0;
            double scale = 0.0;
            for (int equation = 0; equation < numEq<waterEnabled>; ++equation) {
                const double term = weights[0][equation] * storage[equation].derivative(variable);
                BOOST_REQUIRE(std::isfinite(term));
                weightedDerivative += term;
                scale += std::abs(term);
            }
            BOOST_REQUIRE_GT(scale, 0.0);
            if (variable == pressureIdx) {
                BOOST_CHECK_GT(weightedDerivative, 0.0);
            } else {
                BOOST_CHECK_SMALL(weightedDerivative / scale, 1.0e-7);
            }
        }
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(GasStorage)
{
    checkWeights<false>(450.0, {1.0e6, 0.001, 0.998}, HydrocarbonPhase::Gas);
    checkWeights<true>(450.0, {1.0e6, 0.001, 0.998, 0.3}, HydrocarbonPhase::Gas);
}

BOOST_AUTO_TEST_CASE(OilStorage)
{
    checkWeights<false>(300.0, {1.0e7, 0.001, 0.001}, HydrocarbonPhase::Oil);
    checkWeights<true>(300.0, {1.0e7, 0.001, 0.001, 0.3}, HydrocarbonPhase::Oil);
}

BOOST_AUTO_TEST_CASE(TwoPhaseStorage)
{
    checkWeights<false>(300.0, {1.0e7, 0.5, 0.3}, HydrocarbonPhase::TwoPhase);
    checkWeights<true>(300.0, {1.0e7, 0.5, 0.3, 0.0}, HydrocarbonPhase::TwoPhase);
    checkWeights<true>(300.0, {1.0e7, 0.5, 0.3, 0.3}, HydrocarbonPhase::TwoPhase);
}

BOOST_AUTO_TEST_CASE(VanishingComponents)
{
    // Pure CO2 activates the composition floor for both other components.
    checkWeights<false>(300.0, {1.0e7, 1.0, 0.0}, HydrocarbonPhase::Oil);
    checkWeights<true>(300.0, {1.0e7, 1.0, 0.0, 0.3}, HydrocarbonPhase::Oil);
}

BOOST_AUTO_TEST_CASE(NearlyWaterOnlyStorage)
{
    checkWeights<true>(300.0, {1.0e7, 0.5, 0.3, 1.0 - 1.0e-9}, HydrocarbonPhase::TwoPhase);
}

BOOST_AUTO_TEST_CASE(WaterOnlyStorage)
{
    // The hydrocarbon floor must retain derivatives even when the physical cell
    // contains only water, otherwise true-IMPES weight construction is singular.
    checkWeights<true>(300.0, {1.0e7, 0.5, 0.3, 1.0}, HydrocarbonPhase::TwoPhase);
}

bool init_unit_test_func()
{
    Opm::Parameters::Register<Opm::Parameters::FlashTolerance<double>>("Flash tolerance");
    Opm::Parameters::Register<Opm::Parameters::FlashVerbosity>("Flash verbosity");
    Opm::Parameters::Register<Opm::Parameters::FlashTwoPhaseMethod>("Flash method");
    Opm::Parameters::SetDefault<Opm::Parameters::FlashTolerance<double>>(1.0e-8);
    Opm::Parameters::endRegistration();
    initFluidSystem<false>();
    initFluidSystem<true>();
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
