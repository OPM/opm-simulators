#include <config.h>

#include "TestTypeTag.hpp"

#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#define BOOST_TEST_MODULE BlackoilPrimaryVariables
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <tuple>

template<class TypeTag>
class StubBlackoilPrimaryVariablesProblem
{
public:
    using MaterialLawParams = Opm::GetPropType<TypeTag, Opm::Properties::MaterialLawParams>;

    unsigned satnumRegionIndex(unsigned) const
    {
        return 0;
    }

    double temperature(unsigned, unsigned) const
    {
        return 300.0;
    }

    double referencePorosity(unsigned, unsigned) const
    {
        return 0.0;
    }

    const MaterialLawParams& materialLawParams(unsigned) const
    {
        throw std::logic_error("materialLawParams() should not be reached in this test");
    }

    double maxGasDissolutionFactor(unsigned, unsigned) const
    {
        return 1.0;
    }

    double maxOilSaturation(unsigned) const
    {
        return 1.0;
    }

    double maxOilVaporizationFactor(unsigned, unsigned) const
    {
        return 1.0;
    }
};

namespace Opm::Properties {

namespace TTag {

struct ClassicalBlackoilPrimaryVariablesTestTypeTag {
    using InheritsFrom = std::tuple<TestTypeTag>;
};

} // namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::ClassicalBlackoilPrimaryVariablesTestTypeTag>
{
    using type = ::StubBlackoilPrimaryVariablesProblem<TypeTag>;
};

} // namespace Opm::Properties

namespace {

using TypeTag = Opm::Properties::TTag::ClassicalBlackoilPrimaryVariablesTestTypeTag;
using PrimaryVariables = Opm::BlackOilPrimaryVariables<TypeTag>;
using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
using Problem = Opm::GetPropType<TypeTag, Opm::Properties::Problem>;

static_assert(Opm::getPropValue<TypeTag, Opm::Properties::EnergyModuleType>()
              == Opm::EnergyModules::NoTemperature);

struct FluidSystemGuard
{
    FluidSystemGuard()
        : dissolvedGasInWater(FluidSystem::enableDissolvedGasInWater())
        , reservoirTemperature(FluidSystem::reservoirTemperature())
    {
    }

    ~FluidSystemGuard()
    {
        FluidSystem::setReservoirTemperature(reservoirTemperature);
        FluidSystem::setEnableDissolvedGasInWater(dissolvedGasInWater);
    }

    bool dissolvedGasInWater;
    double reservoirTemperature;
};

PrimaryVariables makePrimaryVariables(PrimaryVariables::GasMeaning gasMeaning,
                                      double waterSwitch,
                                      double compositionSwitch)
{
    PrimaryVariables priVars;
    for (std::size_t i = 0; i < priVars.size(); ++i) {
        priVars[i] = 0.0;
    }

    priVars.setPrimaryVarsMeaningWater(PrimaryVariables::WaterMeaning::Sw);
    priVars.setPrimaryVarsMeaningGas(gasMeaning);
    priVars.setPrimaryVarsMeaningPressure(PrimaryVariables::PressureMeaning::Po);
    priVars.setPrimaryVarsMeaningBrine(PrimaryVariables::BrineMeaning::Disabled);
    priVars.setPrimaryVarsMeaningSolvent(PrimaryVariables::SolventMeaning::Disabled);
    priVars[Indices::waterSwitchIdx] = waterSwitch;
    priVars[Indices::compositionSwitchIdx] = compositionSwitch;
    return priVars;
}

} // namespace

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchesToGasSaturationInWaterFilledCells)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setReservoirTemperature(300.0);

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Rs, 1.2, 0.37);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/0.95);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK(priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw);
    BOOST_CHECK_CLOSE(priVars[Indices::waterSwitchIdx], 1.0, 1e-12);
    BOOST_CHECK_SMALL(priVars[Indices::compositionSwitchIdx], 1e-12);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesKeepsGasMeaningWhenAlreadySg)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setReservoirTemperature(300.0);

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg, 1.1, 0.44);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/0.95);

    BOOST_CHECK(!changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK_CLOSE(priVars[Indices::waterSwitchIdx], 1.0, 1e-12);
    BOOST_CHECK_SMALL(priVars[Indices::compositionSwitchIdx], 1e-12);
}