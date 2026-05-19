#include <config.h>

#include "TestTypeTag.hpp"

#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#define BOOST_TEST_MODULE BlackoilPrimaryVariables
#include <boost/test/unit_test.hpp>

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

struct ZeroCapillaryMaterialLaw
{
    struct Params {};

    template<class Container, class FluidState>
    static void capillaryPressures(Container& result, const Params&, const FluidState&)
    {
        for (auto& value : result) {
            value = 0.0;
        }
    }
};

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
        return materialLawParams_;
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

private:
    MaterialLawParams materialLawParams_{};
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

template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::ClassicalBlackoilPrimaryVariablesTestTypeTag>
{
    using type = ::ZeroCapillaryMaterialLaw;
};

} // namespace Opm::Properties

namespace {

using TypeTag = Opm::Properties::TTag::ClassicalBlackoilPrimaryVariablesTestTypeTag;
using PrimaryVariables = Opm::BlackOilPrimaryVariables<TypeTag>;
using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
using Problem = Opm::GetPropType<TypeTag, Opm::Properties::Problem>;
using GasPvt = typename FluidSystem::GasPvt;
using OilPvt = typename FluidSystem::OilPvt;

static_assert(Opm::getPropValue<TypeTag, Opm::Properties::EnergyModuleType>()
              == Opm::EnergyModules::NoTemperature);

struct FluidSystemGuard
{
    FluidSystemGuard()
        : dissolvedGas(FluidSystem::enableDissolvedGas())
        , dissolvedGasInWater(FluidSystem::enableDissolvedGasInWater())
        , vaporizedOil(FluidSystem::enableVaporizedOil())
        , vaporizedWater(FluidSystem::enableVaporizedWater())
        , reservoirTemperature(FluidSystem::reservoirTemperature())
        , gasPvt(FluidSystem::gasPvt())
        , oilPvt(FluidSystem::oilPvt())
    {
    }

    ~FluidSystemGuard()
    {
        FluidSystem::setGasPvt(std::make_shared<GasPvt>(gasPvt));
        FluidSystem::setOilPvt(std::make_shared<OilPvt>(oilPvt));
        FluidSystem::setReservoirTemperature(reservoirTemperature);
        FluidSystem::setEnableVaporizedWater(vaporizedWater);
        FluidSystem::setEnableVaporizedOil(vaporizedOil);
        FluidSystem::setEnableDissolvedGasInWater(dissolvedGasInWater);
        FluidSystem::setEnableDissolvedGas(dissolvedGas);
    }

    bool dissolvedGas;
    bool dissolvedGasInWater;
    bool vaporizedOil;
    bool vaporizedWater;
    double reservoirTemperature;
    GasPvt gasPvt;
    OilPvt oilPvt;
};

PrimaryVariables makePrimaryVariables(PrimaryVariables::GasMeaning gasMeaning,
                                      double waterSwitch,
                                      double compositionSwitch,
                                      PrimaryVariables::WaterMeaning waterMeaning = PrimaryVariables::WaterMeaning::Sw,
                                      PrimaryVariables::PressureMeaning pressureMeaning = PrimaryVariables::PressureMeaning::Po,
                                      double pressure = 0.0)
{
    PrimaryVariables priVars;
    for (std::size_t i = 0; i < priVars.size(); ++i) {
        priVars[i] = 0.0;
    }

    priVars.setPrimaryVarsMeaningWater(waterMeaning);
    priVars.setPrimaryVarsMeaningGas(gasMeaning);
    priVars.setPrimaryVarsMeaningPressure(pressureMeaning);
    priVars.setPrimaryVarsMeaningBrine(PrimaryVariables::BrineMeaning::Disabled);
    priVars.setPrimaryVarsMeaningSolvent(PrimaryVariables::SolventMeaning::Disabled);
    priVars[Indices::waterSwitchIdx] = waterSwitch;
    priVars[Indices::compositionSwitchIdx] = compositionSwitch;
    priVars[Indices::pressureSwitchIdx] = pressure;
    return priVars;
}

std::shared_ptr<OilPvt> makeConstantRsOilPvt(double rsSat)
{
    auto oilPvt = std::make_shared<OilPvt>();
    oilPvt->setApproach(Opm::OilPvtApproach::ConstantRsDeadOil);
    auto& realOilPvt = oilPvt->template getRealPvt<Opm::OilPvtApproach::ConstantRsDeadOil>();
    realOilPvt.setNumRegions(1);
    realOilPvt.setConstantRs(rsSat);
    realOilPvt.setBubblePointPressure(0.0);
    return oilPvt;
}

std::shared_ptr<GasPvt> makeWetGasPvt(double rvSat)
{
    auto gasPvt = std::make_shared<GasPvt>();
    gasPvt->setApproach(Opm::GasPvtApproach::WetGas);
    auto& realGasPvt = gasPvt->template getRealPvt<Opm::GasPvtApproach::WetGas>();
    realGasPvt.setNumRegions(1);
    realGasPvt.setSaturatedGasOilVaporizationFactor(0,
                                                    {{1.0e5, rvSat},
                                                     {2.0e5, rvSat}});
    return gasPvt;
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

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchesFromGasSaturationToRs)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(0.25));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg, 0.2, -0.05);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK(priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Po);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], 0.25, 1e-12);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchesFromRsToGasSaturation)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(0.25));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Rs, 0.2, 0.3);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK(priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Po);
    BOOST_CHECK_SMALL(priVars[Indices::compositionSwitchIdx], 1e-12);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchesFromGasSaturationToRv)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(false);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(true);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setGasPvt(makeWetGasPvt(0.15));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        0.4,
                                        0.7,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv);
    BOOST_CHECK(priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], 0.15, 1e-12);
    BOOST_CHECK_CLOSE(priVars[Indices::pressureSwitchIdx], 1.5e5, 1e-12);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchesFromRvToGasSaturation)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(false);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(true);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setGasPvt(makeWetGasPvt(0.15));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Rv,
                                        0.2,
                                        0.3,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Pg,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK(priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Po);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], 0.8, 1e-12);
    BOOST_CHECK_CLOSE(priVars[Indices::pressureSwitchIdx], 1.5e5, 1e-12);
}