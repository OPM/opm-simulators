#include <config.h>

#include "TestTypeTag.hpp"

#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/blackoil/blackoilnewtonupdate.hpp>

#define BOOST_TEST_MODULE BlackoilPrimaryVariables
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
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

struct SaturationDependentCapillaryMaterialLaw
{
    struct Params {};

    template<class Container, class FluidState>
    static void capillaryPressures(Container& result, const Params&, const FluidState& state)
    {
        result[FluidSystem::waterPhaseIdx] = 1000.0 * state.saturation(FluidSystem::waterPhaseIdx);
        result[FluidSystem::oilPhaseIdx] = 2000.0 * state.saturation(FluidSystem::oilPhaseIdx);
        result[FluidSystem::gasPhaseIdx] = 3000.0 * state.saturation(FluidSystem::gasPhaseIdx);
    }

private:
    using FluidSystem = Opm::BlackOilFluidSystem<double>;
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

    // No DRSDT/DRVDT limits.
    double maxGasDissolutionFactor(unsigned, unsigned) const
    {
        return std::numeric_limits<double>::max();
    }

    double maxOilSaturation(unsigned) const
    {
        return maxOilSaturation_;
    }

    // Historical maximum oil saturation (VAPPARS).
    double maxOilSaturation_ = 1.0;

    double maxOilVaporizationFactor(unsigned, unsigned) const
    {
        return std::numeric_limits<double>::max();
    }

private:
    MaterialLawParams materialLawParams_{};
};

namespace Opm::Properties {

namespace TTag {

struct ClassicalBlackoilPrimaryVariablesTestTypeTag {
    using InheritsFrom = std::tuple<TestTypeTag>;
};

struct CapillaryBlackoilPrimaryVariablesTestTypeTag {
    using InheritsFrom = std::tuple<ClassicalBlackoilPrimaryVariablesTestTypeTag>;
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

template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::CapillaryBlackoilPrimaryVariablesTestTypeTag>
{
    using type = ::SaturationDependentCapillaryMaterialLaw;
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
    FluidSystemGuard(const FluidSystemGuard&) = delete;
    FluidSystemGuard& operator=(const FluidSystemGuard&) = delete;
    FluidSystemGuard(FluidSystemGuard&&) = delete;
    FluidSystemGuard& operator=(FluidSystemGuard&&) = delete;

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

// Typical black-oil values in SI units: Rs and Rv in sm3/sm3, Bo and Bg in
// rm3/sm3 (Bg at a reservoir pressure of about 250 bar).  The PVT objects
// below are independent of pressure.
constexpr double typicalRs = 100.0;
constexpr double typicalRv = 2.0e-4;
constexpr double typicalBo = 1.3;
constexpr double typicalBg = 0.005;

std::shared_ptr<OilPvt> makeConstantRsOilPvt(double rsSat, double bo = typicalBo)
{
    auto oilPvt = std::make_shared<OilPvt>();
    oilPvt->setApproach(Opm::OilPvtApproach::ConstantRsDeadOil);
    auto& realOilPvt = oilPvt->template getRealPvt<Opm::OilPvtApproach::ConstantRsDeadOil>();
    realOilPvt.setNumRegions(1);
    realOilPvt.setConstantRs(rsSat);
    realOilPvt.setBubblePointPressure(0.0);
    using TabulatedFunction = typename std::decay_t<decltype(realOilPvt)>::TabulatedOneDFunction;
    realOilPvt.setInverseOilFormationVolumeFactor(0, TabulatedFunction(2,
                                                                       std::vector<double>{1.0e5, 2.0e5},
                                                                       std::vector<double>{1.0 / bo, 1.0 / bo}));
    return oilPvt;
}

std::shared_ptr<GasPvt> makeWetGasPvt(double rvSat, double bg = typicalBg)
{
    auto gasPvt = std::make_shared<GasPvt>();
    gasPvt->setApproach(Opm::GasPvtApproach::WetGas);
    auto& realGasPvt = gasPvt->template getRealPvt<Opm::GasPvtApproach::WetGas>();
    realGasPvt.setNumRegions(1);
    realGasPvt.setSaturatedGasOilVaporizationFactor(0,
                                                    {{1.0e5, rvSat},
                                                     {2.0e5, rvSat}});

    // Constant gas formation volume factor, independent of pressure and Rv.
    using RealGasPvt = std::decay_t<decltype(realGasPvt)>;
    typename RealGasPvt::TabulatedTwoDFunctionBuilder invBg;
    for (const double p : {1.0e5, 2.0e5}) {
        const auto i = invBg.appendXPos(p);
        invBg.appendSamplePoint(i, 0.0, 1.0 / bg);
        invBg.appendSamplePoint(i, 1.0, 1.0 / bg);
    }
    realGasPvt.setInverseGasFormationVolumeFactor(0, std::move(invBg).build());
    return gasPvt;
}

// Keep the supplied instance's temperature distinct from the static system's.
// A temperature-dependent PVT makes mixing the two systems observable in both
// the resulting Rs value and the decision whether a gas phase appears.
struct TemperatureDependentFluidSystem : FluidSystem
{
    struct TemperatureDependentOilPvt
    {
        double saturatedGasDissolutionFactor(unsigned, double temperature,
                                             double, double, double) const
        {
            return 0.001 * temperature;
        }
    };

    double reservoirTemperature(unsigned = 0) const { return 350.0; }
    TemperatureDependentOilPvt oilPvt() const { return {}; }
};

} // namespace

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesUsesExplicitFluidSystemTemperature)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(0.25));
    Problem problem;
    TemperatureDependentFluidSystem supplied;

    // Disappearing gas must use Rs_sat(350 K) = 0.35 from the supplied system.
    // Reading the static temperature instead produces Rs_sat(300 K) = 0.30.
    auto disappearing = makePrimaryVariables(
        PrimaryVariables::GasMeaning::Sg, 0.2, -0.01,
        PrimaryVariables::WaterMeaning::Sw, PrimaryVariables::PressureMeaning::Po, 1.e5);
    auto legacy = disappearing;
    BOOST_REQUIRE(disappearing.adaptPrimaryVariables(problem, supplied, 0, 1.0, 1.0));
    BOOST_CHECK(disappearing.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_CLOSE(disappearing[Indices::compositionSwitchIdx], 0.35, 1.e-12);

    // The legacy overload continues to use the static system's PVT.
    BOOST_REQUIRE(legacy.adaptPrimaryVariables(problem, 0, 1.0, 1.0));
    BOOST_CHECK_CLOSE(legacy[Indices::compositionSwitchIdx], 0.25, 1.e-12);

    // Rs = 0.325 is below the instance's saturation limit of 0.35, so gas
    // must stay absent. The wrong temperature lowers the limit to 0.30 and
    // incorrectly switches this cell to gas saturation.
    auto undersaturated = makePrimaryVariables(
        PrimaryVariables::GasMeaning::Rs, 0.2, 0.325,
        PrimaryVariables::WaterMeaning::Sw, PrimaryVariables::PressureMeaning::Po, 1.e5);
    BOOST_CHECK(!undersaturated.adaptPrimaryVariables(problem, supplied, 0, 1.0, 1.0));
    BOOST_CHECK(undersaturated.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_CLOSE(undersaturated[Indices::compositionSwitchIdx], 0.325, 1.e-12);
}

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
    const double rsSat = typicalRs;
    const double Bo = typicalBo;
    const double Bg = typicalBg;
    FluidSystem::setOilPvt(makeConstantRsOilPvt(rsSat, Bo));
    FluidSystem::setGasPvt(makeWetGasPvt(0.0, Bg));

    Problem problem;
    const double sw = 0.95;
    const double sg = -5.0e-4;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        sw,
                                        sg,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK(priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Po);

    // At a small oil saturation, Rs is initialised from the gas/oil component
    // ratio, i.e. the saturated oil's dissolved gas reduced by the "negative"
    // free gas.
    const double so = 1.0 - sw - sg;
    const double expectedRs = (rsSat * so / Bo + sg / Bg) / (so / Bo);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], expectedRs, 1e-10);
    BOOST_CHECK(priVars[Indices::compositionSwitchIdx] < rsSat);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsKeepsSaturatedRsWithoutGasDeficit)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(typicalRs));
    FluidSystem::setGasPvt(makeWetGasPvt(0.0));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        0.95,
                                        -1.0e-14,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], typicalRs, 1e-6);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsLimitsRsToZero)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(typicalRs));
    FluidSystem::setGasPvt(makeWetGasPvt(0.0));

    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        0.97,
                                        -0.05,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_SMALL(priVars[Indices::compositionSwitchIdx], 1e-12);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsKeepsSaturatedRsAtLargeOilSaturation)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    const double rsSat = typicalRs;
    FluidSystem::setOilPvt(makeConstantRsOilPvt(rsSat));
    FluidSystem::setGasPvt(makeWetGasPvt(0.0));

    // A Newton overshoot of Sg in a cell with much oil is not taken as a gas
    // deficit; the oil is initialised as saturated.
    Problem problem;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        0.2,
                                        -0.2,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], rsSat, 1e-10);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsAtOilSaturationLimit)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    const double rsSat = typicalRs;
    const double Bo = typicalBo;
    const double Bg = typicalBg;
    FluidSystem::setOilPvt(makeConstantRsOilPvt(rsSat, Bo));
    FluidSystem::setGasPvt(makeWetGasPvt(0.0, Bg));

    // The component-based Rs is used for oil saturations (1 - Sw - Sg) below
    // 0.1, and the saturated Rs from 0.1 on.
    const double sg = -5.0e-4;
    auto rsAfterSwitch = [&](const double so) {
        Problem problem;
        auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                            1.0 - so - sg,
                                            sg,
                                            PrimaryVariables::WaterMeaning::Sw,
                                            PrimaryVariables::PressureMeaning::Po,
                                            1.5e5);
        const bool changed = priVars.adaptPrimaryVariables(problem,
                                                           /*globalDofIdx=*/0,
                                                           /*swMaximum=*/1.0,
                                                           /*thresholdWaterFilledCell=*/1.1);
        BOOST_CHECK(changed);
        BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
        return priVars[Indices::compositionSwitchIdx];
    };

    const double soBelow = 0.0999;
    BOOST_CHECK_CLOSE(rsAfterSwitch(soBelow), (rsSat * soBelow / Bo + sg / Bg) / (soBelow / Bo), 1e-10);
    BOOST_CHECK_CLOSE(rsAfterSwitch(0.1001), rsSat, 1e-10);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsAccountsForVaporizedOil)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(true);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    const double rsSat = typicalRs;
    const double rvSat = typicalRv;
    const double Bo = typicalBo;
    const double Bg = typicalBg;
    FluidSystem::setOilPvt(makeConstantRsOilPvt(rsSat, Bo));
    FluidSystem::setGasPvt(makeWetGasPvt(rvSat, Bg));

    Problem problem;
    const double sw = 0.95;
    const double sg = -5.0e-4;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        sw,
                                        sg,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);

    const double so = 1.0 - sw - sg;
    const double gasComponent = rsSat * so / Bo + sg / Bg;
    const double oilComponent = so / Bo + rvSat * sg / Bg;
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], gasComponent / oilComponent, 1e-10);
}

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesSwitchToRsUsesHistoryDependentRv)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(true);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    const double rsSat = typicalRs;
    const double rvSat = typicalRv;
    const double Bo = typicalBo;
    const double Bg = typicalBg;
    const double vapPar1 = 0.5;
    FluidSystem::setOilPvt(makeConstantRsOilPvt(rsSat, Bo));
    auto gasPvt = makeWetGasPvt(rvSat, Bg);
    gasPvt->template getRealPvt<Opm::GasPvtApproach::WetGas>().setVapPars(vapPar1, 0.0);
    FluidSystem::setGasPvt(gasPvt);

    // The oil saturation is below its historical maximum, so VAPPARS
    // reduces the Rv of the saturated gas as in the intensive quantities.
    Problem problem;
    problem.maxOilSaturation_ = 0.9;
    const double sw = 0.95;
    const double sg = -5.0e-4;
    auto priVars = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                        sw,
                                        sg,
                                        PrimaryVariables::WaterMeaning::Sw,
                                        PrimaryVariables::PressureMeaning::Po,
                                        1.5e5);

    const bool changed = priVars.adaptPrimaryVariables(problem,
                                                       /*globalDofIdx=*/0,
                                                       /*swMaximum=*/1.0,
                                                       /*thresholdWaterFilledCell=*/1.1);

    BOOST_CHECK(changed);
    BOOST_CHECK(priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);

    const double so = 1.0 - sw - sg;
    const double rv = rvSat * std::pow(so / problem.maxOilSaturation_, vapPar1);
    auto componentRatio = [&](const double rvUsed) {
        return (rsSat * so / Bo + sg / Bg) / (so / Bo + rvUsed * sg / Bg);
    };
    BOOST_CHECK_CLOSE(priVars[Indices::compositionSwitchIdx], componentRatio(rv), 1e-10);

    // Not vacuous: ignoring VAPPARS gives a different Rs.
    BOOST_CHECK_GT(std::abs(componentRatio(rv) - componentRatio(rvSat)), 1.0e-6);
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

BOOST_AUTO_TEST_CASE(AdaptPrimaryVariablesConvertsPressureWithNonzeroCapillaryPressure)
{
    using CapillaryTypeTag = Opm::Properties::TTag::CapillaryBlackoilPrimaryVariablesTestTypeTag;
    using CapillaryPrimaryVariables = Opm::BlackOilPrimaryVariables<CapillaryTypeTag>;
    using CapillaryProblem = Opm::GetPropType<CapillaryTypeTag, Opm::Properties::Problem>;

    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(false);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(true);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setGasPvt(makeWetGasPvt(0.15));
    CapillaryProblem problem;

    // Oil disappears: the capillary state has (Sw, So, Sg) = (0.4, 0, 0.7).
    // Pg = Po + Pc_g - Pc_o = 150000 + 2100 - 0.
    CapillaryPrimaryVariables disappearing(makePrimaryVariables(
        PrimaryVariables::GasMeaning::Sg, 0.4, 0.7,
        PrimaryVariables::WaterMeaning::Sw, PrimaryVariables::PressureMeaning::Po, 1.5e5));
    BOOST_REQUIRE(disappearing.adaptPrimaryVariables(problem, 0, 1.0, 1.1));
    BOOST_CHECK(disappearing.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv);
    BOOST_CHECK(disappearing.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg);
    BOOST_CHECK_CLOSE(disappearing[Indices::pressureSwitchIdx], 152100.0, 1e-12);
    BOOST_CHECK_CLOSE(disappearing[Indices::compositionSwitchIdx], 0.15, 1e-12);
    BOOST_CHECK_CLOSE(disappearing[Indices::waterSwitchIdx], 0.4, 1e-12);

    // Oil reappears: (Sw, So, Sg) = (0.4, 0, 0.6) at the phase boundary.
    // Po = Pg + Pc_o - Pc_g = 150000 + 0 - 1800.
    CapillaryPrimaryVariables appearing(makePrimaryVariables(
        PrimaryVariables::GasMeaning::Rv, 0.4, 0.3,
        PrimaryVariables::WaterMeaning::Sw, PrimaryVariables::PressureMeaning::Pg, 1.5e5));
    BOOST_REQUIRE(appearing.adaptPrimaryVariables(problem, 0, 1.0, 1.1));
    BOOST_CHECK(appearing.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK(appearing.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Po);
    BOOST_CHECK_CLOSE(appearing[Indices::pressureSwitchIdx], 148200.0, 1e-12);
    BOOST_CHECK_CLOSE(appearing[Indices::compositionSwitchIdx], 0.6, 1e-12);
    BOOST_CHECK_CLOSE(appearing[Indices::waterSwitchIdx], 0.4, 1e-12);
}

// Exercise the shared update through real primary variables, including the
// switch-history threshold supplied by the nonlinear method.
BOOST_AUTO_TEST_CASE(NewtonUpdateLimitsPressureAndCombinedSaturationChange)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    Problem problem;
    Opm::BlackoilNewtonParams<double> params{};
    params.dpMaxRel_ = 0.3;
    params.dsMax_ = 0.2;
    params.pressMin_ = 0.0;
    params.pressMax_ = 1.e9;
    params.waterSaturationMax_ = 1.0;
    params.waterOnlyThreshold_ = 1.0;

    const auto current = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                              0.3, 0.3,
                                              PrimaryVariables::WaterMeaning::Sw,
                                              PrimaryVariables::PressureMeaning::Po,
                                              1.e5);
    auto next = current;
    Opm::GetPropType<TypeTag, Opm::Properties::EqVector> correction(0.0);
    correction[Indices::pressureSwitchIdx] = 8.e4;
    correction[Indices::waterSwitchIdx] = 0.3;
    correction[Indices::compositionSwitchIdx] = 0.3;
    const bool switched = Opm::BlackOilNewtonUpdate<TypeTag>::update(
        problem, FluidSystem{}, 0, next, current, correction, params, false);

    BOOST_CHECK(!switched);
    BOOST_CHECK_CLOSE(next[Indices::pressureSwitchIdx], 7.e4, 1.e-12);
    // The oil saturation changes by the sum of the other two changes.
    BOOST_CHECK_CLOSE(next[Indices::waterSwitchIdx], 0.2, 1.e-12);
    BOOST_CHECK_CLOSE(next[Indices::compositionSwitchIdx], 0.2, 1.e-12);
    BOOST_CHECK_CLOSE(current[Indices::waterSwitchIdx], 0.3, 1.e-12);
}

BOOST_AUTO_TEST_CASE(NewtonUpdateUsesSwitchHistoryToSuppressOscillation)
{
    FluidSystemGuard guard;
    FluidSystem::setEnableDissolvedGas(true);
    FluidSystem::setEnableDissolvedGasInWater(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setEnableVaporizedWater(false);
    FluidSystem::setReservoirTemperature(300.0);
    FluidSystem::setOilPvt(makeConstantRsOilPvt(0.25));
    Problem problem;
    Opm::BlackoilNewtonParams<double> params{};
    params.dpMaxRel_ = 0.3;
    params.dsMax_ = 0.2;
    params.pressMin_ = 0.0;
    params.pressMax_ = 1.e9;
    params.waterSaturationMax_ = 1.0;
    params.waterOnlyThreshold_ = 1.0;
    params.priVarOscilationThreshold_ = 1.e-5;
    const auto current = makePrimaryVariables(PrimaryVariables::GasMeaning::Sg,
                                              0.2, 1.e-6,
                                              PrimaryVariables::WaterMeaning::Sw,
                                              PrimaryVariables::PressureMeaning::Po,
                                              1.e5);
    Opm::GetPropType<TypeTag, Opm::Properties::EqVector> correction(0.0);
    correction[Indices::compositionSwitchIdx] = 2.e-6;
    auto next = current;
    BOOST_CHECK(Opm::BlackOilNewtonUpdate<TypeTag>::update(
        problem, FluidSystem{}, 0, next, current, correction, params, false));
    BOOST_CHECK(next.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs);
    BOOST_CHECK_CLOSE(next[Indices::compositionSwitchIdx], 0.25, 1.e-12);

    next = current;
    BOOST_CHECK(!Opm::BlackOilNewtonUpdate<TypeTag>::update(
        problem, FluidSystem{}, 0, next, current, correction, params, true));
    BOOST_CHECK(next.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg);
    BOOST_CHECK_CLOSE(next[Indices::compositionSwitchIdx], -1.e-6, 1.e-12);
}
