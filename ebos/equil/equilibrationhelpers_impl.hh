// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <ebos/equil/equilibrationhelpers.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <fmt/format.h>

namespace Opm {
namespace EQUIL {

using FluidSystemSimple = BlackOilFluidSystem<double>;

// Adjust oil pressure according to gas saturation and cap pressure
using SatOnlyFluidState = SimpleModularFluidState<double,
                                                  /*numPhases=*/3,
                                                  /*numComponents=*/3,
                                                  FluidSystemSimple,
                                                  /*storePressure=*/false,
                                                  /*storeTemperature=*/false,
                                                  /*storeComposition=*/false,
                                                  /*storeFugacity=*/false,
                                                  /*storeSaturation=*/true,
                                                  /*storeDensity=*/false,
                                                  /*storeViscosity=*/false,
                                                  /*storeEnthalpy=*/false>;

namespace Miscibility {

template<class FluidSystem>
RsVD<FluidSystem>::RsVD(const int pvtRegionIdx,
                        const std::vector<double>& depth,
                        const std::vector<double>& rs)
    : pvtRegionIdx_(pvtRegionIdx)
    , rsVsDepth_(depth, rs)
{
}

template<class FluidSystem>
double RsVD<FluidSystem>::
operator()(const double depth,
           const double press,
           const double temp,
           const double satGas) const
{
    const auto sat_rs = satRs(press, temp);
    if (satGas > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return sat_rs;
    }
    else {
        if (rsVsDepth_.xMin() > depth)
            return std::min(sat_rs, rsVsDepth_.valueAt(0));
        else if (rsVsDepth_.xMax() < depth)
            return std::min(sat_rs, rsVsDepth_.valueAt(rsVsDepth_.numSamples() - 1));
        return std::min(sat_rs, rsVsDepth_.eval(depth, /*extrapolate=*/false));
    }
}

template<class FluidSystem>
double RsVD<FluidSystem>::satRs(const double press, const double temp) const
{
    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
PBVD<FluidSystem>::PBVD(const int pvtRegionIdx,
                        const std::vector<double>& depth,
                        const std::vector<double>& pbub)
    : pvtRegionIdx_(pvtRegionIdx)
    , pbubVsDepth_(depth, pbub)
{
}

template<class FluidSystem>
double PBVD<FluidSystem>::
operator()(const double depth,
           const double cellPress,
           const double temp,
           const double satGas) const
{
    double press = cellPress;
    if (satGas <= 0.0) {
        if (pbubVsDepth_.xMin() > depth)
            press = pbubVsDepth_.valueAt(0);
        else if (pbubVsDepth_.xMax() < depth)
            press = pbubVsDepth_.valueAt(pbubVsDepth_.numSamples() - 1);
        else
            press = pbubVsDepth_.eval(depth, /*extrapolate=*/false);
    }
    return satRs(std::min(press, cellPress), temp);
}

template<class FluidSystem>
double PBVD<FluidSystem>::
satRs(const double press, const double temp) const
{
    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
PDVD<FluidSystem>::PDVD(const int pvtRegionIdx,
                        const std::vector<double>& depth,
                        const std::vector<double>& pdew)
    : pvtRegionIdx_(pvtRegionIdx)
    , pdewVsDepth_(depth, pdew)
{
}

template<class FluidSystem>
double PDVD<FluidSystem>::
operator()(const double depth,
           const double cellPress,
           const double temp,
           const double satOil) const
{
    double press = cellPress;
    if (satOil <= 0.0) {
        if (pdewVsDepth_.xMin() > depth)
            press = pdewVsDepth_.valueAt(0);
        else if (pdewVsDepth_.xMax() < depth)
            press = pdewVsDepth_.valueAt(pdewVsDepth_.numSamples() - 1);
        else
            press = pdewVsDepth_.eval(depth, /*extrapolate=*/false);
    }
    return satRv(std::min(press, cellPress), temp);
}

template<class FluidSystem>
double PDVD<FluidSystem>::
satRv(const double press, const double temp) const
{
    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
}


template<class FluidSystem>
RvVD<FluidSystem>::RvVD(const int pvtRegionIdx,
                        const std::vector<double>& depth,
                        const std::vector<double>& rv)
    : pvtRegionIdx_(pvtRegionIdx)
    , rvVsDepth_(depth, rv)
{
}

template<class FluidSystem>
double RvVD<FluidSystem>::
operator()(const double depth,
           const double press,
           const double temp,
           const double satOil) const
{
    if (satOil < - std::sqrt(std::numeric_limits<double>::epsilon())) {
        throw std::logic_error {
            "Must not pass negative oil saturation"
        };
    }
    const auto sat_rv = satRv(press, temp);
    if (satOil > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return sat_rv;
    }
    else {
        if (rvVsDepth_.xMin() > depth)
            return std::min(sat_rv, rvVsDepth_.valueAt(0));
        else if (rvVsDepth_.xMax() < depth)
            return std::min(sat_rv, rvVsDepth_.valueAt(rvVsDepth_.numSamples() - 1));
        return std::min(sat_rv, rvVsDepth_.eval(depth, /*extrapolate=*/false));
    }
}

template<class FluidSystem>
double RvVD<FluidSystem>::
satRv(const double press, const double temp) const
{
    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
}


template<class FluidSystem>
RvwVD<FluidSystem>::RvwVD(const int pvtRegionIdx,
                        const std::vector<double>& depth,
                        const std::vector<double>& rvw)
    : pvtRegionIdx_(pvtRegionIdx)
    , rvwVsDepth_(depth, rvw)
{
}

template<class FluidSystem>
double RvwVD<FluidSystem>::
operator()(const double depth,
           const double press,
           const double temp,
           const double satWat) const
{
    if (satWat < - std::sqrt(std::numeric_limits<double>::epsilon())) {
        throw std::logic_error {
            "Must not pass negative water saturation"
        };
    }

    const auto sat_rvw = satRvw(press, temp);
    if (satWat > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return sat_rvw; //saturated Rvw
    }
    else {
        if (rvwVsDepth_.xMin() > depth)
            return std::min(sat_rvw,rvwVsDepth_.valueAt(0));
        else if (rvwVsDepth_.xMax() < depth)
            return std::min(sat_rvw, rvwVsDepth_.valueAt(rvwVsDepth_.numSamples() - 1));
        return std::min(sat_rvw, rvwVsDepth_.eval(depth, /*extrapolate=*/false));
    }
}

template<class FluidSystem>
double RvwVD<FluidSystem>::
satRvw(const double press, const double temp) const
{
    return FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press);
}


template<class FluidSystem>
RsSatAtContact<FluidSystem>::
RsSatAtContact(const int pvtRegionIdx, const double pContact,  const double T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rsSatContact_ = satRs(pContact, T_contact);
}

template<class FluidSystem>
double RsSatAtContact<FluidSystem>::
operator()(const double /* depth */,
           const double press,
           const double temp,
           const double satGas) const
{
    if (satGas > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return satRs(press, temp);
    }
    else {
        return std::min(satRs(press, temp), rsSatContact_);
    }
}

template<class FluidSystem>
double RsSatAtContact<FluidSystem>::
satRs(const double press, const double temp) const
{
    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
RvSatAtContact<FluidSystem>::
RvSatAtContact(const int pvtRegionIdx, const double pContact, const double T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rvSatContact_ = satRv(pContact, T_contact);
}

template<class FluidSystem>
double RvSatAtContact<FluidSystem>::
operator()(const double /*depth*/,
           const double press,
           const double temp,
           const double satOil) const
{
    if (satOil > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return satRv(press, temp);
    }
    else {
        return std::min(satRv(press, temp), rvSatContact_);
    }
}

template<class FluidSystem>
double RvSatAtContact<FluidSystem>::
satRv(const double press, const double temp) const
{
    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);;
}

template<class FluidSystem>
RvwSatAtContact<FluidSystem>::
RvwSatAtContact(const int pvtRegionIdx, const double pContact, const double T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rvwSatContact_ = satRvw(pContact, T_contact);
}

template<class FluidSystem>
double RvwSatAtContact<FluidSystem>::
operator()(const double /*depth*/,
           const double press,
           const double temp,
           const double satWat) const
{
    if (satWat > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return satRvw(press, temp);
    }
    else {
        return std::min(satRvw(press, temp), rvwSatContact_);
    }
}

template<class FluidSystem>
double RvwSatAtContact<FluidSystem>::
satRvw(const double press, const double temp) const
{
    return FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press);;
}

} // namespace Miscibility

EquilReg::EquilReg(const EquilRecord& rec,
                   std::shared_ptr<Miscibility::RsFunction> rs,
                   std::shared_ptr<Miscibility::RsFunction> rv,
                   std::shared_ptr<Miscibility::RsFunction> rvw,
                   const TabulatedFunction& tempVdTable,
                   const TabulatedFunction& saltVdTable,
                   const int pvtIdx)
    : rec_    (rec)
    , rs_     (rs)
    , rv_     (rv)
    , rvw_     (rvw)
    , tempVdTable_ (tempVdTable)
    , saltVdTable_ (saltVdTable)
    , pvtIdx_ (pvtIdx)
{
}

double EquilReg::datum() const
{
    return this->rec_.datumDepth();
}

double EquilReg::pressure() const
{
    return this->rec_.datumDepthPressure();
}

double EquilReg::zwoc() const
{
    return this->rec_.waterOilContactDepth();
}

double EquilReg::pcowWoc() const
{
    return this->rec_.waterOilContactCapillaryPressure();
}

double EquilReg::zgoc() const
{
    return this->rec_.gasOilContactDepth();
}

double EquilReg::pcgoGoc() const
{
    return this->rec_.gasOilContactCapillaryPressure();
}

int EquilReg::equilibrationAccuracy() const
{
    return this->rec_.initializationTargetAccuracy();
}

const EquilReg::CalcDissolution&
EquilReg::dissolutionCalculator() const
{
    return *this->rs_;
}

const EquilReg::CalcEvaporation&
EquilReg::evaporationCalculator() const
{
    return *this->rv_;
}

const EquilReg::CalcWaterEvaporation&
EquilReg::waterEvaporationCalculator() const
{
    return *this->rvw_;
}

const EquilReg::TabulatedFunction&
EquilReg::saltVdTable() const
{
    return saltVdTable_;
}

const EquilReg::TabulatedFunction&
EquilReg::tempVdTable() const
{
    return tempVdTable_;
}

int EquilReg::pvtIdx() const
{
    return this->pvtIdx_;
}

template<class FluidSystem, class MaterialLawManager>
PcEq<FluidSystem,MaterialLawManager>::
PcEq(const MaterialLawManager& materialLawManager,
     const int phase,
     const int cell,
     const double targetPc)
    : materialLawManager_(materialLawManager),
      phase_(phase),
      cell_(cell),
      targetPc_(targetPc)
{
}

template<class FluidSystem, class MaterialLawManager>
double PcEq<FluidSystem,MaterialLawManager>::
operator()(double s) const
{
    const auto& matParams = materialLawManager_.materialLawParams(cell_);
    SatOnlyFluidState fluidState;
    fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
    fluidState.setSaturation(phase_, s);

    std::array<double, FluidSystem::numPhases> pc{0.0};
    using MaterialLaw = typename MaterialLawManager::MaterialLaw;
    MaterialLaw::capillaryPressures(pc, matParams, fluidState);
    double sign = (phase_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    double pcPhase = pc[FluidSystem::oilPhaseIdx] + sign *  pc[phase_];
    return pcPhase - targetPc_;
}

template<class FluidSystem, class MaterialLawManager>
PcEqSum<FluidSystem,MaterialLawManager>::
PcEqSum(const MaterialLawManager& materialLawManager,
        const int phase1,
        const int phase2,
        const int cell,
        const double targetPc)
    : materialLawManager_(materialLawManager),
      phase1_(phase1),
      phase2_(phase2),
      cell_(cell),
      targetPc_(targetPc)
{
}

template<class FluidSystem, class MaterialLawManager>
double PcEqSum<FluidSystem,MaterialLawManager>::
operator()(double s) const
{
    const auto& matParams = materialLawManager_.materialLawParams(cell_);
    SatOnlyFluidState fluidState;
    fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
    fluidState.setSaturation(phase1_, s);
    fluidState.setSaturation(phase2_, 1.0 - s);

    std::array<double, FluidSystem::numPhases> pc {0.0};

    using MaterialLaw = typename MaterialLawManager::MaterialLaw;
    MaterialLaw::capillaryPressures(pc, matParams, fluidState);
    double sign1 = (phase1_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    double pc1 = pc[FluidSystem::oilPhaseIdx] + sign1 *  pc[phase1_];
    double sign2 = (phase2_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    double pc2 = pc[FluidSystem::oilPhaseIdx] + sign2 *  pc[phase2_];
    return pc1 + pc2 - targetPc_;
}

template <class FluidSystem, class MaterialLawManager>
double minSaturations(const MaterialLawManager& materialLawManager,
                      const int phase, const int cell)
{
    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    // Find minimum and maximum saturations.
    switch(phase) {
    case FluidSystem::waterPhaseIdx:
        return scaledDrainageInfo.Swl;

    case FluidSystem::gasPhaseIdx:
        return scaledDrainageInfo.Sgl;

    case FluidSystem::oilPhaseIdx:
        throw std::runtime_error("Min saturation not implemented for oil phase.");

    default:
        throw std::runtime_error("Unknown phaseIdx .");
    }
    return -1.0;
}

template <class FluidSystem, class MaterialLawManager>
double maxSaturations(const MaterialLawManager& materialLawManager,
                      const int phase, const int cell)
{
    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    // Find minimum and maximum saturations.
    switch(phase) {
    case FluidSystem::waterPhaseIdx:
        return scaledDrainageInfo.Swu;

    case FluidSystem::gasPhaseIdx:
        return scaledDrainageInfo.Sgu;

    case FluidSystem::oilPhaseIdx:
        throw std::runtime_error("Max saturation not implemented for oil phase.");

    default:
        throw std::runtime_error("Unknown phaseIdx .");
    }
    return -1.0;
}

template <class FluidSystem, class MaterialLawManager>
double satFromPc(const MaterialLawManager& materialLawManager,
                 const int phase,
                 const int cell,
                 const double targetPc,
                 const bool increasing)
{
    // Find minimum and maximum saturations.
    double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
    double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

    // Create the equation f(s) = pc(s) - targetPc
    const PcEq<FluidSystem, MaterialLawManager> f(materialLawManager, phase, cell, targetPc);
    double f0 = f(s0);
    double f1 = f(s1);
    if (!std::isfinite(f0 + f1))
        throw std::logic_error(fmt::format("The capillary pressure values {} and {} are not finite", f0, f1));

    if (f0 <= 0.0)
        return s0;
    else if (f1 >= 0.0)
        return s1;

    const double tol = 1e-10;
    // should at least converge in 2 times bisection but some safety here:
    const int maxIter = -2*static_cast<int>(std::log2(tol)) + 10;
    int usedIterations = -1;
    const double root = RegulaFalsiBisection<ThrowOnError>::solve(f, s0, s1, maxIter, tol, usedIterations);
    return root;
}

template<class FluidSystem, class MaterialLawManager>
double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                       const int phase1,
                       const int phase2,
                       const int cell,
                       const double targetPc)
{
    // Find minimum and maximum saturations.
    double s0 = minSaturations<FluidSystem>(materialLawManager, phase1, cell);
    double s1 = maxSaturations<FluidSystem>(materialLawManager, phase1, cell);

    // Create the equation f(s) = pc1(s) + pc2(1-s) - targetPc
    const PcEqSum<FluidSystem, MaterialLawManager> f(materialLawManager, phase1, phase2, cell, targetPc);
    double f0 = f(s0);
    double f1 = f(s1);
    if (f0 <= 0.0)
        return s0;
    else if (f1 >= 0.0)
        return s1;

    assert(f0 > 0.0 && f1 < 0.0);
    const double tol = 1e-10;
    // should at least converge in 2 times bisection but some safety here:
    const int maxIter = -2*static_cast<int>(std::log2(tol)) + 10;
    int usedIterations = -1;
    const double root = RegulaFalsiBisection<ThrowOnError>::solve(f, s0, s1, maxIter, tol, usedIterations);
    return root;
}

template<class FluidSystem, class MaterialLawManager>
double satFromDepth(const MaterialLawManager& materialLawManager,
                    const double cellDepth,
                    const double contactDepth,
                    const int phase,
                    const int cell,
                    const bool increasing)
{
    const double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
    const double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

    if (cellDepth < contactDepth) {
        return s0;
    }
    else {
        return s1;
    }
}

template<class FluidSystem, class MaterialLawManager>
bool isConstPc(const MaterialLawManager& materialLawManager,
               const int phase,
               const int cell)
{
    // Create the equation f(s) = pc(s);
    const PcEq<FluidSystem, MaterialLawManager> f(materialLawManager, phase, cell, 0);
    const double f0 = f(minSaturations<FluidSystem>(materialLawManager, phase, cell));
    const double f1 = f(maxSaturations<FluidSystem>(materialLawManager, phase, cell));
    return std::abs(f0 - f1) < std::numeric_limits<double>::epsilon();
}

}
}
