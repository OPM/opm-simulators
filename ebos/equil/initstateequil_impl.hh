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
#ifndef EWOMS_INITSTATEEQUIL_IMPL_HH
#define EWOMS_INITSTATEEQUIL_IMPL_HH


#include <ebos/equil/initstateequil.hh>
#include <ebos/equil/equilibrationhelpers.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/utility/RegionMapping.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RsvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RvvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RvwvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PbvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PdvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>

#include <opm/input/eclipse/EclipseState/Tables/SaltpvdTable.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <dune/grid/common/mcmgmapper.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace Opm {
namespace EQUIL {

namespace Details {

template <typename CellRange, typename Comm>
void verticalExtent(const CellRange&      cells,
                    const std::vector<std::pair<double, double>>& cellZMinMax,
                    const Comm& comm,
                    std::array<double,2>& span)
{
    span[0] = std::numeric_limits<double>::max();
    span[1] = std::numeric_limits<double>::lowest();

    // Define vertical span as
    //
    //   [minimum(node depth(cells)), maximum(node depth(cells))]
    //
    // Note: The implementation of 'RK4IVP<>' implicitly
    // imposes the requirement that cell centroids are all
    // within this vertical span.  That requirement is not
    // checked.
    for (const auto& cell : cells) {
        if (cellZMinMax[cell].first < span[0]) { span[0] = cellZMinMax[cell].first; }
        if (cellZMinMax[cell].second > span[1]) { span[1] = cellZMinMax[cell].second; }
    }
    span[0] = comm.min(span[0]);
    span[1] = comm.max(span[1]);
}

void subdivisionCentrePoints(const double                            left,
                             const double                            right,
                             const int                               numIntervals,
                             std::vector<std::pair<double, double>>& subdiv)
{
    const auto h = (right - left) / numIntervals;

    auto end = left;
    for (auto i = 0*numIntervals; i < numIntervals; ++i) {
        const auto start = end;
        end = left + (i + 1)*h;

        subdiv.emplace_back((start + end) / 2, h);
    }
}

template <typename CellID>
std::vector<std::pair<double, double>>
horizontalSubdivision(const CellID cell,
                      const std::pair<double, double> topbot,
                      const int    numIntervals)
{
    auto subdiv = std::vector<std::pair<double, double>>{};
    subdiv.reserve(2 * numIntervals);

    if (topbot.first > topbot.second) {
        throw std::out_of_range {
            "Negative thickness (inverted top/bottom faces) in cell "
            + std::to_string(cell)
        };
    }

    subdivisionCentrePoints(topbot.first, topbot.second,
                            2*numIntervals, subdiv);

    return subdiv;
}

template <class Element>
double cellCenterDepth(const Element& element)
{
    typedef typename Element::Geometry Geometry;
    static constexpr int zCoord = Element::dimension - 1;
    double zz = 0.0;

    const Geometry& geometry = element.geometry();
    const int corners = geometry.corners();
    for (int i=0; i < corners; ++i)
        zz += geometry.corner(i)[zCoord];

    return zz/corners;
}

template <class Element>
std::pair<double,double> cellZSpan(const Element& element)
{
    typedef typename Element::Geometry Geometry;
    static constexpr int zCoord = Element::dimension - 1;
    double bot = 0.0;
    double top = 0.0;

    const Geometry& geometry = element.geometry();
    const int corners = geometry.corners();
    assert(corners == 8);
    for (int i=0; i < 4; ++i)
        bot += geometry.corner(i)[zCoord];
    for (int i=4; i < corners; ++i)
        top += geometry.corner(i)[zCoord];

    return std::make_pair(bot/4, top/4);
}

template <class Element>
std::pair<double,double> cellZMinMax(const Element& element)
{
    typedef typename Element::Geometry Geometry;
    static constexpr int zCoord = Element::dimension - 1;
    const Geometry& geometry = element.geometry();
    const int corners = geometry.corners();
    assert(corners == 8);
    auto min = std::numeric_limits<double>::max();
    auto max = std::numeric_limits<double>::lowest();


    for (int i=0; i < corners; ++i) {
        min = std::min(min, geometry.corner(i)[zCoord]);
        max = std::max(max, geometry.corner(i)[zCoord]);
    }
    return std::make_pair(min, max);
}

template<class RHS>
RK4IVP<RHS>::RK4IVP(const RHS& f,
                    const std::array<double,2>& span,
                    const double y0,
                    const int N)
    : N_(N)
    , span_(span)
{
    const double h = stepsize();
    const double h2 = h / 2;
    const double h6 = h / 6;

    y_.reserve(N + 1);
    f_.reserve(N + 1);

    y_.push_back(y0);
    f_.push_back(f(span_[0], y0));

    for (int i = 0; i < N; ++i) {
        const double x = span_[0] + i*h;
        const double y = y_.back();

        const double k1 = f_[i];
        const double k2 = f(x + h2, y + h2*k1);
        const double k3 = f(x + h2, y + h2*k2);
        const double k4 = f(x + h, y + h*k3);

        y_.push_back(y + h6*(k1 + 2*(k2 + k3) + k4));
        f_.push_back(f(x + h, y_.back()));
    }

    assert (y_.size() == std::vector<double>::size_type(N + 1));
}

template<class RHS>
double RK4IVP<RHS>::
operator()(const double x) const
{
    // Dense output (O(h**3)) according to Shampine
    // (Hermite interpolation)
    const double h = stepsize();
    int i = (x - span_[0]) / h;
    const double t = (x - (span_[0] + i*h)) / h;

    // Crude handling of evaluation point outside "span_";
    if (i  <  0) { i = 0;      }
    if (N_ <= i) { i = N_ - 1; }

    const double y0 = y_[i], y1 = y_[i + 1];
    const double f0 = f_[i], f1 = f_[i + 1];

    double u = (1 - 2*t) * (y1 - y0);
    u += h * ((t - 1)*f0 + t*f1);
    u *= t * (t - 1);
    u += (1 - t)*y0 + t*y1;

    return u;
}

template<class RHS>
double RK4IVP<RHS>::
stepsize() const
{
    return (span_[1] - span_[0]) / N_;
}

namespace PhasePressODE {

template<class FluidSystem>
Water<FluidSystem>::
Water(const TabulatedFunction& tempVdTable,
      const TabulatedFunction& saltVdTable,
      const int pvtRegionIdx,
      const double normGrav)
    : tempVdTable_(tempVdTable)
    , saltVdTable_(saltVdTable)
    , pvtRegionIdx_(pvtRegionIdx)
    , g_(normGrav)
{
}

template<class FluidSystem>
double Water<FluidSystem>::
operator()(const double depth,
           const double press) const
{
    return this->density(depth, press) * g_;
}

template<class FluidSystem>
double Water<FluidSystem>::
density(const double depth,
        const double press) const
{
    // The initializing algorithm can give depths outside the range due to numerical noise i.e. we extrapolate
    double saltConcentration = saltVdTable_.eval(depth, /*extrapolate=*/true);
    double temp = tempVdTable_.eval(depth, /*extrapolate=*/true);
    double rho = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, 0.0 /*=Rsw*/, saltConcentration);
    rho *= FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx_);
    return rho;
}

template<class FluidSystem, class RS>
Oil<FluidSystem,RS>::
Oil(const TabulatedFunction& tempVdTable,
    const RS& rs,
    const int pvtRegionIdx,
    const double normGrav)
    : tempVdTable_(tempVdTable)
    , rs_(rs)
    , pvtRegionIdx_(pvtRegionIdx)
    , g_(normGrav)
{
}

template<class FluidSystem, class RS>
double Oil<FluidSystem,RS>::
operator()(const double depth,
           const double press) const
{
    return this->density(depth, press) * g_;
}

template<class FluidSystem, class RS>
double Oil<FluidSystem,RS>::
density(const double depth,
        const double press) const
{
    const double temp = tempVdTable_.eval(depth, /*extrapolate=*/true);
    double rs = 0.0;
    if(FluidSystem::enableDissolvedGas())
        rs = rs_(depth, press, temp);

    double bOil = 0.0;
    if (rs >= FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press)) {
        bOil = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp, press);
    }
    else {
        bOil = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, rs);
    }
    double rho = bOil * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx_);
    if (FluidSystem::enableDissolvedGas()) {
        rho += rs * bOil * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
    }

    return rho;
}

template<class FluidSystem, class RV, class RVW>
Gas<FluidSystem,RV,RVW>::
Gas(const TabulatedFunction& tempVdTable,
    const RV& rv,
    const RVW& rvw,
    const int pvtRegionIdx,
    const double normGrav)
    : tempVdTable_(tempVdTable)
    , rv_(rv)
    , rvw_(rvw)
    , pvtRegionIdx_(pvtRegionIdx)
    , g_(normGrav)
{
}

template<class FluidSystem, class RV, class RVW>
double Gas<FluidSystem,RV,RVW>::
operator()(const double depth,
           const double press) const
{
    return this->density(depth, press) * g_;
}

template<class FluidSystem, class RV, class RVW>
double Gas<FluidSystem,RV,RVW>::
density(const double depth,
        const double press) const
{
    const double temp = tempVdTable_.eval(depth, /*extrapolate=*/true);
    double rv = 0.0;
    if (FluidSystem::enableVaporizedOil())
        rv = rv_(depth, press, temp);

    double rvw = 0.0;
    if (FluidSystem::enableVaporizedWater())
        rvw = rvw_(depth, press, temp);

    double bGas = 0.0;

    if (FluidSystem::enableVaporizedOil() && FluidSystem::enableVaporizedWater()) {
        if (rv >= FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press)
            && rvw >= FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press))
        {
            bGas = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp, press);
        } else {
            bGas = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, rv, rvw);
        }
        double rho = bGas * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
        rho += rv * bGas * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx_)
               + rvw * bGas * FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx_);
        return rho;
    }

    if (FluidSystem::enableVaporizedOil()){
        if (rv >= FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press)) {
            bGas = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp, press);
        } else {
            bGas = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, rv, 0.0/*=rvw*/);
        }
        double rho = bGas * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
        rho += rv * bGas * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx_);
        return rho;
     }

    if (FluidSystem::enableVaporizedWater()){
        if (rvw >= FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press)) {
            bGas = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp, press);
        }
        else {
            bGas = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, 0.0/*=rv*/, rvw);
        }
        double rho = bGas * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
        rho += rvw * bGas * FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx_);
        return rho;
    }

    // immiscible gas
    bGas = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp, press, 0.0/*=rv*/, 0.0/*=rvw*/);
    double rho = bGas * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);

    return rho;
}

}

template<class FluidSystem, class Region>
template<class ODE>
PressureTable<FluidSystem,Region>::
PressureFunction<ODE>::PressureFunction(const ODE&      ode,
                                        const InitCond& ic,
                                        const int       nsample,
                                        const VSpan&    span)
    : initial_(ic)
{
    this->value_[Direction::Up] = std::make_unique<Distribution>
        (ode, VSpan {{ ic.depth, span[0] }}, ic.pressure, nsample);

    this->value_[Direction::Down] = std::make_unique<Distribution>
        (ode, VSpan {{ ic.depth, span[1] }}, ic.pressure, nsample);
}

template<class FluidSystem, class Region>
template<class ODE>
PressureTable<FluidSystem,Region>::
PressureFunction<ODE>::PressureFunction(const PressureFunction& rhs)
    : initial_(rhs.initial_)
{
    this->value_[Direction::Up] =
        std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

    this->value_[Direction::Down] =
        std::make_unique<Distribution>(*rhs.value_[Direction::Down]);
}

template<class FluidSystem, class Region>
template<class ODE>
typename PressureTable<FluidSystem,Region>::template PressureFunction<ODE>&
PressureTable<FluidSystem,Region>::
PressureFunction<ODE>::
operator=(const PressureFunction& rhs)
{
    this->initial_ = rhs.initial_;

    this->value_[Direction::Up] =
        std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

    this->value_[Direction::Down] =
        std::make_unique<Distribution>(*rhs.value_[Direction::Down]);

    return *this;
}

template<class FluidSystem, class Region>
template<class ODE>
typename PressureTable<FluidSystem,Region>::template PressureFunction<ODE>&
PressureTable<FluidSystem,Region>::
PressureFunction<ODE>::
operator=(PressureFunction&& rhs)
{
    this->initial_ = rhs.initial_;
    this->value_   = std::move(rhs.value_);

    return *this;
}

template<class FluidSystem, class Region>
template<class ODE>
double
PressureTable<FluidSystem,Region>::
PressureFunction<ODE>::
value(const double depth) const
{
    if (depth < this->initial_.depth) {
        // Value above initial condition depth.
        return (*this->value_[Direction::Up])(depth);
    }
    else if (depth > this->initial_.depth) {
        // Value below initial condition depth.
        return (*this->value_[Direction::Down])(depth);
    }
    else {
        // Value *at* initial condition depth.
        return this->initial_.pressure;
    }
}


template<class FluidSystem, class Region>
template<typename PressFunc>
void PressureTable<FluidSystem,Region>::
checkPtr(const PressFunc*   phasePress,
         const std::string& phaseName) const
{
    if (phasePress != nullptr) { return; }

    throw std::invalid_argument {
        "Phase pressure function for \"" + phaseName
        + "\" most not be null"
    };
}

template<class FluidSystem, class Region>
typename PressureTable<FluidSystem,Region>::Strategy
PressureTable<FluidSystem,Region>::
selectEquilibrationStrategy(const Region& reg) const
{
    if (!this->oilActive()) {
        if (reg.datum() > reg.zwoc()) { // Datum in water zone
            return &PressureTable::equil_WOG;
        }
        return &PressureTable::equil_GOW;
    }

    if (reg.datum() > reg.zwoc()) {      // Datum in water zone
        return &PressureTable::equil_WOG;
    }
    else if (reg.datum() < reg.zgoc()) { // Datum in gas zone
        return &PressureTable::equil_GOW;
    }
    else {                               // Datum in oil zone
        return &PressureTable::equil_OWG;
    }
}

template<class FluidSystem, class Region>
void  PressureTable<FluidSystem,Region>::
copyInPointers(const PressureTable& rhs)
{
    if (rhs.oil_ != nullptr) {
        this->oil_ = std::make_unique<OPress>(*rhs.oil_);
    }

    if (rhs.gas_ != nullptr) {
        this->gas_ = std::make_unique<GPress>(*rhs.gas_);
    }

    if (rhs.wat_ != nullptr) {
        this->wat_ = std::make_unique<WPress>(*rhs.wat_);
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
PhaseSaturations<MaterialLawManager,FluidSystem,Region,CellID>::
PhaseSaturations(MaterialLawManager& matLawMgr,
                 const std::vector<double>& swatInit)
    : matLawMgr_(matLawMgr)
    , swatInit_ (swatInit)
{
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
PhaseSaturations<MaterialLawManager,FluidSystem,Region,CellID>::
PhaseSaturations(const PhaseSaturations& rhs)
        : matLawMgr_(rhs.matLawMgr_)
        , swatInit_ (rhs.swatInit_)
        , sat_      (rhs.sat_)
        , press_    (rhs.press_)
{
    // Note: We don't need to do anything to the 'fluidState_' here.
    this->setEvaluationPoint(*rhs.evalPt_.position,
                             *rhs.evalPt_.region,
                             *rhs.evalPt_.ptable);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
const PhaseQuantityValue&
PhaseSaturations<MaterialLawManager,FluidSystem,Region,CellID>::
deriveSaturations(const Position& x,
                  const Region&   reg,
                  const PTable&   ptable)
{
    this->setEvaluationPoint(x, reg, ptable);
    this->initializePhaseQuantities();

    if (ptable.waterActive()) { this->deriveWaterSat(); }
    if (ptable.gasActive())   { this->deriveGasSat();   }

    if (this->isOverlappingTransition()) {
        this->fixUnphysicalTransition();
    }

    if (ptable.oilActive()) { this->deriveOilSat(); }

    this->accountForScaledSaturations();

    return this->sat_;
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager,FluidSystem,Region,CellID>::
setEvaluationPoint(const Position& x,
                   const Region&   reg,
                   const PTable&   ptable)
{
    this->evalPt_.position = &x;
    this->evalPt_.region   = &reg;
    this->evalPt_.ptable   = &ptable;
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager,FluidSystem,Region,CellID>::
initializePhaseQuantities()
{
    this->sat_.reset();
    this->press_.reset();

    const auto  depth  = this->evalPt_.position->depth;
    const auto& ptable = *this->evalPt_.ptable;

    if (ptable.oilActive()) {
        this->press_.oil = ptable.oil(depth);
    }

    if (ptable.gasActive()) {
        this->press_.gas = ptable.gas(depth);
    }

    if (ptable.waterActive()) {
        this->press_.water = ptable.water(depth);
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::deriveOilSat()
{
    this->sat_.oil = 1.0 - this->sat_.water - this->sat_.gas;
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::deriveGasSat()
{
    auto& sg = this->sat_.gas;

    const auto isIncr = true; // dPcgo/dSg >= 0 for all Sg.
    const auto oilActive = this->evalPt_.ptable->oilActive();

    if (this->isConstCapPress(this->gasPos())) {
        // Sharp interface between phases.  Can derive phase saturation
        // directly from knowing where 'depth' of evaluation point is
        // relative to depth of O/G contact.
        const auto gas_contact = oilActive? this->evalPt_.region->zgoc() : this->evalPt_.region->zwoc();
        sg = this->fromDepthTable(gas_contact,
                                  this->gasPos(), isIncr);
    }
    else {
        // Capillary pressure curve is non-constant, meaning there is a
        // transition zone between the gas and oil phases.  Invert capillary
        // pressure relation
        //
        //    Pcgo(Sg) = Pg - Po
        //
        // Note that Pcgo is defined to be (Pg - Po), not (Po - Pg).
        const auto pw = oilActive? this->press_.oil : this->press_.water;
        const auto pcgo = this->press_.gas - pw;
        sg = this->invertCapPress(pcgo, this->gasPos(), isIncr);
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::deriveWaterSat()
{
    auto& sw = this->sat_.water;

    const auto isIncr = false; // dPcow/dSw <= 0 for all Sw.

    if (this->isConstCapPress(this->waterPos())) {
        // Sharp interface between phases.  Can derive phase saturation
        // directly from knowing where 'depth' of evaluation point is
        // relative to depth of O/W contact.
        sw = this->fromDepthTable(this->evalPt_.region->zwoc(),
                                  this->waterPos(), isIncr);
    }
    else {
        // Capillary pressure curve is non-constant, meaning there is a
        // transition zone between the oil and water phases.  Invert
        // capillary pressure relation
        //
        //    Pcow(Sw) = Po - Pw
        //
        // unless the model uses "SWATINIT".  In the latter case, pick the
        // saturation directly from the SWATINIT array of the pertinent
        // cell.
        const auto pcow = this->press_.oil - this->press_.water;

        if (this->swatInit_.empty()) {
            sw = this->invertCapPress(pcow, this->waterPos(), isIncr);
        }
        else {
            auto [swout, newSwatInit] = this->applySwatInit(pcow);
            if (newSwatInit)
                sw = this->invertCapPress(pcow, this->waterPos(), isIncr);
            else {
                sw = swout;
            }
        }
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
fixUnphysicalTransition()
{
    auto& sg = this->sat_.gas;
    auto& sw = this->sat_.water;

    // Overlapping gas/oil and oil/water transition zones can lead to
    // unphysical phase saturations when individual saturations are derived
    // directly from inverting O/G and O/W capillary pressure curves.
    //
    // Recalculate phase saturations using the implied gas/water capillary
    // pressure: Pg - Pw.
    const auto pcgw = this->press_.gas - this->press_.water;
    if (! this->swatInit_.empty()) {
        // Re-scale Pc to reflect imposed sw for vanishing oil phase.  This
        // seems consistent with ECLIPSE, but fails to honour SWATINIT in
        // case of non-trivial gas/oil capillary pressure.
        auto [swout, newSwatInit] = this->applySwatInit(pcgw, sw);
        if (newSwatInit){
            const auto isIncr = false; // dPcow/dSw <= 0 for all Sw.
            sw = this->invertCapPress(pcgw, this->waterPos(), isIncr);
        }
        else {
            sw = swout;
        }
    }

    sw = satFromSumOfPcs<FluidSystem>
        (this->matLawMgr_, this->waterPos(), this->gasPos(),
         this->evalPt_.position->cell, pcgw);
    sg = 1.0 - sw;

    this->fluidState_.setSaturation(this->oilPos(), 1.0 - sw - sg);
    this->fluidState_.setSaturation(this->gasPos(), sg);
    this->fluidState_.setSaturation(this->waterPos(), this->evalPt_
                                    .ptable->waterActive() ? sw : 0.0);

    // Pcgo = Pg - Po => Po = Pg - Pcgo
    this->computeMaterialLawCapPress();
    this->press_.oil = this->press_.gas - this->materialLawCapPressGasOil();
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
accountForScaledSaturations()
{
    const auto gasActive = this->evalPt_.ptable->gasActive();
    const auto watActive = this->evalPt_.ptable->waterActive();
    const auto oilActive = this->evalPt_.ptable->oilActive();

    auto sg = gasActive? this->sat_.gas : 0.0;
    auto sw = watActive? this->sat_.water : 0.0;
    auto so = oilActive? this->sat_.oil : 0.0;

    this->fluidState_.setSaturation(this->waterPos(), sw);
    this->fluidState_.setSaturation(this->oilPos(), so);
    this->fluidState_.setSaturation(this->gasPos(), sg);

    const auto& scaledDrainageInfo = this->matLawMgr_
        .oilWaterScaledEpsInfoDrainage(this->evalPt_.position->cell);

    const auto thresholdSat = 1.0e-6;
    if (watActive && ((sw + thresholdSat) > scaledDrainageInfo.Swu)) {
        // Water saturation exceeds maximum possible value.  Reset oil phase
        // pressure to that which corresponds to maximum possible water
        // saturation value.
        this->fluidState_.setSaturation(this->waterPos(), scaledDrainageInfo.Swu);
        if (oilActive) {
            this->fluidState_.setSaturation(this->oilPos(), so + sw - scaledDrainageInfo.Swu);
        } else if (gasActive) {
            this->fluidState_.setSaturation(this->gasPos(), sg + sw - scaledDrainageInfo.Swu);
        }
        sw = scaledDrainageInfo.Swu;
        this->computeMaterialLawCapPress();

        if (oilActive) {
            // Pcow = Po - Pw => Po = Pw + Pcow
            this->press_.oil = this->press_.water + this->materialLawCapPressOilWater();
        } else {
            // Pcgw = Pg - Pw => Pg = Pw + Pcgw
            this->press_.gas = this->press_.water + this->materialLawCapPressGasWater();
        }

    }
    if (gasActive && ((sg + thresholdSat) > scaledDrainageInfo.Sgu)) {
        // Gas saturation exceeds maximum possible value.  Reset oil phase
        // pressure to that which corresponds to maximum possible gas
        // saturation value.
        this->fluidState_.setSaturation(this->gasPos(), scaledDrainageInfo.Sgu);
        if (oilActive) {
            this->fluidState_.setSaturation(this->oilPos(), so + sg - scaledDrainageInfo.Sgu);
        } else if (watActive) {
            this->fluidState_.setSaturation(this->waterPos(), sw + sg - scaledDrainageInfo.Sgu);
        }
        sg = scaledDrainageInfo.Sgu;
        this->computeMaterialLawCapPress();

        if (oilActive) {
            // Pcgo = Pg - Po => Po = Pg - Pcgo
            this->press_.oil = this->press_.gas - this->materialLawCapPressGasOil();
        } else {
            // Pcgw = Pg - Pw => Pw = Pg - Pcgw
            this->press_.water =  this->press_.gas - this->materialLawCapPressGasWater();
        }
    }

    if (watActive && ((sw - thresholdSat) < scaledDrainageInfo.Swl)) {
        // Water saturation less than minimum possible value in cell.  Reset
        // water phase pressure to that which corresponds to minimum
        // possible water saturation value.
        this->fluidState_.setSaturation(this->waterPos(), scaledDrainageInfo.Swl);
        if (oilActive) {
            this->fluidState_.setSaturation(this->oilPos(), so + sw - scaledDrainageInfo.Swl);
        } else if (gasActive) {
            this->fluidState_.setSaturation(this->gasPos(), sg + sw - scaledDrainageInfo.Swl);
        }
        sw = scaledDrainageInfo.Swl;
        this->computeMaterialLawCapPress();

        if (oilActive) {
            // Pcwo = Po - Pw => Pw = Po - Pcow
            this->press_.water = this->press_.oil - this->materialLawCapPressOilWater();
        } else {
            // Pcgw = Pg - Pw => Pw = Pg - Pcgw
            this->press_.water = this->press_.gas - this->materialLawCapPressGasWater();
        }
    }

    if (gasActive && ((sg - thresholdSat) < scaledDrainageInfo.Sgl)) {
        // Gas saturation less than minimum possible value in cell.  Reset
        // gas phase pressure to that which corresponds to minimum possible
        // gas saturation.
        this->fluidState_.setSaturation(this->gasPos(), scaledDrainageInfo.Sgl);
        if (oilActive) {
            this->fluidState_.setSaturation(this->oilPos(), so + sg - scaledDrainageInfo.Sgl);
        } else if (watActive) {
            this->fluidState_.setSaturation(this->waterPos(), sw + sg - scaledDrainageInfo.Sgl);
        }
        sg = scaledDrainageInfo.Sgl;
        this->computeMaterialLawCapPress();

        if (oilActive) {
            // Pcgo = Pg - Po => Pg = Po + Pcgo
            this->press_.gas = this->press_.oil + this->materialLawCapPressGasOil();
        } else {
            // Pcgw = Pg - Pw => Pg = Pw + Pcgw
            this->press_.gas = this->press_.water + this->materialLawCapPressGasWater();
        }
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
std::pair<double, bool>
PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
applySwatInit(const double pcow)
{
    return this->applySwatInit(pcow, this->swatInit_[this->evalPt_.position->cell]);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
std::pair<double, bool>
PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
applySwatInit(const double pcow, const double sw)
{
    return this->matLawMgr_.applySwatinit(this->evalPt_.position->cell, pcow, sw);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
void PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
computeMaterialLawCapPress()
{
    const auto& matParams = this->matLawMgr_
        .materialLawParams(this->evalPt_.position->cell);

    this->matLawCapPress_.fill(0.0);
    MaterialLaw::capillaryPressures(this->matLawCapPress_,
                                    matParams, this->fluidState_);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
materialLawCapPressGasOil() const
{
    return this->matLawCapPress_[this->oilPos()]
        + this->matLawCapPress_[this->gasPos()];
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
materialLawCapPressOilWater() const
{
    return this->matLawCapPress_[this->oilPos()]
        - this->matLawCapPress_[this->waterPos()];
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
materialLawCapPressGasWater() const
{
    return this->matLawCapPress_[this->gasPos()]
        - this->matLawCapPress_[this->waterPos()];
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
bool PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
isConstCapPress(const PhaseIdx phaseIdx) const
{
    return isConstPc<FluidSystem>
        (this->matLawMgr_, phaseIdx, this->evalPt_.position->cell);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
bool PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
isOverlappingTransition() const
{
    return this->evalPt_.ptable->gasActive()
        && this->evalPt_.ptable->waterActive()
        && ((this->sat_.gas + this->sat_.water) > 1.0);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
fromDepthTable(const double   contactdepth,
               const PhaseIdx phasePos,
               const bool     isincr) const
{
    return satFromDepth<FluidSystem>
        (this->matLawMgr_, this->evalPt_.position->depth,
         contactdepth, static_cast<int>(phasePos),
         this->evalPt_.position->cell, isincr);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
invertCapPress(const double   pc,
               const PhaseIdx phasePos,
               const bool     isincr) const
{
    return satFromPc<FluidSystem>
        (this->matLawMgr_, static_cast<int>(phasePos),
         this->evalPt_.position->cell, pc, isincr);
}

template<class FluidSystem, class Region>
PressureTable<FluidSystem,Region>::
PressureTable(const double gravity,
              const int    samplePoints)
    : gravity_(gravity)
    , nsample_(samplePoints)
{
}

template <class FluidSystem, class Region>
PressureTable<FluidSystem,Region>::
PressureTable(const PressureTable<FluidSystem,Region>& rhs)
    : gravity_(rhs.gravity_)
    , nsample_(rhs.nsample_)
{
    this->copyInPointers(rhs);
}

template <class FluidSystem, class Region>
PressureTable<FluidSystem,Region>::
PressureTable(PressureTable<FluidSystem,Region>&& rhs)
    : gravity_(rhs.gravity_)
    , nsample_(rhs.nsample_)
    , oil_    (std::move(rhs.oil_))
    , gas_    (std::move(rhs.gas_))
    , wat_    (std::move(rhs.wat_))
{
}

template <class FluidSystem, class Region>
PressureTable<FluidSystem,Region>&
PressureTable<FluidSystem,Region>::
operator=(const PressureTable<FluidSystem,Region>& rhs)
{
    this->gravity_ = rhs.gravity_;
    this->nsample_ = rhs.nsample_;
    this->copyInPointers(rhs);

    return *this;
}

template <class FluidSystem, class Region>
PressureTable<FluidSystem,Region>&
PressureTable<FluidSystem,Region>::
operator=(PressureTable<FluidSystem,Region>&& rhs)
{
    this->gravity_ = rhs.gravity_;
    this->nsample_ = rhs.nsample_;

    this->oil_ = std::move(rhs.oil_);
    this->gas_ = std::move(rhs.gas_);
    this->wat_ = std::move(rhs.wat_);

    return *this;
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem,Region>::
equilibrate(const Region& reg,
            const VSpan&  span)
{
    // One of the PressureTable::equil_*() member functions.
    auto equil = this->selectEquilibrationStrategy(reg);

    (this->*equil)(reg, span);
}

template <class FluidSystem, class Region>
bool PressureTable<FluidSystem,Region>::
oilActive() const
{
    return FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
}

template <class FluidSystem, class Region>
bool PressureTable<FluidSystem,Region>::
gasActive() const
{
    return FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
}

template <class FluidSystem, class Region>
bool PressureTable<FluidSystem,Region>::
waterActive() const
{
    return FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
}

template <class FluidSystem, class Region>
double PressureTable<FluidSystem,Region>::
oil(const double depth) const
{
    this->checkPtr(this->oil_.get(), "OIL");

    return this->oil_->value(depth);
}

template <class FluidSystem, class Region>
double PressureTable<FluidSystem,Region>::
gas(const double depth) const
{
    this->checkPtr(this->gas_.get(), "GAS");

    return this->gas_->value(depth);
}


template <class FluidSystem, class Region>
double PressureTable<FluidSystem,Region>::
water(const double depth) const
{
    this->checkPtr(this->wat_.get(), "WATER");

    return this->wat_->value(depth);
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
equil_WOG(const Region& reg, const VSpan& span)
{
    // Datum depth in water zone.  Calculate phase pressure for water first,
    // followed by oil and gas if applicable.

    if (! this->waterActive()) {
        throw std::invalid_argument {
            "Don't know how to interpret EQUIL datum depth in "
            "WATER zone in model without active water phase"
        };
    }

    {
        const auto ic = typename WPress::InitCond {
            reg.datum(), reg.pressure()
        };

        this->makeWatPressure(ic, reg, span);
    }

    if (this->oilActive()) {
        // Pcow = Po - Pw => Po = Pw + Pcow
        const auto ic = typename OPress::InitCond {
            reg.zwoc(),
            this->water(reg.zwoc()) + reg.pcowWoc()
        };

        this->makeOilPressure(ic, reg, span);
    }

    if (this->gasActive() && this->oilActive()) {
        // Pcgo = Pg - Po => Pg = Po + Pcgo
        const auto ic = typename GPress::InitCond {
            reg.zgoc(),
            this->oil(reg.zgoc()) + reg.pcgoGoc()
        };

        this->makeGasPressure(ic, reg, span);
    } else if (this->gasActive() && !this->oilActive()) {
        // No oil phase set Pg = Pw + Pcgw
        const auto ic = typename GPress::InitCond {
            reg.zwoc(), // The WOC is really the GWC for gas/water cases
            this->water(reg.zwoc()) + reg.pcowWoc() // Pcow(WOC) is really Pcgw(GWC) for gas/water cases
        };
        this->makeGasPressure(ic, reg, span);
    }
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
equil_GOW(const Region& reg, const VSpan& span)
{
    // Datum depth in gas zone.  Calculate phase pressure for gas first,
    // followed by oil and water if applicable.

    if (! this->gasActive()) {
        throw std::invalid_argument {
            "Don't know how to interpret EQUIL datum depth in "
            "GAS zone in model without active gas phase"
        };
    }

    {
        const auto ic = typename GPress::InitCond {
            reg.datum(), reg.pressure()
        };

        this->makeGasPressure(ic, reg, span);
    }

    if (this->oilActive()) {
        // Pcgo = Pg - Po => Po = Pg - Pcgo
        const auto ic = typename OPress::InitCond {
            reg.zgoc(),
            this->gas(reg.zgoc()) - reg.pcgoGoc()
        };
        this->makeOilPressure(ic, reg, span);
    }

    if (this->waterActive() && this->oilActive()) {
        // Pcow = Po - Pw => Pw = Po - Pcow
        const auto ic = typename WPress::InitCond {
            reg.zwoc(),
            this->oil(reg.zwoc()) - reg.pcowWoc()
        };

        this->makeWatPressure(ic, reg, span);
    } else if (this->waterActive() && !this->oilActive()) {
        // No oil phase set Pw = Pg - Pcgw
        const auto ic = typename WPress::InitCond {
            reg.zwoc(), // The WOC is really the GWC for gas/water cases
            this->gas(reg.zwoc()) - reg.pcowWoc() // Pcow(WOC) is really Pcgw(GWC) for gas/water cases
        };
        this->makeWatPressure(ic, reg, span);
    }
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
equil_OWG(const Region& reg, const VSpan& span)
{
    // Datum depth in oil zone.  Calculate phase pressure for oil first,
    // followed by gas and water if applicable.

    if (! this->oilActive()) {
        throw std::invalid_argument {
            "Don't know how to interpret EQUIL datum depth in "
            "OIL zone in model without active oil phase"
        };
    }

    {
        const auto ic = typename OPress::InitCond {
            reg.datum(), reg.pressure()
        };

        this->makeOilPressure(ic, reg, span);
    }

    if (this->waterActive()) {
        // Pcow = Po - Pw => Pw = Po - Pcow
        const auto ic = typename WPress::InitCond {
            reg.zwoc(),
            this->oil(reg.zwoc()) - reg.pcowWoc()
        };

        this->makeWatPressure(ic, reg, span);
    }

    if (this->gasActive()) {
        // Pcgo = Pg - Po => Pg = Po + Pcgo
        const auto ic = typename GPress::InitCond {
            reg.zgoc(),
            this->oil(reg.zgoc()) + reg.pcgoGoc()
        };
        this->makeGasPressure(ic, reg, span);
    }
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
makeOilPressure(const typename OPress::InitCond& ic,
                const Region&                    reg,
                const VSpan&                     span)
{
    const auto drho = OilPressODE {
        reg.tempVdTable(), reg.dissolutionCalculator(),
        reg.pvtIdx(), this->gravity_
    };

    this->oil_ = std::make_unique<OPress>(drho, ic, this->nsample_, span);
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
makeGasPressure(const typename GPress::InitCond& ic,
                const Region&                    reg,
                const VSpan&                     span)
{
    const auto drho = GasPressODE {
        reg.tempVdTable(), reg.evaporationCalculator(), reg.waterEvaporationCalculator(),
        reg.pvtIdx(), this->gravity_
    };

    this->gas_ = std::make_unique<GPress>(drho, ic, this->nsample_, span);
}

template <class FluidSystem, class Region>
void PressureTable<FluidSystem, Region>::
makeWatPressure(const typename WPress::InitCond& ic,
                const Region&                    reg,
                const VSpan&                     span)
{
    const auto drho = WatPressODE {
        reg.tempVdTable(), reg.saltVdTable(), reg.pvtIdx(), this->gravity_
    };

    this->wat_ = std::make_unique<WPress>(drho, ic, this->nsample_, span);
}

}

namespace DeckDependent {

std::vector<EquilRecord>
getEquil(const EclipseState& state)
{
    const auto& init = state.getInitConfig();

    if(!init.hasEquil()) {
        throw std::domain_error("Deck does not provide equilibration data.");
    }

    const auto& equil = init.getEquil();
    return { equil.begin(), equil.end() };
}

template<class GridView>
std::vector<int>
equilnum(const EclipseState& eclipseState,
         const GridView& gridview)
{
    std::vector<int> eqlnum(gridview.size(0), 0);

    if (eclipseState.fieldProps().has_int("EQLNUM")) {
        const auto& e = eclipseState.fieldProps().get_int("EQLNUM");
        std::transform(e.begin(), e.end(), eqlnum.begin(), [](int n){ return n - 1;});
    }
    OPM_BEGIN_PARALLEL_TRY_CATCH();
    const int num_regions = eclipseState.getTableManager().getEqldims().getNumEquilRegions();
    if ( std::any_of(eqlnum.begin(), eqlnum.end(), [num_regions](int n){return n >= num_regions;}) ) {
        throw std::runtime_error("Values larger than maximum Equil regions " + std::to_string(num_regions) + " provided in EQLNUM");
    }
    if ( std::any_of(eqlnum.begin(), eqlnum.end(), [](int n){return n < 0;}) ) {
        throw std::runtime_error("zero or negative values provided in EQLNUM");
    }
    OPM_END_PARALLEL_TRY_CATCH("Invalied EQLNUM numbers: ", gridview.comm());

    return eqlnum;
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class MaterialLawManager>
InitialStateComputer<FluidSystem,
                     Grid,
                     GridView,
                     ElementMapper,
                     CartesianIndexMapper>::
InitialStateComputer(MaterialLawManager& materialLawManager,
                     const EclipseState& eclipseState,
                     const Grid& grid,
                     const GridView& gridView,
                     const CartesianIndexMapper& cartMapper,
                     const double grav,
                     const int num_pressure_points,
                     const bool applySwatInit)
    : temperature_(grid.size(/*codim=*/0), eclipseState.getTableManager().rtemp()),
      saltConcentration_(grid.size(/*codim=*/0)),
      saltSaturation_(grid.size(/*codim=*/0)),
      pp_(FluidSystem::numPhases,
          std::vector<double>(grid.size(/*codim=*/0))),
      sat_(FluidSystem::numPhases,
           std::vector<double>(grid.size(/*codim=*/0))),
      rs_(grid.size(/*codim=*/0)),
      rv_(grid.size(/*codim=*/0)),
      rvw_(grid.size(/*codim=*/0)),
      cartesianIndexMapper_(cartMapper),
      num_pressure_points_(num_pressure_points)
{
    //Check for presence of kw SWATINIT
    if (applySwatInit) {
        if (eclipseState.fieldProps().has_double("SWATINIT")) {
            swatInit_ = eclipseState.fieldProps().get_double("SWATINIT");
        }
    }

    // Querry cell depth, cell top-bottom.
    // numerical aquifer cells might be specified with different depths.
    const auto& num_aquifers = eclipseState.aquifer().numericalAquifers();
    updateCellProps_(gridView, num_aquifers);

    // Get the equilibration records.
    const std::vector<EquilRecord> rec = getEquil(eclipseState);
    const auto& tables = eclipseState.getTableManager();
    // Create (inverse) region mapping.
    const RegionMapping<> eqlmap(equilnum(eclipseState, grid));
    const int invalidRegion = -1;
    regionPvtIdx_.resize(rec.size(), invalidRegion);
    setRegionPvtIdx(eclipseState, eqlmap);

    // Create Rs functions.
    rsFunc_.reserve(rec.size());
    if (FluidSystem::enableDissolvedGas()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rsFunc_.push_back(std::shared_ptr<Miscibility::RsVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = regionPvtIdx_[i];
            if (!rec[i].liveOilInitConstantRs()) {
                const TableContainer& rsvdTables = tables.getRsvdTables();
                const TableContainer& pbvdTables = tables.getPbvdTables();
                if (rsvdTables.size() > 0) {

                    const RsvdTable& rsvdTable = rsvdTables.getTable<RsvdTable>(i);
                    std::vector<double> depthColumn = rsvdTable.getColumn("DEPTH").vectorCopy();
                    std::vector<double> rsColumn = rsvdTable.getColumn("RS").vectorCopy();
                    rsFunc_.push_back(std::make_shared<Miscibility::RsVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rsColumn));
                } else if (pbvdTables.size() > 0) {
                    const PbvdTable& pbvdTable = pbvdTables.getTable<PbvdTable>(i);
                    std::vector<double> depthColumn = pbvdTable.getColumn("DEPTH").vectorCopy();
                    std::vector<double> pbubColumn = pbvdTable.getColumn("PBUB").vectorCopy();
                    rsFunc_.push_back(std::make_shared<Miscibility::PBVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, pbubColumn));

                } else {
                    throw std::runtime_error("Cannot initialise: RSVD or PBVD table not available.");
                }

            }
            else {
                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                    throw std::runtime_error("Cannot initialise: when no explicit RSVD table is given, \n"
                                             "datum depth must be at the gas-oil-contact. "
                                             "In EQUIL region "+std::to_string(i + 1)+"  (counting from 1), this does not hold.");
                }
                const double pContact = rec[i].datumDepthPressure();
                const double TContact = 273.15 + 20; // standard temperature for now
                rsFunc_.push_back(std::make_shared<Miscibility::RsSatAtContact<FluidSystem>>(pvtIdx, pContact, TContact));
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rsFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
        }
    }

    rvFunc_.reserve(rec.size());
    if (FluidSystem::enableVaporizedOil()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rvFunc_.push_back(std::shared_ptr<Miscibility::RvVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = regionPvtIdx_[i];
            if (!rec[i].wetGasInitConstantRv()) {
                const TableContainer& rvvdTables = tables.getRvvdTables();
                const TableContainer& pdvdTables = tables.getPdvdTables();

                if (rvvdTables.size() > 0) {
                    const RvvdTable& rvvdTable = rvvdTables.getTable<RvvdTable>(i);
                    std::vector<double> depthColumn = rvvdTable.getColumn("DEPTH").vectorCopy();
                    std::vector<double> rvColumn = rvvdTable.getColumn("RV").vectorCopy();
                    rvFunc_.push_back(std::make_shared<Miscibility::RvVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rvColumn));
                } else if (pdvdTables.size() > 0) {
                    const PdvdTable& pdvdTable = pdvdTables.getTable<PdvdTable>(i);
                    std::vector<double> depthColumn = pdvdTable.getColumn("DEPTH").vectorCopy();
                    std::vector<double> pdewColumn = pdvdTable.getColumn("PDEW").vectorCopy();
                    rvFunc_.push_back(std::make_shared<Miscibility::PDVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, pdewColumn));
                } else {
                    throw std::runtime_error("Cannot initialise: RVVD or PDCD table not available.");
                }
            }
            else {
                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                    throw std::runtime_error(
                              "Cannot initialise: when no explicit RVVD table is given, \n"
                              "datum depth must be at the gas-oil-contact. "
                              "In EQUIL region "+std::to_string(i + 1)+" (counting from 1), this does not hold.");
                }
                const double pContact = rec[i].datumDepthPressure() + rec[i].gasOilContactCapillaryPressure();
                const double TContact = 273.15 + 20; // standard temperature for now
                rvFunc_.push_back(std::make_shared<Miscibility::RvSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rvFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
        }
    }

 rvwFunc_.reserve(rec.size());
    if (FluidSystem::enableVaporizedWater()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rvwFunc_.push_back(std::shared_ptr<Miscibility::RvwVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = regionPvtIdx_[i];
            if (!rec[i].humidGasInitConstantRvw()) {
                const TableContainer& rvwvdTables = tables.getRvwvdTables();

                if (rvwvdTables.size() > 0) {
                    const RvwvdTable& rvwvdTable = rvwvdTables.getTable<RvwvdTable>(i);
                    std::vector<double> depthColumn = rvwvdTable.getColumn("DEPTH").vectorCopy();
                    std::vector<double> rvwvdColumn = rvwvdTable.getColumn("RVWVD").vectorCopy();
                    rvwFunc_.push_back(std::make_shared<Miscibility::RvwVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rvwvdColumn));
                } else {
                    throw std::runtime_error("Cannot initialise: RVWVD table not available.");
                }
            }
            else {
                const auto oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
                if (oilActive){
                    if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                        rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
                        const auto msg = "No explicit RVWVD table is given for EQUIL region " + std::to_string(i + 1) +". \n"
                                        "and datum depth is not at the gas-oil-contact. \n"
                                        "Rvw is set to 0.0 in all cells. \n";
                        OpmLog::warning(msg);
                    } else {
                        // pg = po + Pcgo = po + (pg - po)
                        // for gas-condensate with initial no oil zone: water-oil contact depth (OWC) equal gas-oil contact depth (GOC)
                        const double pContact = rec[i].datumDepthPressure() + rec[i].gasOilContactCapillaryPressure();
                        const double TContact = 273.15 + 20; // standard temperature for now
                        rvwFunc_.push_back(std::make_shared<Miscibility::RvwSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
                    }
                }
                else {
                     // two-phase gas-water sytem:  water-oil contact depth is taken equal to gas-water contact depth (GWC)
                     // and water-oil capillary pressure (Pcwo) is taken equal to gas-water capillary pressure (Pcgw) at GWC
                     if (rec[i].waterOilContactDepth() != rec[i].datumDepth()) {
                        rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
                        const auto msg = "No explicit RVWVD table is given for EQUIL region " + std::to_string(i + 1) +". \n"
                                        "and datum depth is not at the gas-water-contact. \n"
                                        "Rvw is set to 0.0 in all cells. \n";
                        OpmLog::warning(msg);
                    } else {
                        // pg = pw + Pcgw = pw + (pg - pw)
                        const double pContact = rec[i].datumDepthPressure() + rec[i].waterOilContactCapillaryPressure();
                        const double TContact = 273.15 + 20; // standard temperature for now
                        rvwFunc_.push_back(std::make_shared<Miscibility::RvwSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
                    }
                }
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
        }
    }


    // EXTRACT the initial temperature
    updateInitialTemperature_(eclipseState, eqlmap);

    // EXTRACT the initial salt concentration
    updateInitialSaltConcentration_(eclipseState, eqlmap);

    // EXTRACT the initial salt saturation
    updateInitialSaltSaturation_(eclipseState, eqlmap);

    // Compute pressures, saturations, rs and rv factors.
    const auto& comm = grid.comm();
    calcPressSatRsRv(eqlmap, rec, materialLawManager, comm, grav);

    // modify the pressure and saturation for numerical aquifer cells
    applyNumericalAquifers_(gridView, num_aquifers, eclipseState.runspec().co2Storage() || eclipseState.runspec().h2Storage());

    // Modify oil pressure in no-oil regions so that the pressures of present phases can
    // be recovered from the oil pressure and capillary relations.
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
updateInitialTemperature_(const EclipseState& eclState, const RMap& reg)
{
    const int numEquilReg = rsFunc_.size();
    tempVdTable_.resize(numEquilReg);
    const auto& tables = eclState.getTableManager();
    if (!tables.hasTables("RTEMPVD")) {
        std::vector<double> x = {0.0,1.0};
        std::vector<double> y = {tables.rtemp(),tables.rtemp()};
        for (auto& table : this->tempVdTable_) {
            table.setXYContainers(x, y);
        }
    } else {
        const TableContainer& tempvdTables = tables.getRtempvdTables();
        for (std::size_t i = 0; i < tempvdTables.size(); ++i) {
            const RtempvdTable& tempvdTable = tempvdTables.getTable<RtempvdTable>(i);
            tempVdTable_[i].setXYContainers(tempvdTable.getDepthColumn(), tempvdTable.getTemperatureColumn());
            const auto& cells = reg.cells(i);
            for (const auto& cell : cells) {
                const double depth = cellCenterDepth_[cell];
                this->temperature_[cell] = tempVdTable_[i].eval(depth, /*extrapolate=*/true);
            }
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
updateInitialSaltConcentration_(const EclipseState& eclState, const RMap& reg)
{
    const int numEquilReg = rsFunc_.size();
    saltVdTable_.resize(numEquilReg);
    const auto& tables = eclState.getTableManager();
    const TableContainer& saltvdTables = tables.getSaltvdTables();

    // If no saltvd table is given, we create a trivial table for the density calculations
    if (saltvdTables.empty()) {
        std::vector<double> x = {0.0,1.0};
        std::vector<double> y = {0.0,0.0};
        for (auto& table : this->saltVdTable_) {
            table.setXYContainers(x, y);
        }
    } else {
        for (std::size_t i = 0; i < saltvdTables.size(); ++i) {
            const SaltvdTable& saltvdTable = saltvdTables.getTable<SaltvdTable>(i);
            saltVdTable_[i].setXYContainers(saltvdTable.getDepthColumn(), saltvdTable.getSaltColumn());

            const auto& cells = reg.cells(i);
            for (const auto& cell : cells) {
                const double depth = cellCenterDepth_[cell];
                this->saltConcentration_[cell] = saltVdTable_[i].eval(depth, /*extrapolate=*/true);
            }
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
updateInitialSaltSaturation_(const EclipseState& eclState, const RMap& reg)
{
    const int numEquilReg = rsFunc_.size();
    saltpVdTable_.resize(numEquilReg);
    const auto& tables = eclState.getTableManager();
    const TableContainer& saltpvdTables = tables.getSaltpvdTables();

    for (std::size_t i = 0; i < saltpvdTables.size(); ++i) {
        const SaltpvdTable& saltpvdTable = saltpvdTables.getTable<SaltpvdTable>(i);
        saltpVdTable_[i].setXYContainers(saltpvdTable.getDepthColumn(), saltpvdTable.getSaltpColumn());

        const auto& cells = reg.cells(i);
        for (const auto& cell : cells) {
            const double depth = cellCenterDepth_[cell];
            this->saltSaturation_[cell] = saltpVdTable_[i].eval(depth, /*extrapolate=*/true);
        }
    }
}


template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
updateCellProps_(const GridView& gridView,
                 const NumericalAquifers& aquifer)
{
    ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
    int numElements = gridView.size(/*codim=*/0);
    cellCenterDepth_.resize(numElements);
    cellZSpan_.resize(numElements);
    cellZMinMax_.resize(numElements);

    auto elemIt = gridView.template begin</*codim=*/0>();
    const auto& elemEndIt = gridView.template end</*codim=*/0>();
    const auto num_aqu_cells = aquifer.allAquiferCells();
    for (; elemIt != elemEndIt; ++elemIt) {
        const Element& element = *elemIt;
        const unsigned int elemIdx = elemMapper.index(element);
        cellCenterDepth_[elemIdx] = Details::cellCenterDepth(element);
        const auto cartIx = cartesianIndexMapper_.cartesianIndex(elemIdx);
        cellZSpan_[elemIdx] = Details::cellZSpan(element);
        cellZMinMax_[elemIdx] = Details::cellZMinMax(element);
        if (!num_aqu_cells.empty()) {
            const auto search = num_aqu_cells.find(cartIx);
            if (search != num_aqu_cells.end()) {
                const auto* aqu_cell = num_aqu_cells.at(cartIx);
                const double depth_change_num_aqu = aqu_cell->depth - cellCenterDepth_[elemIdx];
                cellCenterDepth_[elemIdx] += depth_change_num_aqu;
                cellZSpan_[elemIdx].first += depth_change_num_aqu;
                cellZSpan_[elemIdx].second += depth_change_num_aqu;
                cellZMinMax_[elemIdx].first += depth_change_num_aqu;
                cellZMinMax_[elemIdx].second += depth_change_num_aqu;
            }
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
applyNumericalAquifers_(const GridView& gridView,
                        const NumericalAquifers& aquifer,
                        const bool co2store_or_h2store)
{
    const auto num_aqu_cells = aquifer.allAquiferCells();
    if (num_aqu_cells.empty()) return;

    // Check if water phase is active, or in the case of CO2STORE and H2STORE, water is modelled as oil phase
    bool oil_as_brine = co2store_or_h2store && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    const auto watPos =  oil_as_brine? FluidSystem::oilPhaseIdx : FluidSystem::waterPhaseIdx;
    if (!FluidSystem::phaseIsActive(watPos)){
        throw std::logic_error  { "Water phase has to be active for numerical aquifer case" };
    }

    ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
    auto elemIt = gridView.template begin</*codim=*/0>();
    const auto& elemEndIt = gridView.template end</*codim=*/0>();
    const auto oilPos = FluidSystem::oilPhaseIdx;
    const auto gasPos = FluidSystem::gasPhaseIdx;
    for (; elemIt != elemEndIt; ++elemIt) {
        const Element& element = *elemIt;
        const unsigned int elemIdx = elemMapper.index(element);
        const auto cartIx = cartesianIndexMapper_.cartesianIndex(elemIdx);
        const auto search = num_aqu_cells.find(cartIx);
        if (search != num_aqu_cells.end()) {
            // numerical aquifer cells are filled with water initially
            this->sat_[watPos][elemIdx] = 1.;

            if (!co2store_or_h2store && FluidSystem::phaseIsActive(oilPos)) {
                this->sat_[oilPos][elemIdx] = 0.;
            }

            if (FluidSystem::phaseIsActive(gasPos)) {
                this->sat_[gasPos][elemIdx] = 0.;
            }
            const auto* aqu_cell = num_aqu_cells.at(cartIx);
            const auto msg = fmt::format("FOR AQUIFER CELL AT ({}, {}, {}) OF NUMERICAL "
                                         "AQUIFER {}, WATER SATURATION IS SET TO BE UNITY",
                                         aqu_cell->I+1, aqu_cell->J+1, aqu_cell->K+1, aqu_cell->aquifer_id);
            OpmLog::info(msg);

            // if pressure is specified for numerical aquifers, we use these pressure values
            // for numerical aquifer cells
            if (aqu_cell->init_pressure) {
                const double pres = *(aqu_cell->init_pressure);
                this->pp_[watPos][elemIdx] = pres;
                if (FluidSystem::phaseIsActive(gasPos)) {
                    this->pp_[gasPos][elemIdx] = pres;
                }
                if (FluidSystem::phaseIsActive(oilPos)) {
                    this->pp_[oilPos][elemIdx] = pres;
                }
            }
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
setRegionPvtIdx(const EclipseState& eclState, const RMap& reg)
{
    const auto& pvtnumData = eclState.fieldProps().get_int("PVTNUM");

    for (const auto& r : reg.activeRegions()) {
        const auto& cells = reg.cells(r);
        regionPvtIdx_[r] = pvtnumData[*cells.begin()] - 1;
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap, class MaterialLawManager, class Comm>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
calcPressSatRsRv(const RMap& reg,
                 const std::vector<EquilRecord>& rec,
                 MaterialLawManager& materialLawManager,
                 const Comm& comm,
                 const double grav)
{
    using PhaseSat = Details::PhaseSaturations<
        MaterialLawManager, FluidSystem, EquilReg, typename RMap::CellId
    >;
    
    auto ptable = Details::PressureTable<FluidSystem, EquilReg>{ grav, this->num_pressure_points_ };
    auto psat   = PhaseSat { materialLawManager, this->swatInit_ };
    auto vspan  = std::array<double, 2>{};

    std::vector<int> regionIsEmpty(rec.size(), 0);
    for (std::size_t r = 0; r < rec.size(); ++r) {
        const auto& cells = reg.cells(r);

        Details::verticalExtent(cells, cellZMinMax_, comm, vspan);

        const auto acc = rec[r].initializationTargetAccuracy();
        if (acc > 0) {
            throw std::runtime_error {
                "Cannot initialise model: Positive item 9 is not supported "
                "in EQUIL keyword, record " + std::to_string(r + 1)
            };
        }

        if (cells.empty()) {
            regionIsEmpty[r] = 1;
            continue;
        }

        const auto eqreg = EquilReg {
            rec[r], this->rsFunc_[r], this->rvFunc_[r], this->rvwFunc_[r], this->tempVdTable_[r], this->saltVdTable_[r], this->regionPvtIdx_[r]
        };

        // Ensure gas/oil and oil/water contacts are within the span for the
        // phase pressure calculation.
        vspan[0] = std::min(vspan[0], std::min(eqreg.zgoc(), eqreg.zwoc()));
        vspan[1] = std::max(vspan[1], std::max(eqreg.zgoc(), eqreg.zwoc()));

        ptable.equilibrate(eqreg, vspan);

        if (acc == 0) {
            // Centre-point method
            this->equilibrateCellCentres(cells, eqreg, ptable, psat);
        }
        else if (acc < 0) {
            // Horizontal subdivision
            this->equilibrateHorizontal(cells, eqreg, -acc,
                                        ptable, psat);
        } else {
            // Horizontal subdivision with titled fault blocks
            // the simulator throw a few line above for the acc > 0 case
            // i.e. we should not reach here.
            assert(false);
        }
    }
    comm.min(regionIsEmpty.data(),regionIsEmpty.size());
    if (comm.rank() == 0) {
        for (std::size_t r = 0; r < rec.size(); ++r) {
            if (regionIsEmpty[r]) //region is empty on all partitions
                OpmLog::warning("Equilibration region " + std::to_string(r + 1)
                                 + " has no active cells");
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class EquilibrationMethod>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
cellLoop(const CellRange&      cells,
         EquilibrationMethod&& eqmethod)
{
    const auto oilPos = FluidSystem::oilPhaseIdx;
    const auto gasPos = FluidSystem::gasPhaseIdx;
    const auto watPos = FluidSystem::waterPhaseIdx;

    const auto oilActive = FluidSystem::phaseIsActive(oilPos);
    const auto gasActive = FluidSystem::phaseIsActive(gasPos);
    const auto watActive = FluidSystem::phaseIsActive(watPos);

    auto pressures   = Details::PhaseQuantityValue{};
    auto saturations = Details::PhaseQuantityValue{};
    auto Rs          = 0.0;
    auto Rv          = 0.0;
    auto Rvw         = 0.0;

    for (const auto& cell : cells) {
        eqmethod(cell, pressures, saturations, Rs, Rv, Rvw);

        if (oilActive) {
            this->pp_ [oilPos][cell] = pressures.oil;
            this->sat_[oilPos][cell] = saturations.oil;
        }

        if (gasActive) {
            this->pp_ [gasPos][cell] = pressures.gas;
            this->sat_[gasPos][cell] = saturations.gas;
        }

        if (watActive) {
            this->pp_ [watPos][cell] = pressures.water;
            this->sat_[watPos][cell] = saturations.water;
        }

        if (oilActive && gasActive) {
            this->rs_[cell] = Rs;
            this->rv_[cell] = Rv;
        }

        if (watActive && gasActive) {
            this->rvw_[cell] = Rvw;
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
equilibrateCellCentres(const CellRange&         cells,
                       const EquilReg&          eqreg,
                       const PressTable&        ptable,
                       PhaseSat&                psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;
    this->cellLoop(cells, [this, &eqreg,  &ptable, &psat]
        (const CellID                 cell,
         Details::PhaseQuantityValue& pressures,
         Details::PhaseQuantityValue& saturations,
         double&                      Rs,
         double&                      Rv,
         double&                      Rvw) -> void
    {
        const auto pos = CellPos {
            cell, cellCenterDepth_[cell]
        };

        saturations = psat.deriveSaturations(pos, eqreg, ptable);
        pressures   = psat.correctedPhasePressures();

        const auto temp = this->temperature_[cell];

        Rs = eqreg.dissolutionCalculator()
            (pos.depth, pressures.oil, temp, saturations.gas);

        Rv = eqreg.evaporationCalculator()
            (pos.depth, pressures.gas, temp, saturations.oil);

        Rvw = eqreg.waterEvaporationCalculator()
            (pos.depth, pressures.gas, temp, saturations.water);
    });
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
equilibrateHorizontal(const CellRange&  cells,
                      const EquilReg&   eqreg,
                      const int         acc,
                      const PressTable& ptable,
                      PhaseSat&         psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;

    this->cellLoop(cells, [this, acc, &eqreg, &ptable, &psat]
        (const CellID                 cell,
         Details::PhaseQuantityValue& pressures,
         Details::PhaseQuantityValue& saturations,
         double&                      Rs,
         double&                      Rv,
         double&                      Rvw) -> void
    {
        pressures  .reset();
        saturations.reset();

        auto totfrac = 0.0;
        for (const auto& [depth, frac] : Details::horizontalSubdivision(cell, cellZSpan_[cell], acc)) {
            const auto pos = CellPos { cell, depth };

            saturations.axpy(psat.deriveSaturations(pos, eqreg, ptable), frac);
            pressures  .axpy(psat.correctedPhasePressures(), frac);

            totfrac += frac;
        }

        if (totfrac > 0.) {
            saturations /= totfrac;
            pressures /= totfrac;
        } else {
            // Fall back to centre point method for zero-thickness cells.
            const auto pos = CellPos {
                    cell, cellCenterDepth_[cell]
            };

            saturations = psat.deriveSaturations(pos, eqreg, ptable);
            pressures   = psat.correctedPhasePressures();
        }

        const auto temp = this->temperature_[cell];
        const auto cz   = cellCenterDepth_[cell];

        Rs = eqreg.dissolutionCalculator()
            (cz, pressures.oil, temp, saturations.gas);

        Rv = eqreg.evaporationCalculator()
            (cz, pressures.gas, temp, saturations.oil);

        Rvw = eqreg.waterEvaporationCalculator()
            (cz, pressures.gas, temp, saturations.water);
    });
}

}
} // namespace EQUIL
} // namespace Opm
#endif
