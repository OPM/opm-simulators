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
/**
 * \file
 *
 * \brief Routines that actually solve the ODEs that emerge from the hydrostatic
 *        equilibrium problem
 */
#ifndef EWOMS_INITSTATEEQUIL_HH
#define EWOMS_INITSTATEEQUIL_HH

#include "equilibrationhelpers.hh"
#include "opm/grid/utility/RegionMapping.hpp"

#include <opm/models/utils/propertysystem.hh>

#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RsvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RvvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PbvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PdvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltpvdTable.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <fmt/format.h>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <string>

namespace Opm {

/**
 * Types and routines that collectively implement a basic
 * ECLIPSE-style equilibration-based initialisation scheme.
 *
 * This namespace is intentionally nested to avoid name clashes
 * with other parts of OPM.
 */
namespace EQUIL {
namespace Details {
template <class RHS>
class RK4IVP {
public:
    RK4IVP(const RHS& f,
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

    double
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

private:
    int N_;
    std::array<double,2> span_;
    std::vector<double>  y_;
    std::vector<double>  f_;

    double
    stepsize() const { return (span_[1] - span_[0]) / N_; }
};

namespace PhasePressODE {
template <class FluidSystem>
class Water
{
using TabulatedFunction = Tabulated1DFunction<double>;
public:
    Water(const double temp,
          const TabulatedFunction& saltVdTable,
          const int pvtRegionIdx,
          const double normGrav)
        : temp_(temp)
        , saltVdTable_(saltVdTable)
        , pvtRegionIdx_(pvtRegionIdx)
        , g_(normGrav)
    {}

    double
    operator()(const double depth,
               const double press) const
    {
        return this->density(depth, press) * g_;
    }

private:
    const double temp_;
    const TabulatedFunction& saltVdTable_;
    const int pvtRegionIdx_;
    const double g_;

    double
    density(const double depth,
            const double press) const
    {
        // The initializing algorithm can give depths outside the range due to numerical noise i.e. we extrapolate
        double saltConcentration = saltVdTable_.eval(depth, /*extrapolate=*/true);
        double rho = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp_, press, saltConcentration);
        rho *= FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx_);
        return rho;
    }
};

template <class FluidSystem, class RS>
class Oil
{
public:
    Oil(const double temp,
        const RS& rs,
        const int pvtRegionIdx,
        const double normGrav)
        : temp_(temp)
        , rs_(rs)
        , pvtRegionIdx_(pvtRegionIdx)
        , g_(normGrav)
    {}

    double
    operator()(const double depth,
               const double press) const
    {
        return this->density(depth, press) * g_;
    }

private:
    const double temp_;
    const RS& rs_;
    const int pvtRegionIdx_;
    const double g_;

    double
    density(const double depth,
            const double press) const
    {
        double rs = 0.0;
        if(FluidSystem::enableDissolvedGas())
            rs = rs_(depth, press, temp_);

        double bOil = 0.0;
        if (rs >= FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp_, press)) {
            bOil = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
        }
        else {
            bOil = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp_, press, rs);
        }
        double rho = bOil * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx_);
        if (FluidSystem::enableDissolvedGas()) {
            rho += rs * bOil * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
        }

        return rho;
    }
};

template <class FluidSystem, class RV>
class Gas
{
public:
    Gas(const double temp,
        const RV& rv,
        const int pvtRegionIdx,
        const double normGrav)
        : temp_(temp)
        , rv_(rv)
        , pvtRegionIdx_(pvtRegionIdx)
        , g_(normGrav)
    {}

    double
    operator()(const double depth,
               const double press) const
    {
        return this->density(depth, press) * g_;
    }

private:
    const double temp_;
    const RV& rv_;
    const int pvtRegionIdx_;
    const double g_;

    double
    density(const double depth,
            const double press) const
    {
        double rv = 0.0;
        if (FluidSystem::enableVaporizedOil())
            rv = rv_(depth, press, temp_);

        double bGas = 0.0;
        if (rv >= FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp_, press)) {
            bGas = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
        }
        else {
            bGas = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp_, press, rv);
        }
        double rho = bGas * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx_);
        if (FluidSystem::enableVaporizedOil()) {
            rho += rv * bGas * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx_);
        }

        return rho;
    }
};
} // namespace PhasePressODE

template <class FluidSystem, class Region>
class PressureTable
{
public:
    using VSpan = std::array<double, 2>;

    /// Constructor
    ///
    /// \param[in] gravity Norm of gravity vector (acceleration strength due
    ///    to gravity).  Normally the standardised value at Tellus equator
    ///    (9.80665 m/s^2).
    ///
    /// \param[in] samplePoints Number of equally spaced depth sample points
    ///    in each internal phase pressure table.
    explicit PressureTable(const double gravity,
                           const int    samplePoints = 2000)
        : gravity_(gravity)
        , nsample_(samplePoints)
    {}

    /// Copy constructor
    ///
    /// \param[in] rhs Source object for copy initialization.
    PressureTable(const PressureTable& rhs)
        : gravity_(rhs.gravity)
        , nsample_(rhs.nsample_)
    {
        this->copyInPointers(rhs);
    }

    /// Move constructor
    ///
    /// \param[in,out] rhs Source object for move initialization.  On output,
    ///    left in a moved-from ("valid but unspecified") state.  Internal
    ///    pointers in \p rhs are null (\c unique_ptr guarantee).
    PressureTable(PressureTable&& rhs)
        : gravity_(rhs.gravity_)
        , nsample_(rhs.nsample_)
        , oil_    (std::move(rhs.oil_))
        , gas_    (std::move(rhs.gas_))
        , wat_    (std::move(rhs.wat_))
    {}

    /// Assignment operator
    ///
    /// \param[in] rhs Source object.
    ///
    /// \return \code *this \endcode.
    PressureTable& operator=(const PressureTable& rhs)
    {
        this->gravity_ = rhs.gravity_;
        this->nsample_ = rhs.nsample_;
        this->copyInPointers(rhs);

        return *this;
    }

    /// Move-assignment operator
    ///
    /// \param[in] rhs Source object.  On output, left in a moved-from ("valid
    ///    but unspecified") state.  Internal pointers in \p rhs are null (\c
    ///    unique_ptr guarantee).
    ///
    /// \return \code *this \endcode.
    PressureTable& operator=(PressureTable&& rhs)
    {
        this->gravity_ = rhs.gravity_;
        this->nsample_ = rhs.nsample_;

        this->oil_ = std::move(rhs.oil_);
        this->gas_ = std::move(rhs.gas_);
        this->wat_ = std::move(rhs.wat_);

        return *this;
    }

    void equilibrate(const Region& reg,
                     const VSpan&  span)
    {
        // One of the PressureTable::equil_*() member functions.
        auto equil = this->selectEquilibrationStrategy(reg);

        (this->*equil)(reg, span);
    }

    /// Predicate for whether or not oil is an active phase
    bool oilActive() const
    {
        return FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    }

    /// Predicate for whether or not gas is an active phase
    bool gasActive() const
    {
        return FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    }

    /// Predicate for whether or not water is an active phase
    bool waterActive() const
    {
        return FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    }

    /// Evaluate oil phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Oil phase pressure at specified depth.
    double oil(const double depth) const
    {
        this->checkPtr(this->oil_.get(), "OIL");

        return this->oil_->value(depth);
    }

    /// Evaluate gas phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Gas phase pressure at specified depth.
    double gas(const double depth) const
    {
        this->checkPtr(this->gas_.get(), "GAS");

        return this->gas_->value(depth);
    }

    /// Evaluate water phase pressure at specified depth.
    ///
    /// \param[in] depth Depth of evaluation point.  Should generally be
    ///    within the \c span from the previous call to \code equilibrate()
    ///    \endcode.
    ///
    /// \return Water phase pressure at specified depth.
    double water(const double depth) const
    {
        this->checkPtr(this->wat_.get(), "WATER");

        return this->wat_->value(depth);
    }

private:
    template <class ODE>
    class PressureFunction
    {
    public:
        struct InitCond {
            double depth;
            double pressure;
        };

        explicit PressureFunction(const ODE&      ode,
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

        PressureFunction(const PressureFunction& rhs)
            : initial_(rhs.initial_)
        {
            this->value_[Direction::Up] =
                std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

            this->value_[Direction::Down] =
                std::make_unique<Distribution>(*rhs.value_[Direction::Down]);
        }

        PressureFunction(PressureFunction&& rhs) = default;

        PressureFunction& operator=(const PressureFunction& rhs)
        {
            this->initial_ = rhs.initial_;

            this->value_[Direction::Up] =
                std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

            this->value_[Direction::Down] =
                std::make_unique<Distribution>(*rhs.value_[Direction::Down]);

            return *this;
        }

        PressureFunction& operator=(PressureFunction&& rhs)
        {
            this->initial_ = rhs.initial_;
            this->value_   = std::move(rhs.value_);

            return *this;
        }

        double value(const double depth) const
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

    private:
        enum Direction : std::size_t { Up, Down, NumDir };

        using Distribution = Details::RK4IVP<ODE>;
        using DistrPtr = std::unique_ptr<Distribution>;

        InitCond initial_;
        std::array<DistrPtr, Direction::NumDir> value_;
    };

    using OilPressODE = PhasePressODE::Oil<
        FluidSystem, typename Region::CalcDissolution
    >;

    using GasPressODE = PhasePressODE::Gas<
        FluidSystem, typename Region::CalcEvaporation
    >;

    using WatPressODE = PhasePressODE::Water<FluidSystem>;

    using OPress = PressureFunction<OilPressODE>;
    using GPress = PressureFunction<GasPressODE>;
    using WPress = PressureFunction<WatPressODE>;

    using Strategy = void (PressureTable::*)
        (const Region&, const VSpan&);

    double gravity_;
    int    nsample_;
    double temperature_{ 273.15 + 20 };

    std::unique_ptr<OPress> oil_{};
    std::unique_ptr<GPress> gas_{};
    std::unique_ptr<WPress> wat_{};

    template <typename PressFunc>
    void checkPtr(const PressFunc*   phasePress,
                  const std::string& phaseName) const
    {
        if (phasePress != nullptr) { return; }

        throw std::invalid_argument {
            "Phase pressure function for \"" + phaseName
            + "\" most not be null"
        };
    }

    Strategy selectEquilibrationStrategy(const Region& reg) const
    {
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

    void copyInPointers(const PressureTable& rhs)
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

    void equil_WOG(const Region& reg, const VSpan& span);
    void equil_GOW(const Region& reg, const VSpan& span);
    void equil_OWG(const Region& reg, const VSpan& span);

    void makeOilPressure(const typename OPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);

    void makeGasPressure(const typename GPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);

    void makeWatPressure(const typename WPress::InitCond& ic,
                         const Region&                    reg,
                         const VSpan&                     span);
};

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
    // Datum depth in gas zone.  Calculate phase pressure for gas first,
    // followed by oil and water if applicable.

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
        this->temperature_, reg.dissolutionCalculator(),
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
        this->temperature_, reg.evaporationCalculator(),
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
        this->temperature_, reg.saltVdTable(), reg.pvtIdx(), this->gravity_
    };

    this->wat_ = std::make_unique<WPress>(drho, ic, this->nsample_, span);
}

// ===========================================================================

/// Simple set of per-phase (named by primary component) quantities.
struct PhaseQuantityValue {
    double oil{0.0};
    double gas{0.0};
    double water{0.0};

    PhaseQuantityValue& axpy(const PhaseQuantityValue& rhs, const double a)
    {
        this->oil   += a * rhs.oil;
        this->gas   += a * rhs.gas;
        this->water += a * rhs.water;

        return *this;
    }

    PhaseQuantityValue& operator/=(const double x)
    {
        this->oil   /= x;
        this->gas   /= x;
        this->water /= x;

        return *this;
    }

    void reset()
    {
        this->oil = this->gas = this->water = 0.0;
    }
};

/// Calculator for phase saturations
///
/// Computes saturation values at arbitrary depths.
///
/// \tparam MaterialLawManager Container for material laws.  Typically a
///    specialization of the \code Opm::EclMaterialLawManager<> \endcode
///    template.
///
/// \tparam FluidSystem An OPM fluid system type.  Typically a
///    specialization of the \code Opm::BlackOilFluidSystem<> \endcode
///    template.
///
/// \tparam Region Representation of an equilibration region.  Typically
///    \code Opm::EQUIL::EquilReg \endcode from the equilibrationhelpers.
///
/// \tparam CellID Representation an equilibration region's cell IDs.
///    Typically \code std::size_t \endcode.
template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
class PhaseSaturations
{
public:
    /// Evaluation point within a model geometry.
    ///
    /// Associates a particular depth to specific cell.
    struct Position {
        CellID cell;
        double depth;
    };

    /// Convenience type alias
    using PTable = PressureTable<FluidSystem, Region>;

    /// Constructor
    ///
    /// \param[in,out] matLawMgr Read/write reference to a material law
    ///    container.  Mutated by member functions.
    ///
    /// \param[in] swatInit Initial water saturation array (from SWATINIT
    ///    data).  Empty if SWATINIT is not used in this simulation model.
    explicit PhaseSaturations(MaterialLawManager&        matLawMgr,
                              const std::vector<double>& swatInit)
        : matLawMgr_(matLawMgr)
        , swatInit_ (swatInit)
    {}

    /// Copy constructor.
    ///
    /// \param[in] rhs Source object.
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

    /// Disabled assignment operator.
    PhaseSaturations& operator=(const PhaseSaturations&) = delete;

    /// Disabled move-assignment operator.
    PhaseSaturations& operator=(PhaseSaturations&&) = delete;

    /// Calculate phase saturations at particular point of the simulation
    /// model geometry.
    ///
    /// \param[in] x Specific geometric point (depth within a specific cell).
    ///
    /// \param[in] reg Equilibration information for a single equilibration
    ///    region; notably contact depths.
    ///
    /// \param[in] ptable Previously equilibrated phase pressure table
    ///    pertaining to the equilibration region \p reg.
    ///
    /// \return Set of phase saturation values defined at particular point.
    const PhaseQuantityValue&
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

    /// Retrieve saturation-corrected phase pressures
    ///
    /// Values associated with evaluation point of previous call to \code
    /// deriveSaturations() \endcode.
    const PhaseQuantityValue& correctedPhasePressures() const
    {
        return this->press_;
    }

private:
    /// Convenience amalgamation of the deriveSaturations() input state.
    /// These values are almost always used in concert.
    struct EvaluationPoint {
        const Position* position{nullptr};
        const Region*   region  {nullptr};
        const PTable*   ptable  {nullptr};
    };

    /// Simplified fluid state object that contains only the pieces of
    /// information needed to calculate the capillary pressure values from
    /// the current set of material laws.
    using FluidState = ::Opm::
        SimpleModularFluidState<double, /*numPhases=*/3, /*numComponents=*/3,
                                FluidSystem,
                                /*storePressure=*/false,
                                /*storeTemperature=*/false,
                                /*storeComposition=*/false,
                                /*storeFugacity=*/false,
                                /*storeSaturation=*/true,
                                /*storeDensity=*/false,
                                /*storeViscosity=*/false,
                                /*storeEnthalpy=*/false>;

    /// Convenience type alias.
    using MaterialLaw = typename MaterialLawManager::MaterialLaw;

    /// Fluid system's representation of phase indices.
    using PhaseIdx = std::remove_cv_t<
        std::remove_reference_t<decltype(FluidSystem::oilPhaseIdx)>
    >;

    /// Read/write reference to client's material law container.
    MaterialLawManager& matLawMgr_;

    /// Client's SWATINIT data.
    const std::vector<double>& swatInit_;

    /// Evaluated phase saturations.
    PhaseQuantityValue sat_;

    /// Saturation-corrected phase pressure values.
    PhaseQuantityValue press_;

    /// Current evaluation point.
    EvaluationPoint evalPt_;

    /// Capillary pressure fluid state.
    FluidState fluidState_;

    /// Evaluated capillary pressures from current set of material laws.
    std::array<double, FluidSystem::numPhases> matLawCapPress_;

    /// Capture the input evaluation point information in internal state.
    ///
    /// \param[in] x Specific geometric point (depth within a specific cell).
    ///
    /// \param[in] reg Equilibration information for a single equilibration
    ///    region; notably contact depths.
    ///
    /// \param[in] ptable Previously equilibrated phase pressure table
    ///    pertaining to the equilibration region \p reg.
    void setEvaluationPoint(const Position& x,
                            const Region&   reg,
                            const PTable&   ptable)
    {
        this->evalPt_.position = &x;
        this->evalPt_.region   = &reg;
        this->evalPt_.ptable   = &ptable;
    }

    /// Initialize phase saturation and phase pressure values.
    ///
    /// Looks up phase pressure values from the input pressure table.
    void initializePhaseQuantities()
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

    /// Derive phase saturation for oil.
    ///
    /// Calculated as 1 - Sw - Sg.
    void deriveOilSat();

    /// Derive phase saturation for gas.
    ///
    /// Inverts capillary pressure curve if non-constant or uses a simple
    /// depth consideration with respect to G/O contact depth otherwise.
    void deriveGasSat();

    /// Derive phase saturation for water.
    ///
    /// Uses input data if simulation model is defined in terms of SWATINIT.
    /// Otherwise, inverts capillary pressure curve if non-constant or uses
    /// a simple depth consideration with respect to the O/W contact depth
    /// if capillary pressure curve is constant within the current cell.
    void deriveWaterSat();

    /// Correct phase saturation and pressure values to account for
    /// overlapping transition zones between G/O and O/W systems.
    void fixUnphysicalTransition();

    /// Re-adjust phase pressure values to account for phase saturations
    /// outside permissible ranges.
    void accountForScaledSaturations();

    // --------------------------------------------------------------------
    // Note: Function 'applySwatInit' is non-const because the overload set
    // needs to mutate the 'matLawMgr_'.
    // --------------------------------------------------------------------

    /// Derive water saturation from SWATINIT data.
    ///
    /// Uses SWATINIT array data from current cell directly.  Also updates
    /// the material law container's internal notion of the maximum
    /// attainable O/W capillary pressure value.
    ///
    /// \param[in] pcow O/W capillary pressure value (Po - Pw).
    ///
    /// \return Water saturation value.
    double applySwatInit(const double pcow);

    /// Derive water saturation from SWATINIT data.
    ///
    /// Uses explicitly passed-in saturation value.  Also updates the
    /// material law container's internal notion of the maximum attainable
    /// O/W capillary pressure value.
    ///
    /// \param[in] pc x/W capillary pressure value (Px - Pw; x in {O, G}).
    ///
    /// \param[in] sw Water saturation value.
    ///
    /// \return Water saturation value.  Input value, possibly mollified by
    ///    current set of material laws.
    double applySwatInit(const double pc, const double sw);

    /// Invoke material law container's capillary pressure calculator on
    /// current fluid state.
    void computeMaterialLawCapPress();

    /// Extract gas/oil capillary pressure value (Pg - Po) from current
    /// fluid state.
    double materialLawCapPressGasOil() const;

    /// Extract oil/water capillary pressure value (Po - Pw) from current
    /// fluid state.
    double materialLawCapPressOilWater() const;

    /// Predicate for whether specific phase has constant capillary pressure
    /// curve in current cell.
    ///
    /// \param[in] phaseIdx Phase.  Typically gas or water.
    ///
    /// \return Whether or not \p phaseIdx has constant capillary pressure
    /// curve in current cell.
    bool isConstCapPress(const PhaseIdx phaseIdx) const;

    /// Predicate for whether or not the G/O and O/W transition zones
    /// overlap in the current cell.
    ///
    /// This is the case when inverting the capillary pressure curves
    /// produces a negative oil saturation--i.e., when Sg + Sw > 1.
    bool isOverlappingTransition() const;

    /// Derive phase saturation value from simple depth consideration.
    ///
    /// Assumes that the pertinent capillary pressure curve is constant
    /// (typically zero) in the current cell--i.e., that there is a sharp
    /// interface between the two phases.
    ///
    /// \param[in] contactdepth Depth of relevant phase separation contact.
    ///
    /// \param[in] Position of phase in three-phase enumeration.  Typically
    ///    \code gasPos() \endcode or \code waterPos() \endcode.
    ///
    /// \param[in] isincr Whether the capillary pressure curve is normally
    ///    increasing as a function of phase saturation (e.g., Pcgo(Sg) = Pg
    ///    - Po) or if the curve is normally decreasing as a function of
    ///    increasing phase saturation (e.g., Pcow(Sw) = Po - Pw).  True for
    ///    capillary pressure functions that are normally increasing as a
    ///    function of phase saturation.
    ///
    /// \return Phase saturation.
    double fromDepthTable(const double   contactdepth,
                          const PhaseIdx phasePos,
                          const bool     isincr) const;

    /// Derive phase saturation by inverting non-constant capillary pressure
    /// curve.
    ///
    /// \param[in] pc Target capillary pressure value.
    ///
    /// \param[in] Position of phase in three-phase enumeration.  Typically
    ///    \code gasPos() \endcode or \code waterPos() \endcode.
    ///
    /// \param[in] isincr Whether the capillary pressure curve is normally
    ///    increasing as a function of phase saturation (e.g., Pcgo(Sg) = Pg
    ///    - Po) or if the curve is normally decreasing as a function of
    ///    increasing phase saturation (e.g., Pcow(Sw) = Po - Pw).  True for
    ///    capillary pressure functions that are normally increasing as a
    ///    function of phase saturation.
    ///
    /// \return Phase saturation at which capillary pressure attains target
    ///    value.
    double invertCapPress(const double   pc,
                          const PhaseIdx phasePos,
                          const bool     isincr) const;

    /// Position of oil in fluid system's three-phase enumeration.
    PhaseIdx oilPos() const
    {
        return FluidSystem::oilPhaseIdx;
    }

    /// Position of gas in fluid system's three-phase enumeration.
    PhaseIdx gasPos() const
    {
        return FluidSystem::gasPhaseIdx;
    }

    /// Position of water in fluid system's three-phase enumeration.
    PhaseIdx waterPos() const
    {
        return FluidSystem::waterPhaseIdx;
    }
};

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

    if (this->isConstCapPress(this->gasPos())) {
        // Sharp interface between phases.  Can derive phase saturation
        // directly from knowing where 'depth' of evaluation point is
        // relative to depth of O/G contact.
        sg = this->fromDepthTable(this->evalPt_.region->zgoc(),
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
        const auto pcgo = this->press_.gas - this->press_.oil;

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

        sw = this->swatInit_.empty()
            ? this->invertCapPress(pcow, this->waterPos(), isIncr)
            : this->applySwatInit(pcow);
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
        sw = this->applySwatInit(pcgw, sw);
    }

    sw = satFromSumOfPcs<FluidSystem, MaterialLaw>
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

    const auto& scaledDrainageInfo = this->matLawMgr_
        .oilWaterScaledEpsInfoDrainage(this->evalPt_.position->cell);

    const auto sg = this->sat_.gas;
    const auto sw = this->sat_.water;

    {
        auto so = 1.0;

        if (watActive) {
            const auto swu = scaledDrainageInfo.Swu;
            so -= swu;

            this->fluidState_.setSaturation(this->waterPos(), swu);
        }

        if (gasActive) {
            const auto sgu = scaledDrainageInfo.Sgu;
            so -= sgu;

            this->fluidState_.setSaturation(this->gasPos(), sgu);
        }

        this->fluidState_.setSaturation(this->oilPos(), so);
    }

    const auto thresholdSat = 1.0e-6;
    if (watActive && ((sw + thresholdSat) > scaledDrainageInfo.Swu)) {
        // Water saturation exceeds maximum possible value.  Reset oil phase
        // pressure to that which corresponds to maximum possible water
        // saturation value.
        this->fluidState_.setSaturation(this->waterPos(), scaledDrainageInfo.Swu);
        this->computeMaterialLawCapPress();

        // Pcow = Po - Pw => Po = Pw + Pcow
        this->press_.oil = this->press_.water + this->materialLawCapPressOilWater();
    }
    else if (gasActive && ((sg + thresholdSat) > scaledDrainageInfo.Sgu)) {
        // Gas saturation exceeds maximum possible value.  Reset oil phase
        // pressure to that which corresponds to maximum possible gas
        // saturation value.
        this->fluidState_.setSaturation(this->gasPos(), scaledDrainageInfo.Sgu);
        this->computeMaterialLawCapPress();

        // Pcgo = Pg - Po => Po = Pg - Pcgo
        this->press_.oil = this->press_.gas - this->materialLawCapPressGasOil();
    }

    if (gasActive && ((sg - thresholdSat) < scaledDrainageInfo.Sgl)) {
        // Gas saturation less than minimum possible value in cell.  Reset
        // gas phase pressure to that which corresponds to minimum possible
        // gas saturation.
        this->fluidState_.setSaturation(this->gasPos(), scaledDrainageInfo.Sgl);
        this->computeMaterialLawCapPress();

        // Pcgo = Pg - Po => Pg = Po + Pcgo
        this->press_.gas = this->press_.oil + this->materialLawCapPressGasOil();
    }

    if (watActive && ((sw - thresholdSat) < scaledDrainageInfo.Swl)) {
        // Water saturation less than minimum possible value in cell.  Reset
        // water phase pressure to that which corresponds to minimum
        // possible water saturation value.
        this->fluidState_.setSaturation(this->waterPos(), scaledDrainageInfo.Swl);
        this->computeMaterialLawCapPress();

        // Pcwo = Po - Pw => Pw = Po - Pcow
        this->press_.water = this->press_.oil - this->materialLawCapPressOilWater();
    }
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
applySwatInit(const double pcow)
{
    return this->applySwatInit(pcow, this->swatInit_[this->evalPt_.position->cell]);
}

template <class MaterialLawManager, class FluidSystem, class Region, typename CellID>
double PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
applySwatInit(const double pcow, const double sw)
{
    return this->matLawMgr_
        .applySwatinit(this->evalPt_.position->cell, pcow, sw);
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
bool PhaseSaturations<MaterialLawManager, FluidSystem, Region, CellID>::
isConstCapPress(const PhaseIdx phaseIdx) const
{
    return isConstPc<FluidSystem, MaterialLaw>
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
    return satFromDepth<FluidSystem, MaterialLaw>
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
    return satFromPc<FluidSystem, MaterialLaw>
        (this->matLawMgr_, static_cast<int>(phasePos),
         this->evalPt_.position->cell, pc, isincr);
}

// ===========================================================================

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

inline
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

} // namespace Details

namespace DeckDependent {
inline std::vector<EquilRecord>
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

    return eqlnum;
}

template<class TypeTag>
class InitialStateComputer
{
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;

public:
    template<class MaterialLawManager>
    InitialStateComputer(MaterialLawManager& materialLawManager,
                         const EclipseState& eclipseState,
                         const GridView& gridView,
                         const CartesianIndexMapper& cartMapper,
                         const double grav = unit::gravity,
                         const bool applySwatInit = true)
        : temperature_(gridView.size(/*codim=*/0)),
          saltConcentration_(gridView.size(/*codim=*/0)),
          saltSaturation_(gridView.size(/*codim=*/0)),
          pp_(FluidSystem::numPhases,
          std::vector<double>(gridView.size(/*codim=*/0))),
          sat_(FluidSystem::numPhases,
          std::vector<double>(gridView.size(/*codim=*/0))),
          rs_(gridView.size(/*codim=*/0)),
          rv_(gridView.size(/*codim=*/0)),
          cartesianIndexMapper_(cartMapper)
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
        const RegionMapping<> eqlmap(equilnum(eclipseState, gridView));
        const int invalidRegion = -1;
        regionPvtIdx_.resize(rec.size(), invalidRegion);
        setRegionPvtIdx(eclipseState, eqlmap);

        // Create Rs functions.
        rsFunc_.reserve(rec.size());
        if (FluidSystem::enableDissolvedGas()) {
            for (size_t i = 0; i < rec.size(); ++i) {
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
            for (size_t i = 0; i < rec.size(); ++i) {
                rsFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
            }
        }

        rvFunc_.reserve(rec.size());
        if (FluidSystem::enableVaporizedOil()) {
            for (size_t i = 0; i < rec.size(); ++i) {
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
            for (size_t i = 0; i < rec.size(); ++i) {
                rvFunc_.push_back(std::make_shared<Miscibility::NoMixing>());
            }
        }

        // EXTRACT the initial temperature
        updateInitialTemperature_(eclipseState);

        // EXTRACT the initial salt concentration
        updateInitialSaltConcentration_(eclipseState, eqlmap);

        // EXTRACT the initial salt saturation
        updateInitialSaltSaturation_(eclipseState, eqlmap);

        // Compute pressures, saturations, rs and rv factors.
        const auto& comm = gridView.comm();
        calcPressSatRsRv(eqlmap, rec, materialLawManager, comm, grav);

        // modify the pressure and saturation for numerical aquifer cells
        applyNumericalAquifers_(gridView, num_aquifers, eclipseState.runspec().co2Storage());

        // Modify oil pressure in no-oil regions so that the pressures of present phases can
        // be recovered from the oil pressure and capillary relations.
    }

    typedef std::vector<double> Vec;
    typedef std::vector<Vec>    PVec; // One per phase.

    const Vec& temperature() const { return temperature_; }
    const Vec& saltConcentration() const { return saltConcentration_; }
    const Vec& saltSaturation() const { return saltSaturation_; }
    const PVec& press() const { return pp_; }
    const PVec& saturation() const { return sat_; }
    const Vec& rs() const { return rs_; }
    const Vec& rv() const { return rv_; }
    const Vec& rvw() const { return rvw_; }

private:

    void updateInitialTemperature_(const EclipseState& eclState)
    {
        this->temperature_ = eclState.fieldProps().get_double("TEMPI");
    }

    template <class RMap>
    void updateInitialSaltConcentration_(const EclipseState& eclState, const RMap& reg)
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
            for (size_t i = 0; i < saltvdTables.size(); ++i) {
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
    template <class RMap>
    void updateInitialSaltSaturation_(const EclipseState& eclState, const RMap& reg)
    {
        const int numEquilReg = rsFunc_.size();
        saltpVdTable_.resize(numEquilReg);
        const auto& tables = eclState.getTableManager();
        const TableContainer& saltpvdTables = tables.getSaltpvdTables();

        for (size_t i = 0; i < saltpvdTables.size(); ++i) {
            const SaltpvdTable& saltpvdTable = saltpvdTables.getTable<SaltpvdTable>(i);
            saltpVdTable_[i].setXYContainers(saltpvdTable.getDepthColumn(), saltpvdTable.getSaltpColumn());

            const auto& cells = reg.cells(i);
            for (const auto& cell : cells) {
                const double depth = cellCenterDepth_[cell];
                this->saltSaturation_[cell] = saltpVdTable_[i].eval(depth, /*extrapolate=*/true);
            }
        }
    }


    std::vector< std::shared_ptr<Miscibility::RsFunction> > rsFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction> > rvFunc_;
    using TabulatedFunction = Tabulated1DFunction<double>;
    std::vector<TabulatedFunction> saltVdTable_;
    std::vector<TabulatedFunction> saltpVdTable_;
    std::vector<int> regionPvtIdx_;
    Vec temperature_;
    Vec saltConcentration_;
    Vec saltSaturation_;
    PVec pp_;
    PVec sat_;
    Vec rs_;
    Vec rv_;
    Vec rvw_;
    const CartesianIndexMapper& cartesianIndexMapper_;
    Vec swatInit_;
    Vec cellCenterDepth_;
    std::vector<std::pair<double,double>> cellZSpan_;
    std::vector<std::pair<double,double>> cellZMinMax_;

    void updateCellProps_(const GridView& gridView,
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
            if (!num_aqu_cells.empty()) {
                const auto search = num_aqu_cells.find(cartIx);
                if (search != num_aqu_cells.end()) {
                    const auto* aqu_cell = num_aqu_cells.at(cartIx);
                    cellCenterDepth_[elemIdx] = aqu_cell->depth;
                }
            }
            cellZSpan_[elemIdx] = Details::cellZSpan(element);
            cellZMinMax_[elemIdx] = Details::cellZMinMax(element);
        }
    }

    void applyNumericalAquifers_(const GridView& gridView,
                                 const NumericalAquifers& aquifer,
                                 const bool co2store)
    {
        const auto num_aqu_cells = aquifer.allAquiferCells();
        if (num_aqu_cells.empty()) return;

        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);
            const auto cartIx = cartesianIndexMapper_.cartesianIndex(elemIdx);
            const auto search = num_aqu_cells.find(cartIx);
            if (search != num_aqu_cells.end()) {
                // numerical aquifer cells are filled with water initially
                // for co2store the oilphase is used for brine
                const auto watPos = co2store? FluidSystem::oilPhaseIdx : FluidSystem::waterPhaseIdx;
                if (FluidSystem::phaseIsActive(watPos)) {
                    this->sat_[watPos][elemIdx] = 1.;
                } else {
                    throw std::logic_error  { "Water phase has to be active for numerical aquifer case" };
                }

                const auto oilPos = FluidSystem::oilPhaseIdx;
                if (!co2store && FluidSystem::phaseIsActive(oilPos)) {
                    this->sat_[oilPos][elemIdx] = 0.;
                }

                const auto gasPos = FluidSystem::gasPhaseIdx;
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
    template<class RMap>
    void setRegionPvtIdx(const EclipseState& eclState, const RMap& reg)
    {
        const auto& pvtnumData = eclState.fieldProps().get_int("PVTNUM");

        for (const auto& r : reg.activeRegions()) {
            const auto& cells = reg.cells(r);
            regionPvtIdx_[r] = pvtnumData[*cells.begin()] - 1;
        }
    }

    template <class RMap, class MaterialLawManager, class Comm>
    void calcPressSatRsRv(const RMap& reg,
                          const std::vector<EquilRecord>& rec,
                          MaterialLawManager& materialLawManager,
                          const Comm& comm,
                          const double grav)
    {
        using PhaseSat = Details::PhaseSaturations<
            MaterialLawManager, FluidSystem, EquilReg, typename RMap::CellId
        >;

        auto ptable = Details::PressureTable<FluidSystem, EquilReg>{ grav };
        auto psat   = PhaseSat { materialLawManager, this->swatInit_ };
        auto vspan  = std::array<double, 2>{};

        std::vector<int> regionIsEmpty(rec.size(), 0);
        for (size_t r = 0; r < rec.size(); ++r) {
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
                rec[r], this->rsFunc_[r], this->rvFunc_[r], this->saltVdTable_[r], this->regionPvtIdx_[r]
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
            for (size_t r = 0; r < rec.size(); ++r) {
                if (regionIsEmpty[r]) //region is empty on all partitions
                    OpmLog::warning("Equilibration region " + std::to_string(r + 1)
                                     + " has no active cells");
            }
        }
    }

    template <class CellRange, class EquilibrationMethod>
    void cellLoop(const CellRange&      cells,
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

        for (const auto& cell : cells) {
            eqmethod(cell, pressures, saturations, Rs, Rv);

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
        }
    }

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateCellCentres(const CellRange&         cells,
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
             double&                      Rv) -> void
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
        });
    }

    template <class CellRange, class PressTable, class PhaseSat>
    void equilibrateHorizontal(const CellRange&  cells,
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
             double&                      Rv) -> void
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

            saturations /= totfrac;
            pressures   /= totfrac;

            const auto temp = this->temperature_[cell];
            const auto cz   = cellCenterDepth_[cell];

            Rs = eqreg.dissolutionCalculator()
                (cz, pressures.oil, temp, saturations.gas);

            Rv = eqreg.evaporationCalculator()
                (cz, pressures.gas, temp, saturations.oil);
        });
    }
};
} // namespace DeckDependent
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
