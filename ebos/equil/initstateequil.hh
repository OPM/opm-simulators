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

#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RsvdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RvvdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PbvdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PdvdTable.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

BEGIN_PROPERTIES

NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(FluidSystem);

END_PROPERTIES

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
public:
    Water(const double temp,
          const int pvtRegionIdx,
          const double normGrav)
        : temp_(temp)
        , pvtRegionIdx_(pvtRegionIdx)
        , g_(normGrav)
    {}

    double
    operator()(const double /* depth */,
               const double press) const
    {
        return this->density(press) * g_;
    }

private:
    const double temp_;
    const int pvtRegionIdx_;
    const double g_;

    double
    density(const double press) const
    {
        double rho = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
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
        double rs = rs_(depth, press, temp_);
        double bOil = 0.0;
        if (!FluidSystem::enableDissolvedGas() || rs >= FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp_, press)) {
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
        double rv = rv_(depth, press, temp_);
        double bGas = 0.0;
        if (!FluidSystem::enableVaporizedOil() || rv >= FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp_, press)) {
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

    if (this->waterActive()) {
        // Pcow = Po - Pw => Pw = Po - Pcow
        const auto ic = typename WPress::InitCond {
            reg.zwoc(),
            this->oil(reg.zwoc()) - reg.pcowWoc()
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
        this->temperature_, reg.pvtIdx(), this->gravity_
    };

    this->wat_ = std::make_unique<WPress>(drho, ic, this->nsample_, span);
}

template <class FluidSystem,
          class Grid,
          class Region,
          class CellRange>
void equilibrateOWG(const Grid& grid,
                    const Region& reg,
                    const double grav,
                    const std::array<double,2>& span,
                    const CellRange& cells,
                    std::vector< std::vector<double> >& press)
{
    using PhaseIX = std::vector<std::vector<double>>::size_type;
    auto ptable   = PressureTable<FluidSystem, Region>{ grav };

    ptable.equilibrate(reg, span);

    const auto oilpos   = static_cast<PhaseIX>(FluidSystem::oilPhaseIdx);
    const auto gaspos   = static_cast<PhaseIX>(FluidSystem::gasPhaseIdx);
    const auto waterpos = static_cast<PhaseIX>(FluidSystem::waterPhaseIdx);

    auto ix = std::vector<double>::size_type{0};
    for (const auto& cell : cells) {
        const auto depth = Opm::UgGridHelpers::cellCenterDepth(grid, cell);

        if (ptable.oilActive())   { press[oilpos]  [ix] = ptable.oil  (depth); }
        if (ptable.gasActive())   { press[gaspos]  [ix] = ptable.gas  (depth); }
        if (ptable.waterActive()) { press[waterpos][ix] = ptable.water(depth); }

        ++ix;
    }
}

template <typename Grid, typename CellRange>
void verticalExtent(const Grid&           grid,
                    const CellRange&      cells,
                    int&                  cellcount,
                    std::array<double,2>& span)
{
    // This code is only supported in three space dimensions
    assert (Grid::dimensionworld == 3);

    span[0] = std::numeric_limits<double>::max();
    span[1] = std::numeric_limits<double>::lowest();
    cellcount = 0;

    const int nd = Grid::dimensionworld;

    // Define vertical span as
    //
    //   [minimum(node depth(cells)), maximum(node depth(cells))]
    //
    // Note: We use a sledgehammer approach--looping all
    // the nodes of all the faces of all the 'cells'--to
    // compute those bounds.  This necessarily entails
    // visiting some nodes (and faces) multiple times.
    //
    // Note: The implementation of 'RK4IVP<>' implicitly
    // imposes the requirement that cell centroids are all
    // within this vertical span.  That requirement is not
    // checked.
    auto cell2Faces = Opm::UgGridHelpers::cell2Faces(grid);
    auto faceVertices = Opm::UgGridHelpers::face2Vertices(grid);

    for (typename CellRange::const_iterator
             ci = cells.begin(), ce = cells.end();
         ci != ce; ++ci, ++cellcount)
    {
        for (auto fi = cell2Faces[*ci].begin(),
                  fe = cell2Faces[*ci].end();
             fi != fe; ++fi)
        {
            for (auto i = faceVertices[*fi].begin(),
                      e = faceVertices[*fi].end();
                 i != e; ++i)
            {
                const double z = Opm::UgGridHelpers::
                    vertexCoordinates(grid, *i)[nd - 1];

                if (z < span[0]) { span[0] = z; }
                if (z > span[1]) { span[1] = z; }
            }
        }
    }
}

} // namespace Details

/**
 * Compute initial phase pressures by means of equilibration.
 *
 * This function uses the information contained in an
 * equilibration record (i.e., depths and pressurs) as well as
 * a density calculator and related data to vertically
 * integrate the phase pressure ODE
 * \f[
 * \frac{\mathrm{d}p_{\alpha}}{\mathrm{d}z} =
 * \rho_{\alpha}(z,p_{\alpha})\cdot g
 * \f]
 * in which \f$\rho_{\alpha}$ denotes the fluid density of
 * fluid phase \f$\alpha\f$, \f$p_{\alpha}\f$ is the
 * corresponding phase pressure, \f$z\f$ is the depth and
 * \f$g\f$ is the acceleration due to gravity (assumed
 * directed downwords, in the positive \f$z\f$ direction).
 *
 * \tparam Region Type of an equilibration region information
 *                base.  Typically an instance of the EquilReg
 *                class template.
 *
 * \tparam CellRange Type of cell range that demarcates the
 *                cells pertaining to the current
 *                equilibration region.  Must implement
 *                methods begin() and end() to bound the range
 *                as well as provide an inner type,
 *                const_iterator, to traverse the range.
 *
 * \param[in] grid     Grid.
 * \param[in] reg   Current equilibration region.
 * \param[in] cells Range that spans the cells of the current
 *                  equilibration region.
 * \param[in] grav  Acceleration of gravity.
 *
 * \return Phase pressures, one vector for each active phase,
 * of pressure values in each cell in the current
 * equilibration region.
 */
template <class FluidSystem, class Grid, class Region, class CellRange>
std::vector< std::vector<double>>
phasePressures(const Grid&      grid,
               const Region&    reg,
               const CellRange& cells,
               const double     grav = Opm::unit::gravity)
{
    int ncell = 0;
    auto span = std::array<double, 2>{};
    Details::verticalExtent(grid, cells, ncell, span);

    using pval = std::vector<double>;
    auto press = std::vector<pval>
        (FluidSystem::numPhases, pval(ncell, 0.0));

    // Ensure gas/oil and oil/water contacts are within the span for the
    // phase pressure calculation.
    span[0] = std::min(span[0], std::min(reg.zgoc(), reg.zwoc()));
    span[1] = std::max(span[1], std::max(reg.zgoc(), reg.zwoc()));

    Details::equilibrateOWG<FluidSystem>(grid, reg, grav, span, cells, press);

    return press;
}

/**
 * Compute initial phase saturations by means of equilibration.
 *
 * \tparam FluidSystem  The FluidSystem from opm-material
 *                      Must be initialized before used.
 *
 * \tparam Grid   Type of the grid
 *
 * \tparam Region Type of an equilibration region information
 *                base.  Typically an instance of the EquilReg
 *                class template.
 *
 * \tparam CellRange Type of cell range that demarcates the
 *                cells pertaining to the current
 *                equilibration region.  Must implement
 *                methods begin() and end() to bound the range
 *                as well as provide an inner type,
 *                const_iterator, to traverse the range.
 *
 * \tparam MaterialLawManager The MaterialLawManager from opm-material
 *
 * \param[in] grid               Grid.
 * \param[in] reg             Current equilibration region.
 * \param[in] cells           Range that spans the cells of the current
 *                            equilibration region.
 * \param[in] materialLawManager   The MaterialLawManager from opm-material
 * \param[in] swatInit       A vector of initial water saturations.
 *                            The capillary pressure is scaled to fit these values
 * \param[in] phasePressures Phase pressures, one vector for each active phase,
 *                            of pressure values in each cell in the current
 *                            equilibration region.
 * \return                    Phase saturations, one vector for each phase, each containing
 *                            one saturation value per cell in the region.
 */
template <class FluidSystem, class Grid, class Region, class CellRange, class MaterialLawManager>
std::vector< std::vector<double>>
phaseSaturations(const Grid& grid,
                 const Region& reg,
                 const CellRange& cells,
                 MaterialLawManager& materialLawManager,
                 const std::vector<double> swatInit,
                 std::vector< std::vector<double> >& phasePressures)
{
    if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        throw std::runtime_error("Cannot initialise: not handling water-gas cases.");
    }

    std::vector< std::vector<double> > phaseSaturations = phasePressures; // Just to get the right size.

    // Adjust oil pressure according to gas saturation and cap pressure
    typedef Opm::SimpleModularFluidState<double,
                                         /*numPhases=*/3,
                                         /*numComponents=*/3,
                                         FluidSystem,
                                         /*storePressure=*/false,
                                         /*storeTemperature=*/false,
                                         /*storeComposition=*/false,
                                         /*storeFugacity=*/false,
                                         /*storeSaturation=*/true,
                                         /*storeDensity=*/false,
                                         /*storeViscosity=*/false,
                                         /*storeEnthalpy=*/false> MySatOnlyFluidState;

    MySatOnlyFluidState fluidState;
    typedef typename MaterialLawManager::MaterialLaw MaterialLaw;

    const bool water = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    const bool gas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const int oilpos = FluidSystem::oilPhaseIdx;
    const int waterpos = FluidSystem::waterPhaseIdx;
    const int gaspos = FluidSystem::gasPhaseIdx;
    std::vector<double>::size_type localIndex = 0;
    for (typename CellRange::const_iterator ci = cells.begin(); ci != cells.end(); ++ci, ++localIndex) {
        const int cell = *ci;
        const auto& scaledDrainageInfo =
            materialLawManager.oilWaterScaledEpsInfoDrainage(cell);
        const auto& matParams = materialLawManager.materialLawParams(cell);

        // Find saturations from pressure differences by
        // inverting capillary pressure functions.
        double sw = 0.0;
        if (water) {
            if (isConstPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,FluidSystem::waterPhaseIdx, cell)){
                const double cellDepth = Opm::UgGridHelpers::cellCenterDepth(grid,
                                                                             cell);
                sw = satFromDepth<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,cellDepth,reg.zwoc(),waterpos,cell,false);
                phaseSaturations[waterpos][localIndex] = sw;
            }
            else {
                const double pcov = phasePressures[oilpos][localIndex] - phasePressures[waterpos][localIndex];
                if (swatInit.empty()) { // Invert Pc to find sw
                    sw = satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, waterpos, cell, pcov);
                    phaseSaturations[waterpos][localIndex] = sw;
                }
                else { // Scale Pc to reflect imposed sw
                    sw = swatInit[cell];
                    sw = materialLawManager.applySwatinit(cell, pcov, sw);
                    phaseSaturations[waterpos][localIndex] = sw;
                }
            }
        }
        double sg = 0.0;
        if (gas) {
            if (isConstPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,FluidSystem::gasPhaseIdx,cell)){
                const double cellDepth = Opm::UgGridHelpers::cellCenterDepth(grid,
                                                                             cell);
                sg = satFromDepth<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,cellDepth,reg.zgoc(),gaspos,cell,true);
                phaseSaturations[gaspos][localIndex] = sg;
            }
            else {
                // Note that pcog is defined to be (pg - po), not (po - pg).
                const double pcog = phasePressures[gaspos][localIndex] - phasePressures[oilpos][localIndex];
                const double increasing = true; // pcog(sg) expected to be increasing function
                sg = satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, gaspos, cell, pcog, increasing);
                phaseSaturations[gaspos][localIndex] = sg;
            }
        }
        if (gas && water && (sg + sw > 1.0)) {
            // Overlapping gas-oil and oil-water transition
            // zones can lead to unphysical saturations when
            // treated as above. Must recalculate using gas-water
            // capillary pressure.
            const double pcgw = phasePressures[gaspos][localIndex] - phasePressures[waterpos][localIndex];
            if (! swatInit.empty()) {
                // Re-scale Pc to reflect imposed sw for vanishing oil phase.
                // This seems consistent with ecl, and fails to honour
                // swatInit in case of non-trivial gas-oil cap pressure.
                sw = materialLawManager.applySwatinit(cell, pcgw, sw);
            }
            sw = satFromSumOfPcs<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, waterpos, gaspos, cell, pcgw);
            sg = 1.0 - sw;
            phaseSaturations[waterpos][localIndex] = sw;
            phaseSaturations[gaspos][localIndex] = sg;
            if (water) {
                fluidState.setSaturation(FluidSystem::waterPhaseIdx, sw);
            }
            else {
                fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
            }
            fluidState.setSaturation(FluidSystem::oilPhaseIdx, 1.0 - sw - sg);
            fluidState.setSaturation(FluidSystem::gasPhaseIdx, sg);

            double pC[/*numPhases=*/3] = { 0.0, 0.0, 0.0 };
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            double pcGas = pC[FluidSystem::oilPhaseIdx] + pC[FluidSystem::gasPhaseIdx];
            phasePressures[oilpos][localIndex] = phasePressures[gaspos][localIndex] - pcGas;
        }
        phaseSaturations[oilpos][localIndex] = 1.0 - sw - sg;

        // Adjust phase pressures for max and min saturation ...
        double thresholdSat = 1.0e-6;

        double so = 1.0;
        double pC[FluidSystem::numPhases] = { 0.0, 0.0, 0.0 };
        if (water) {
            double swu = scaledDrainageInfo.Swu;
            fluidState.setSaturation(FluidSystem::waterPhaseIdx, swu);
            so -= swu;
        }
        if (gas) {
            double sgu = scaledDrainageInfo.Sgu;
            fluidState.setSaturation(FluidSystem::gasPhaseIdx, sgu);
            so-= sgu;
        }
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, so);

        if (water && sw > scaledDrainageInfo.Swu-thresholdSat) {
            fluidState.setSaturation(FluidSystem::waterPhaseIdx, scaledDrainageInfo.Swu);
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            double pcWat = pC[FluidSystem::oilPhaseIdx] - pC[FluidSystem::waterPhaseIdx];
            phasePressures[oilpos][localIndex] = phasePressures[waterpos][localIndex] + pcWat;
        }
        else if (gas && sg > scaledDrainageInfo.Sgu-thresholdSat) {
            fluidState.setSaturation(FluidSystem::gasPhaseIdx, scaledDrainageInfo.Sgu);
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            double pcGas = pC[FluidSystem::oilPhaseIdx] + pC[FluidSystem::gasPhaseIdx];
            phasePressures[oilpos][localIndex] = phasePressures[gaspos][localIndex] - pcGas;
        }
        if (gas && sg < scaledDrainageInfo.Sgl+thresholdSat) {
            fluidState.setSaturation(FluidSystem::gasPhaseIdx, scaledDrainageInfo.Sgl);
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            double pcGas = pC[FluidSystem::oilPhaseIdx] + pC[FluidSystem::gasPhaseIdx];
            phasePressures[gaspos][localIndex] = phasePressures[oilpos][localIndex] + pcGas;
        }
        if (water && sw < scaledDrainageInfo.Swl+thresholdSat) {
            fluidState.setSaturation(FluidSystem::waterPhaseIdx, scaledDrainageInfo.Swl);
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            double pcWat = pC[FluidSystem::oilPhaseIdx] - pC[FluidSystem::waterPhaseIdx];
            phasePressures[waterpos][localIndex] = phasePressures[oilpos][localIndex] - pcWat;
        }
    }
    return phaseSaturations;
}

/**
 * Compute initial Rs values.
 *
 * \tparam CellRangeType Type of cell range that demarcates the
 *                cells pertaining to the current
 *                equilibration region.  Must implement
 *                methods begin() and end() to bound the range
 *                as well as provide an inner type,
 *                const_iterator, to traverse the range.
 *
 * \param[in] grid            Grid.
 * \param[in] cells           Range that spans the cells of the current
 *                            equilibration region.
 * \param[in] oilPressure    Oil pressure for each cell in range.
 * \param[in] temperature     Temperature for each cell in range.
 * \param[in] rsFunc         Rs as function of pressure and depth.
 * \return                    Rs values, one for each cell in the 'cells' range.
 */
template <class Grid, class CellRangeType>
std::vector<double> computeRs(const Grid& grid,
                              const CellRangeType& cells,
                              const std::vector<double> oilPressure,
                              const std::vector<double>& temperature,
                              const Miscibility::RsFunction& rsFunc,
                              const std::vector<double> gasSaturation)
{
    assert(Grid::dimensionworld == 3);
    std::vector<double> rs(cells.size());
    int count = 0;
    for (auto it = cells.begin(); it != cells.end(); ++it, ++count) {
        const double depth = Opm::UgGridHelpers::cellCenterDepth(grid, *it);
        rs[count] = rsFunc(depth, oilPressure[count], temperature[count], gasSaturation[count]);
    }
    return rs;
}

namespace DeckDependent {
inline std::vector<Opm::EquilRecord>
getEquil(const Opm::EclipseState& state)
{
    const auto& init = state.getInitConfig();

    if(!init.hasEquil()) {
        throw std::domain_error("Deck does not provide equilibration data.");
    }

    const auto& equil = init.getEquil();
    return { equil.begin(), equil.end() };
}

template<class Grid>
std::vector<int>
equilnum(const Opm::EclipseState& eclipseState,
         const Grid& grid)
{
    std::vector<int> eqlnum(grid.size(0), 0);

    if (eclipseState.fieldProps().has_int("EQLNUM")) {
        const int nc = grid.size(/*codim=*/0);
        eqlnum.resize(nc);

        const auto& e = eclipseState.fieldProps().get_int("EQLNUM");
        std::transform(e.begin(), e.end(), eqlnum.begin(), [](int n){ return n - 1;});
    }
    return eqlnum;
}

template<class TypeTag>
class InitialStateComputer
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

public:
    template<class MaterialLawManager>
    InitialStateComputer(MaterialLawManager& materialLawManager,
                         const Opm::EclipseState& eclipseState,
                         const Grid& grid,
                         const double grav = Opm::unit::gravity,
                         const bool applySwatInit = true)
        : temperature_(grid.size(/*codim=*/0)),
          pp_(FluidSystem::numPhases,
              std::vector<double>(grid.size(/*codim=*/0))),
          sat_(FluidSystem::numPhases,
               std::vector<double>(grid.size(/*codim=*/0))),
          rs_(grid.size(/*codim=*/0)),
          rv_(grid.size(/*codim=*/0))
    {
        //Check for presence of kw SWATINIT
        if (applySwatInit) {
            if (eclipseState.fieldProps().has_double("SWATINIT")) {
                swatInit_ = eclipseState.fieldProps().get_double("SWATINIT");
            }
        }


        // Get the equilibration records.
        const std::vector<Opm::EquilRecord> rec = getEquil(eclipseState);
        const auto& tables = eclipseState.getTableManager();
        // Create (inverse) region mapping.
        const Opm::RegionMapping<> eqlmap(equilnum(eclipseState, grid));
        const int invalidRegion = -1;
        regionPvtIdx_.resize(rec.size(), invalidRegion);
        setRegionPvtIdx(grid, eclipseState, eqlmap);

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
                    const Opm::TableContainer& rsvdTables = tables.getRsvdTables();
                    const Opm::TableContainer& pbvdTables = tables.getPbvdTables();
                    if (rsvdTables.size() > 0) {

                        const Opm::RsvdTable& rsvdTable = rsvdTables.getTable<Opm::RsvdTable>(i);
                        std::vector<double> depthColumn = rsvdTable.getColumn("DEPTH").vectorCopy();
                        std::vector<double> rsColumn = rsvdTable.getColumn("RS").vectorCopy();
                        rsFunc_.push_back(std::make_shared<Miscibility::RsVD<FluidSystem>>(pvtIdx,
                                                                                           depthColumn, rsColumn));
                    } else if (pbvdTables.size() > 0) {
                        const Opm::PbvdTable& pbvdTable = pbvdTables.getTable<Opm::PbvdTable>(i);
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
                    const Opm::TableContainer& rvvdTables = tables.getRvvdTables();
                    const Opm::TableContainer& pdvdTables = tables.getPdvdTables();

                    if (rvvdTables.size() > 0) {
                        const Opm::RvvdTable& rvvdTable = rvvdTables.getTable<Opm::RvvdTable>(i);
                        std::vector<double> depthColumn = rvvdTable.getColumn("DEPTH").vectorCopy();
                        std::vector<double> rvColumn = rvvdTable.getColumn("RV").vectorCopy();
                        rvFunc_.push_back(std::make_shared<Miscibility::RvVD<FluidSystem>>(pvtIdx,
                                                                                           depthColumn, rvColumn));
                    } else if (pdvdTables.size() > 0) {
                        const Opm::PdvdTable& pdvdTable = pdvdTables.getTable<Opm::PdvdTable>(i);
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

        // Compute pressures, saturations, rs and rv factors.
        calcPressSatRsRv(eclipseState, eqlmap, rec, materialLawManager, grid, grav);

        // Modify oil pressure in no-oil regions so that the pressures of present phases can
        // be recovered from the oil pressure and capillary relations.
    }

    typedef std::vector<double> Vec;
    typedef std::vector<Vec>    PVec; // One per phase.

    const Vec& temperature() const { return temperature_; }
    const PVec& press() const { return pp_; }
    const PVec& saturation() const { return sat_; }
    const Vec& rs() const { return rs_; }
    const Vec& rv() const { return rv_; }

private:
    void updateInitialTemperature_(const Opm::EclipseState& eclState)
    {
        // Get the initial temperature data
        std::vector<double> tempiData = eclState.fieldProps().get_double("TEMPI");
        temperature_ = tempiData;
    }

    typedef EquilReg EqReg;
    std::vector< std::shared_ptr<Miscibility::RsFunction> > rsFunc_;
    std::vector< std::shared_ptr<Miscibility::RsFunction> > rvFunc_;
    std::vector<int> regionPvtIdx_;
    Vec temperature_;
    PVec pp_;
    PVec sat_;
    Vec rs_;
    Vec rv_;
    Vec swatInit_;

    template<class RMap>
    void setRegionPvtIdx(const Grid& grid, const Opm::EclipseState& eclState, const RMap& reg)
    {
        size_t numCompressed = grid.size(/*codim=*/0);
        std::vector<int> cellPvtRegionIdx(numCompressed);

        //Get the PVTNUM data
        const auto& pvtnumData = eclState.fieldProps().get_int("PVTNUM");
        // Save pvt indices of regions. Remember
        // that Eclipse uses Fortran-style indices which start at 1 instead of 0, so we
        // need to subtract 1.
        std::transform(pvtnumData.begin(), pvtnumData.end(),
                       cellPvtRegionIdx.begin(), [](int n){ return n - 1;});

        for (const auto& r : reg.activeRegions()) {
            const auto& cells = reg.cells(r);
            const int cell = *(cells.begin());
            regionPvtIdx_[r] = pvtnumData[cell] - 1;
        }
    }

    template <class RMap, class MaterialLawManager>
    void calcPressSatRsRv(const Opm::EclipseState& eclState OPM_UNUSED,
                          const RMap& reg,
                          const std::vector< Opm::EquilRecord >& rec,
                          MaterialLawManager& materialLawManager,
                          const Grid& grid,
                          const double grav)
    {
        for (const auto& r : reg.activeRegions()) {
            const auto& cells = reg.cells(r);
            if (cells.empty()) {
                Opm::OpmLog::warning("Equilibration region " + std::to_string(r + 1)
                                     + " has no active cells");
                continue;
            }

            const EqReg eqreg(rec[r], rsFunc_[r], rvFunc_[r], regionPvtIdx_[r]);

            PVec pressures = phasePressures<FluidSystem>(grid, eqreg, cells, grav);
            const PVec sat = phaseSaturations<FluidSystem>(grid, eqreg, cells, materialLawManager, swatInit_, pressures);

            const int np = FluidSystem::numPhases;
            for (int p = 0; p < np; ++p) {
                copyFromRegion(pressures[p], cells, pp_[p]);
                copyFromRegion(sat[p], cells, sat_[p]);
            }
            const bool oil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
            const bool gas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
            if (oil && gas) {
                const int oilpos = FluidSystem::oilPhaseIdx;
                const int gaspos = FluidSystem::gasPhaseIdx;
                std::vector<double> regionTemperature(cells.size());
                copyToRegion(temperature_, cells, regionTemperature);
                const Vec rsVals = computeRs(grid, cells, pressures[oilpos], regionTemperature, *(rsFunc_[r]), sat[gaspos]);
                const Vec rvVals = computeRs(grid, cells, pressures[gaspos], regionTemperature, *(rvFunc_[r]), sat[oilpos]);
                copyFromRegion(rsVals, cells, rs_);
                copyFromRegion(rvVals, cells, rv_);
            }
        }
    }

    template <class CellRangeType>
    void copyFromRegion(const Vec& source,
                        const CellRangeType& cells,
                        Vec& destination)
    {
        auto s = source.begin();
        auto c = cells.begin();
        const auto e = cells.end();
        for (; c != e; ++c, ++s) {
            destination[*c] =*s;
        }
    }

    template <class CellRangeType>
    void copyToRegion(const Vec& source,
                        const CellRangeType& cells,
                        Vec& destination)
    {
        auto d = destination.begin();
        auto c = cells.begin();
        const auto e = cells.end();
        for (; c != e; ++c, ++d) {
            *d = source[*c];
        }
    }
};
} // namespace DeckDependent
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
