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
#include "regionmapping.hh"

#include <ewoms/common/propertysystem.hh>

#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
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

#include <array>
#include <cassert>
#include <utility>
#include <vector>

BEGIN_PROPERTIES

NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(FluidSystem);

END_PROPERTIES

namespace Ewoms {

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


namespace PhasePressure {
template <class Grid,
          class PressFunction,
          class CellRange>
void assign(const Grid& grid,
            const std::array<PressFunction, 2>& f    ,
            const double split,
            const CellRange& cells,
            std::vector<double>& p)
{

    enum { up = 0, down = 1 };

    std::vector<double>::size_type c = 0;
    for (typename CellRange::const_iterator
             ci = cells.begin(), ce = cells.end();
         ci != ce; ++ci, ++c)
    {
        assert (c < p.size());

        const double z = Opm::UgGridHelpers::cellCenterDepth(grid, *ci);
        p[c] = (z < split) ? f[up](z) : f[down](z);
    }
}

template <class FluidSystem,
          class Grid,
          class Region,
          class CellRange>
void water(const Grid& grid,
           const Region& reg,
           const std::array<double,2>& span  ,
           const double grav,
           double& poWoc,
           const CellRange& cells,
           std::vector<double>& press)
{
    using PhasePressODE::Water;
    typedef Water<FluidSystem> ODE;

    const double T = 273.15 + 20; // standard temperature for now
    ODE drho(T, reg.pvtIdx() , grav);

    double z0;
    double p0;
    if (reg.datum() > reg.zwoc()) {//Datum in water zone
        z0 = reg.datum();
        p0 = reg.pressure();
    }
    else {
        z0 = reg.zwoc();
        p0 = poWoc - reg.pcowWoc(); // Water pressure at contact
    }

    std::array<double,2> up = {{ z0, span[0] }};
    std::array<double,2> down = {{ z0, span[1] }};

    typedef Details::RK4IVP<ODE> WPress;
    std::array<WPress,2> wpress = {
        {
            WPress(drho, up  , p0, 2000)
            ,
            WPress(drho, down, p0, 2000)
        }
    };

    assign(grid, wpress, z0, cells, press);

    if (reg.datum() > reg.zwoc()) {
        // Return oil pressure at contact
        poWoc = wpress[0](reg.zwoc()) + reg.pcowWoc();
    }
}

template <class FluidSystem,
          class Grid,
          class Region,
          class CellRange>
void oil(const Grid& grid,
         const Region& reg,
         const std::array<double,2>& span  ,
         const double grav,
         const CellRange& cells,
         std::vector<double>& press,
         double& poWoc,
         double& poGoc)
{
    using PhasePressODE::Oil;
    typedef Oil<FluidSystem, typename Region::CalcDissolution> ODE;

    const double T = 273.15 + 20; // standard temperature for now
    ODE drho(T, reg.dissolutionCalculator(),
             reg.pvtIdx(), grav);

    double z0;
    double p0;
    if (reg.datum() > reg.zwoc()) {//Datum in water zone, poWoc given
        z0 = reg.zwoc();
        p0 = poWoc;
    }
    else if (reg.datum() < reg.zgoc()) {//Datum in gas zone, poGoc given
        z0 = reg.zgoc();
        p0 = poGoc;
    }
    else { //Datum in oil zone
        z0 = reg.datum();
        p0 = reg.pressure();
    }

    std::array<double,2> up = {{ z0, span[0] }};
    std::array<double,2> down = {{ z0, span[1] }};

    typedef Details::RK4IVP<ODE> OPress;
    std::array<OPress,2> opress = {
        {
            OPress(drho, up  , p0, 2000)
            ,
            OPress(drho, down, p0, 2000)
        }
    };

    assign(grid, opress, z0, cells, press);

    const double woc = reg.zwoc();
    if      (z0 > woc) { poWoc = opress[0](woc); } // WOC above datum
    else if (z0 < woc) { poWoc = opress[1](woc); } // WOC below datum
    else               { poWoc = p0;             } // WOC *at*  datum

    const double goc = reg.zgoc();
    if      (z0 > goc) { poGoc = opress[0](goc); } // GOC above datum
    else if (z0 < goc) { poGoc = opress[1](goc); } // GOC below datum
    else               { poGoc = p0;             } // GOC *at*  datum
}

template <class FluidSystem,
          class Grid,
          class Region,
          class CellRange>
void gas(const Grid& grid,
         const Region& reg,
         const std::array<double,2>& span  ,
         const double grav,
         double& poGoc,
         const CellRange& cells,
         std::vector<double>& press)
{
    using PhasePressODE::Gas;
    typedef Gas<FluidSystem, typename Region::CalcEvaporation> ODE;

    const double T = 273.15 + 20; // standard temperature for now
    ODE drho(T, reg.evaporationCalculator(),
             reg.pvtIdx(), grav);

    double z0;
    double p0;
    if (reg.datum() < reg.zgoc()) {//Datum in gas zone
        z0 = reg.datum();
        p0 = reg.pressure();
    }
    else {
        z0 = reg.zgoc();
        p0 = poGoc + reg.pcgoGoc(); // Gas pressure at contact
    }

    std::array<double,2> up = {{ z0, span[0] }};
    std::array<double,2> down = {{ z0, span[1] }};

    typedef Details::RK4IVP<ODE> GPress;
    std::array<GPress,2> gpress = {
        {
            GPress(drho, up  , p0, 2000)
            ,
            GPress(drho, down, p0, 2000)
        }
    };

    assign(grid, gpress, z0, cells, press);

    if (reg.datum() < reg.zgoc()) {
        // Return oil pressure at contact
        poGoc = gpress[1](reg.zgoc()) - reg.pcgoGoc();
    }
}
} // namespace PhasePressure

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
    const bool water = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    const bool oil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    const bool gas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const int oilpos = FluidSystem::oilPhaseIdx;
    const int waterpos = FluidSystem::waterPhaseIdx;
    const int gaspos = FluidSystem::gasPhaseIdx;

    if (reg.datum() > reg.zwoc()) { // Datum in water zone
        double poWoc = -1;
        double poGoc = -1;

        if (water) {
            PhasePressure::water<FluidSystem>(grid, reg, span, grav, poWoc,
                                              cells, press[waterpos]);
        }

        if (oil) {
            PhasePressure::oil<FluidSystem>(grid, reg, span, grav, cells,
                                            press[oilpos], poWoc, poGoc);
        }

        if (gas) {
            PhasePressure::gas<FluidSystem>(grid, reg, span, grav, poGoc,
                                            cells, press[gaspos]);
        }
    }
    else if (reg.datum() < reg.zgoc()) { // Datum in gas zone
        double poWoc = -1;
        double poGoc = -1;

        if (gas) {
            PhasePressure::gas<FluidSystem>(grid, reg, span, grav, poGoc,
                                            cells, press[gaspos]);
        }

        if (oil) {
            PhasePressure::oil<FluidSystem>(grid, reg, span, grav, cells,
                                            press[oilpos], poWoc, poGoc);
        }

        if (water) {
            PhasePressure::water<FluidSystem>(grid, reg, span, grav, poWoc,
                                              cells, press[waterpos]);
        }
    }
    else { // Datum in oil zone
        double poWoc = -1;
        double poGoc = -1;

        if (oil) {
            PhasePressure::oil<FluidSystem>(grid, reg, span, grav, cells,
                                            press[oilpos], poWoc, poGoc);
        }

        if (water) {
            PhasePressure::water<FluidSystem>(grid, reg, span, grav, poWoc,
                                              cells, press[waterpos]);
        }

        if (gas) {
            PhasePressure::gas<FluidSystem>(grid, reg, span, grav, poGoc,
                                            cells, press[gaspos]);
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
phasePressures(const Grid& grid,
               const Region& reg,
               const CellRange& cells,
               const double grav = Opm::unit::gravity)
{
    std::array<double,2> span =
        {{  std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max() }}; // Symm. about 0.

    int ncell = 0;
    {
        // This code is only supported in three space dimensions
        assert (Grid::dimensionworld == 3);

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
             ci != ce; ++ci, ++ncell)
        {
            for (auto fi=cell2Faces[*ci].begin(),
                     fe=cell2Faces[*ci].end();
                 fi != fe;
                 ++fi)
            {
                for (auto i = faceVertices[*fi].begin(), e = faceVertices[*fi].end();
                     i != e; ++i)
                {
                    const double z = Opm::UgGridHelpers::vertexCoordinates(grid, *i)[nd-1];

                    if (z < span[0]) { span[0] = z; }
                    if (z > span[1]) { span[1] = z; }
                }
            }
        }
    }
    const int np = FluidSystem::numPhases;  //reg.phaseUsage().numPhases;

    typedef std::vector<double> pval;
    std::vector<pval> press(np, pval(ncell, 0.0));

    const double zwoc = reg.zwoc ();
    const double zgoc = reg.zgoc ();

    // make sure goc and woc is within the span for the phase pressure calculation
    span[0] = std::min(span[0],zgoc);
    span[1] = std::max(span[1],zwoc);

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
    std::vector<int> eqlnum;
    if (eclipseState.get3DProperties().hasDeckIntGridProperty("EQLNUM")) {
        const int nc = grid.size(/*codim=*/0);
        eqlnum.resize(nc);
        const std::vector<int>& e =
            eclipseState.get3DProperties().getIntGridProperty("EQLNUM").getData();
        const int* gc = Opm::UgGridHelpers::globalCell(grid);
        for (int cell = 0; cell < nc; ++cell) {
            const int deckPos = (gc == NULL) ? cell : gc[cell];
            eqlnum[cell] = e[deckPos] - 1;
        }
    }
    else {
        // No explicit equilibration region.
        // All cells in region zero.
        eqlnum.assign(grid.size(/*codim=*/0), 0);
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
        if (eclipseState.get3DProperties().hasDeckDoubleGridProperty("SWATINIT") && applySwatInit) {
            const std::vector<double>& swatInitEcl = eclipseState.
                get3DProperties().getDoubleGridProperty("SWATINIT").getData();
            const int nc = grid.size(/*codim=*/0);
            swatInit_.resize(nc);
            const int* gc = Opm::UgGridHelpers::globalCell(grid);
            for (int c = 0; c < nc; ++c) {
                const int deckPos = (gc == NULL) ? c : gc[c];
                swatInit_[c] = swatInitEcl[deckPos];
            }
        }
        // Get the equilibration records.
        const std::vector<Opm::EquilRecord> rec = getEquil(eclipseState);
        const auto& tables = eclipseState.getTableManager();
        // Create (inverse) region mapping.
        const Ewoms::RegionMapping<> eqlmap(equilnum(eclipseState, grid));
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

        // extract the initial temperature
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
        const std::vector<double>& tempiData =
            eclState.get3DProperties().getDoubleGridProperty("TEMPI").getData();

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
        const auto* globalCell = Opm::UgGridHelpers::globalCell(grid);
        std::vector<int> cellPvtRegionIdx(numCompressed);

        //Get the PVTNUM data
        const std::vector<int>& pvtnumData = eclState.get3DProperties().getIntGridProperty("PVTNUM").getData();

        // Convert PVTNUM data into an array of indices for compressed cells. Remember
        // that Eclipse uses Fortran-style indices which start at 1 instead of 0, so we
        // need to subtract 1.
        for (size_t cellIdx = 0; cellIdx < numCompressed; ++ cellIdx) {
            size_t cartesianCellIdx = globalCell[cellIdx];
            assert(cartesianCellIdx < pvtnumData.size());
            size_t pvtRegionIdx = pvtnumData[cartesianCellIdx] - 1;
            cellPvtRegionIdx[cellIdx] = pvtRegionIdx;
        }

        for (const auto& r : reg.activeRegions()) {
            const auto& cells = reg.cells(r);
            const int cell = *(cells.begin());
            regionPvtIdx_[r] = cellPvtRegionIdx[cell];
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
                const Vec rsVals = computeRs(grid, cells, pressures[oilpos], temperature_, *(rsFunc_[r]), sat[gaspos]);
                const Vec rvVals = computeRs(grid, cells, pressures[gaspos], temperature_, *(rvFunc_[r]), sat[oilpos]);
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

};
} // namespace DeckDependent
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
