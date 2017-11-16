/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#ifndef OPM_INITSTATEEQUIL_IMPL_HEADER_INCLUDED
#define OPM_INITSTATEEQUIL_IMPL_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/grid/GridHelpers.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

namespace Opm
{
    namespace Details {



        template <class RHS>
        class RK4IVP {
        public:
            RK4IVP(const RHS&                  f   ,
                   const std::array<double,2>& span,
                   const double                y0  ,
                   const int                   N   )
                : N_(N)
                , span_(span)
            {
                const double h  = stepsize();
                const double h2 = h / 2;
                const double h6 = h / 6;

                y_.reserve(N + 1);
                f_.reserve(N + 1);

                y_.push_back(y0);
                f_.push_back(f(span_[0], y0));

                for (int i = 0; i < N; ++i) {
                    const double x  = span_[0] + i*h;
                    const double y  = y_.back();

                    const double k1 = f_[i];
                    const double k2 = f(x + h2, y + h2*k1);
                    const double k3 = f(x + h2, y + h2*k2);
                    const double k4 = f(x + h , y + h*k3);

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
                int          i = (x - span_[0]) / h;
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
            int                  N_;
            std::array<double,2> span_;
            std::vector<double>  y_;
            std::vector<double>  f_;

            double
            stepsize() const { return (span_[1] - span_[0]) / N_; }
        };

         namespace PhasePressODE {
          template <class FluidSystem>
            class Water {
            public:
                Water(const double   temp,
                      const int      pvtRegionIdx,
                      const double   norm_grav)
                    : temp_(temp)
                    , pvtRegionIdx_(pvtRegionIdx)
                    , g_(norm_grav)
                {
                }

                double
                operator()(const double /* depth */,
                           const double press) const
                {
                    return this->density(press) * g_;
                }

            private:
                const double        temp_;
                const int           pvtRegionIdx_;
                const double        g_;

                double
                density(const double press) const
                {
                    double rho = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
                    rho *= FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx_);
                    return rho;
                }
            };

            template <class FluidSystem, class RS>
            class Oil {
            public:
                Oil(const double   temp,
                    const RS&      rs,
                    const int      pvtRegionIdx,
                    const double   norm_grav)
                    : temp_(temp)
                    , rs_(rs)
                    , pvtRegionIdx_(pvtRegionIdx)
                    , g_(norm_grav)
                {
                }

                double
                operator()(const double depth,
                           const double press) const
                {
                    return this->density(depth, press) * g_;
                }

            private:
                const double                temp_;
                const RS&                   rs_;
                const int                   pvtRegionIdx_;
                const double                g_;

                double
                density(const double depth,
                        const double press) const
                {
                    double rs = rs_(depth, press, temp_);
                    double bOil = 0.0;
                    if ( !FluidSystem::enableDissolvedGas() || rs >= FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp_, press) ) {
                        bOil = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
                    } else {
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
            class Gas {
            public:
                Gas(const double   temp,
                    const RV&      rv,
                    const int      pvtRegionIdx,
                    const double   norm_grav)
                    : temp_(temp)
                    , rv_(rv)
                    , pvtRegionIdx_(pvtRegionIdx)
                    , g_(norm_grav)
                {
                }

                double
                operator()(const double depth,
                           const double press) const
                {
                    return this->density(depth, press) * g_;
                }

            private:
                const double                temp_;
                const RV&                   rv_;
                const int                   pvtRegionIdx_;
                const double                g_;

                double
                density(const double depth,
                        const double press) const
                {
                    double rv = rv_(depth, press, temp_);
                    double bGas = 0.0;
                    if ( !FluidSystem::enableVaporizedOil() || rv >= FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp_, press) ) {
                        bGas = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx_, temp_, press);
                    } else {
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

        namespace PhaseUsed {
            inline bool
            water(const PhaseUsage& pu)
            {
                return bool(pu.phase_used[ Opm::BlackoilPhases::Aqua ]);
            }

            inline bool
            oil(const PhaseUsage& pu)
            {
                return bool(pu.phase_used[ Opm::BlackoilPhases::Liquid ]);
            }

            inline bool
            gas(const PhaseUsage& pu)
            {
                return bool(pu.phase_used[ Opm::BlackoilPhases::Vapour ]);
            }
        } // namespace PhaseUsed

        namespace PhaseIndex {
            inline int
            water(const PhaseUsage& pu)
            {
                int i = -1;
                if (PhaseUsed::water(pu)) {
                    i = pu.phase_pos[ Opm::BlackoilPhases::Aqua ];
                }

                return i;
            }

            inline int
            oil(const PhaseUsage& pu)
            {
                int i = -1;
                if (PhaseUsed::oil(pu)) {
                    i = pu.phase_pos[ Opm::BlackoilPhases::Liquid ];
                }

                return i;
            }

            inline int
            gas(const PhaseUsage& pu)
            {
                int i = -1;
                if (PhaseUsed::gas(pu)) {
                    i = pu.phase_pos[ Opm::BlackoilPhases::Vapour ];
                }

                return i;
            }
        } // namespace PhaseIndex

        namespace PhasePressure {
            template <class Grid,
                      class PressFunction,
                      class CellRange>
            void
            assign(const Grid&                         G    ,
                   const std::array<PressFunction, 2>& f    ,
                   const double                        split,
                   const CellRange&                    cells,
                   std::vector<double>&                p    )
            {

                enum { up = 0, down = 1 };

                std::vector<double>::size_type c = 0;
                for (typename CellRange::const_iterator
                         ci = cells.begin(), ce = cells.end();
                     ci != ce; ++ci, ++c)
                {
                    assert (c < p.size());

                    const double z = UgGridHelpers::cellCenterDepth(G, *ci);
                    p[c] = (z < split) ? f[up](z) : f[down](z);
                }
            }

            template <                     class FluidSystem,
                                           class Grid,
                                           class Region,
                                           class CellRange>
            void
            water(const Grid&                 G     ,
                  const Region&               reg   ,
                  const std::array<double,2>& span  ,
                  const double                grav  ,
                  double&                     po_woc,
                  const CellRange&            cells ,
                  std::vector<double>&        press )
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
                } else {
                    z0 = reg.zwoc();
                    p0 = po_woc - reg.pcow_woc(); // Water pressure at contact
                }

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> WPress;
                std::array<WPress,2> wpress = {
                    {
                        WPress(drho, up  , p0, 2000)
                        ,
                        WPress(drho, down, p0, 2000)
                    }
                };

                assign(G, wpress, z0, cells, press);

                if (reg.datum() > reg.zwoc()) {
                    // Return oil pressure at contact
                    po_woc = wpress[0](reg.zwoc()) + reg.pcow_woc();
                }
            }

            template <class FluidSystem,
                      class Grid,
                      class Region,
                      class CellRange>

            void
            oil(const Grid&                 G     ,
                const Region&               reg   ,
                const std::array<double,2>& span  ,
                const double                grav  ,
                const CellRange&            cells ,
                std::vector<double>&        press ,
                double&                     po_woc,
                double&                     po_goc)
            {
                using PhasePressODE::Oil;
                typedef Oil<FluidSystem, typename Region::CalcDissolution> ODE;

                const double T = 273.15 + 20; // standard temperature for now

                ODE drho(T, reg.dissolutionCalculator(),
                         reg.pvtIdx(), grav);

                double z0;
                double p0;
                if (reg.datum() > reg.zwoc()) {//Datum in water zone, po_woc given
                    z0 = reg.zwoc();
                    p0 = po_woc;
                } else if (reg.datum() < reg.zgoc()) {//Datum in gas zone, po_goc given
                    z0 = reg.zgoc();
                    p0 = po_goc;
                } else { //Datum in oil zone
                    z0 = reg.datum();
                    p0 = reg.pressure();
                }

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> OPress;
                std::array<OPress,2> opress = {
                    {
                        OPress(drho, up  , p0, 2000)
                        ,
                        OPress(drho, down, p0, 2000)
                    }
                };

                assign(G, opress, z0, cells, press);

                const double woc = reg.zwoc();
                if      (z0 > woc) { po_woc = opress[0](woc); } // WOC above datum
                else if (z0 < woc) { po_woc = opress[1](woc); } // WOC below datum
                else               { po_woc = p0;             } // WOC *at*  datum

                const double goc = reg.zgoc();
                if      (z0 > goc) { po_goc = opress[0](goc); } // GOC above datum
                else if (z0 < goc) { po_goc = opress[1](goc); } // GOC below datum
                else               { po_goc = p0;             } // GOC *at*  datum
            }

            template <class FluidSystem,
                      class Grid,
                      class Region,
                      class CellRange>
            void
            gas(const Grid&                 G     ,
                const Region&               reg   ,
                const std::array<double,2>& span  ,
                const double                grav  ,
                double&                     po_goc,
                const CellRange&            cells ,
                std::vector<double>&        press )
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
                } else {
                    z0 = reg.zgoc();
                    p0 = po_goc + reg.pcgo_goc(); // Gas pressure at contact
                }

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> GPress;
                std::array<GPress,2> gpress = {
                    {
                        GPress(drho, up  , p0, 2000)
                        ,
                        GPress(drho, down, p0, 2000)
                    }
                };

                assign(G, gpress, z0, cells, press);

                if (reg.datum() < reg.zgoc()) {
                    // Return oil pressure at contact
                    po_goc = gpress[1](reg.zgoc()) - reg.pcgo_goc();
                }
            }
        } // namespace PhasePressure

        template <class FluidSystem,
                  class Grid,
                  class Region,
                  class CellRange>
        void
        equilibrateOWG(const Grid&                         G,
                       const Region&                       reg,
                       const double                        grav,
                       const std::array<double,2>&         span,
                       const CellRange&                    cells,
                       std::vector< std::vector<double> >& press)
        {
            const PhaseUsage& pu = reg.phaseUsage();

            if (reg.datum() > reg.zwoc()) { // Datum in water zone
                double po_woc = -1;
                double po_goc = -1;

                if (PhaseUsed::water(pu)) {
                    const int wix = PhaseIndex::water(pu);
                    PhasePressure::water<FluidSystem>(G, reg, span, grav, po_woc,
                                         cells, press[ wix ]);
                }

                if (PhaseUsed::oil(pu)) {
                    const int oix = PhaseIndex::oil(pu);
                    PhasePressure::oil<FluidSystem>(G, reg, span, grav, cells,
                                       press[ oix ], po_woc, po_goc);
                }

                if (PhaseUsed::gas(pu)) {
                    const int gix = PhaseIndex::gas(pu);
                    PhasePressure::gas<FluidSystem>(G, reg, span, grav, po_goc,
                                       cells, press[ gix ]);
                }
            } else if (reg.datum() < reg.zgoc()) { // Datum in gas zone
                double po_woc = -1;
                double po_goc = -1;

                if (PhaseUsed::gas(pu)) {
                    const int gix = PhaseIndex::gas(pu);
                    PhasePressure::gas<FluidSystem>(G, reg, span, grav, po_goc,
                                       cells, press[ gix ]);
                }

                if (PhaseUsed::oil(pu)) {
                    const int oix = PhaseIndex::oil(pu);
                    PhasePressure::oil<FluidSystem>(G, reg, span, grav, cells,
                                       press[ oix ], po_woc, po_goc);
                }

                if (PhaseUsed::water(pu)) {
                    const int wix = PhaseIndex::water(pu);
                    PhasePressure::water<FluidSystem>(G, reg, span, grav, po_woc,
                                         cells, press[ wix ]);
                }
            } else { // Datum in oil zone
                double po_woc = -1;
                double po_goc = -1;

                if (PhaseUsed::oil(pu)) {
                    const int oix = PhaseIndex::oil(pu);
                    PhasePressure::oil<FluidSystem>(G, reg, span, grav, cells,
                                                    press[ oix ], po_woc, po_goc);
                }

                if (PhaseUsed::water(pu)) {
                    const int wix = PhaseIndex::water(pu);
                    PhasePressure::water<FluidSystem>(G, reg, span, grav, po_woc,
                                         cells, press[ wix ]);
                }

                if (PhaseUsed::gas(pu)) {
                    const int gix = PhaseIndex::gas(pu);
                    PhasePressure::gas<FluidSystem>(G, reg, span, grav, po_goc,
                                       cells, press[ gix ]);
                }
            }
        }
    } // namespace Details


    namespace EQUIL {


        template <class FluidSystem,
                  class Grid,
                  class Region,
                  class CellRange>
        std::vector< std::vector<double> >
        phasePressures(const Grid&             G,
                       const Region&           reg,
                       const CellRange&        cells,
                       const double            grav)
        {
            std::array<double,2> span =
                {{  std::numeric_limits<double>::max() ,
                   -std::numeric_limits<double>::max() }}; // Symm. about 0.

            int ncell = 0;
            {
                // This code is only supported in three space dimensions
                assert (UgGridHelpers::dimensions(G) == 3);

                const int nd = UgGridHelpers::dimensions(G);

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
                auto cell2Faces = UgGridHelpers::cell2Faces(G);
                auto faceVertices = UgGridHelpers::face2Vertices(G);

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
                            const double z = UgGridHelpers::vertexCoordinates(G, *i)[nd-1];

                            if (z < span[0]) { span[0] = z; }
                            if (z > span[1]) { span[1] = z; }
                        }
                    }
                }
            }
            const int np = reg.phaseUsage().num_phases;

            typedef std::vector<double> pval;
            std::vector<pval> press(np, pval(ncell, 0.0));

            const double zwoc = reg.zwoc ();
            const double zgoc = reg.zgoc ();

            // make sure goc and woc is within the span for the phase pressure calculation
            span[0] = std::min(span[0],zgoc);
            span[1] = std::max(span[1],zwoc);

            Details::equilibrateOWG<FluidSystem>(G, reg, grav, span, cells, press);

            return press;
        }

        template <class Grid,
                  class Region,
                  class CellRange>
        std::vector<double>
        temperature(const Grid&             /* G */,
                    const Region&           /* reg */,
                    const CellRange&        cells)
        {
            // use the standard temperature for everything for now
            return std::vector<double>(cells.size(), 273.15 + 20.0);
        }

        template <class FluidSystem, class Grid, class Region, class CellRange, class MaterialLawManager>
        std::vector< std::vector<double> >
        phaseSaturations(const Grid&             G,
                         const Region&           reg,
                         const CellRange&        cells,
                         std::shared_ptr<MaterialLawManager> materialLawManager,
                         const std::vector<double> swat_init,
                         std::vector< std::vector<double> >& phase_pressures)
        {
            if (!reg.phaseUsage().phase_used[BlackoilPhases::Liquid]) {
                OPM_THROW(std::runtime_error, "Cannot initialise: not handling water-gas cases.");
            }

            std::vector< std::vector<double> > phase_saturations = phase_pressures; // Just to get the right size.

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
                    /*storeEnthalpy=*/false> SatOnlyFluidState;

            SatOnlyFluidState fluidState;
            typedef typename MaterialLawManager::MaterialLaw MaterialLaw;

            const bool water = reg.phaseUsage().phase_used[BlackoilPhases::Aqua];
            const bool gas = reg.phaseUsage().phase_used[BlackoilPhases::Vapour];
            const int oilpos = reg.phaseUsage().phase_pos[BlackoilPhases::Liquid];
            const int waterpos = reg.phaseUsage().phase_pos[BlackoilPhases::Aqua];
            const int gaspos = reg.phaseUsage().phase_pos[BlackoilPhases::Vapour];
            std::vector<double>::size_type local_index = 0;
            for (typename CellRange::const_iterator ci = cells.begin(); ci != cells.end(); ++ci, ++local_index) {
                const int cell = *ci;
                const auto& scaledDrainageInfo =
                    materialLawManager->oilWaterScaledEpsInfoDrainage(cell);
                const auto& matParams = materialLawManager->materialLawParams(cell);

                // Find saturations from pressure differences by
                // inverting capillary pressure functions.
                double sw = 0.0;
                if (water) {
                    if (isConstPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,FluidSystem::waterPhaseIdx, cell)){
                        const double cellDepth  =  UgGridHelpers::cellCenterDepth(G,
                                                                            cell);
                        sw = satFromDepth<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,cellDepth,reg.zwoc(),waterpos,cell,false);
                        phase_saturations[waterpos][local_index] = sw;
                    }
                    else{
                        const double pcov = phase_pressures[oilpos][local_index] - phase_pressures[waterpos][local_index];
                        if (swat_init.empty()) { // Invert Pc to find sw
                            sw = satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, waterpos, cell, pcov);
                            phase_saturations[waterpos][local_index] = sw;
                        } else { // Scale Pc to reflect imposed sw
                            sw = swat_init[cell];
                            sw = materialLawManager->applySwatinit(cell, pcov, sw);
                            phase_saturations[waterpos][local_index] = sw;
                        }
                    }
                }
                double sg = 0.0;
                if (gas) {
                    if (isConstPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,FluidSystem::gasPhaseIdx,cell)){
                        const double cellDepth  = UgGridHelpers::cellCenterDepth(G,
                                                                                        cell);
                        sg = satFromDepth<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager,cellDepth,reg.zgoc(),gaspos,cell,true);
                        phase_saturations[gaspos][local_index] = sg;
                    }
                    else{
                        // Note that pcog is defined to be (pg - po), not (po - pg).
                        const double pcog = phase_pressures[gaspos][local_index] - phase_pressures[oilpos][local_index];
                        const double increasing = true; // pcog(sg) expected to be increasing function
                        sg = satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, gaspos, cell, pcog, increasing);
                        phase_saturations[gaspos][local_index] = sg;
                    }
                }
                if (gas && water && (sg + sw > 1.0)) {
                    // Overlapping gas-oil and oil-water transition
                    // zones can lead to unphysical saturations when
                    // treated as above. Must recalculate using gas-water
                    // capillary pressure.
                    const double pcgw = phase_pressures[gaspos][local_index] - phase_pressures[waterpos][local_index];
                    if (! swat_init.empty()) { 
                        // Re-scale Pc to reflect imposed sw for vanishing oil phase.
                        // This seems consistent with ecl, and fails to honour 
                        // swat_init in case of non-trivial gas-oil cap pressure.
                        sw = materialLawManager->applySwatinit(cell, pcgw, sw);
                    }
                    sw = satFromSumOfPcs<FluidSystem, MaterialLaw, MaterialLawManager>(materialLawManager, waterpos, gaspos, cell, pcgw);
                    sg = 1.0 - sw;
                    phase_saturations[waterpos][local_index] = sw;
                    phase_saturations[gaspos][local_index] = sg;

                    if ( water ) {
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
                    phase_pressures[oilpos][local_index] = phase_pressures[gaspos][local_index] - pcGas;
                }
                phase_saturations[oilpos][local_index] = 1.0 - sw - sg;
                
                // Adjust phase pressures for max and min saturation ...
                double threshold_sat = 1.0e-6;

                double so = 1.0;
                double pC[/*numPhases=*/BlackoilPhases::MaxNumPhases] = { 0.0, 0.0, 0.0 };

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

                  if (water && sw > scaledDrainageInfo.Swu-threshold_sat ) {
                       fluidState.setSaturation(FluidSystem::waterPhaseIdx, scaledDrainageInfo.Swu);
                       MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                       double pcWat = pC[FluidSystem::oilPhaseIdx] - pC[FluidSystem::waterPhaseIdx];
                      phase_pressures[oilpos][local_index] = phase_pressures[waterpos][local_index] + pcWat;
                  } else if (gas && sg > scaledDrainageInfo.Sgu-threshold_sat) {
                      fluidState.setSaturation(FluidSystem::gasPhaseIdx, scaledDrainageInfo.Sgu);
                      MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                      double pcGas = pC[FluidSystem::oilPhaseIdx] + pC[FluidSystem::gasPhaseIdx];
                      phase_pressures[oilpos][local_index] = phase_pressures[gaspos][local_index] - pcGas;
                  }
                  if (gas && sg < scaledDrainageInfo.Sgl+threshold_sat) {
                      fluidState.setSaturation(FluidSystem::gasPhaseIdx, scaledDrainageInfo.Sgl);
                      MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                      double pcGas = pC[FluidSystem::oilPhaseIdx] + pC[FluidSystem::gasPhaseIdx];
                      phase_pressures[gaspos][local_index] = phase_pressures[oilpos][local_index] + pcGas;
                  }
                  if (water && sw < scaledDrainageInfo.Swl+threshold_sat) {
                      fluidState.setSaturation(FluidSystem::waterPhaseIdx, scaledDrainageInfo.Swl);
                      MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                      double pcWat = pC[FluidSystem::oilPhaseIdx] - pC[FluidSystem::waterPhaseIdx];
                      phase_pressures[waterpos][local_index] = phase_pressures[oilpos][local_index] - pcWat;
                  }
            }
            return phase_saturations;
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
         * \param[in] oil_pressure    Oil pressure for each cell in range.
         * \param[in] temperature     Temperature for each cell in range.
         * \param[in] rs_func         Rs as function of pressure and depth.
         * \return                    Rs values, one for each cell in the 'cells' range.
         */
        template <class Grid, class CellRangeType>
        std::vector<double> computeRs(const Grid& grid,
                                      const CellRangeType& cells,
                                      const std::vector<double> oil_pressure,
                                      const std::vector<double>& temperature,
                                      const Miscibility::RsFunction& rs_func,
                                      const std::vector<double> gas_saturation)
        {
            assert(UgGridHelpers::dimensions(grid) == 3);
            std::vector<double> rs(cells.size());
            int count = 0;
            for (auto it = cells.begin(); it != cells.end(); ++it, ++count) {
                const double depth = UgGridHelpers::cellCenterDepth(grid, *it);
                rs[count] = rs_func(depth, oil_pressure[count], temperature[count], gas_saturation[count]);
            }
            return rs;
        }

    } // namespace Equil


    namespace Details
    {
        /// Convert saturations from a vector of individual phase saturation vectors
        /// to an interleaved format where all values for a given cell come before all
        /// values for the next cell, all in a single vector.
        inline std::vector<double>
        convertSats(const std::vector< std::vector<double> >& sat)
        {
            const auto np = sat.size();
            const auto nc = sat[0].size();

            std::vector<double> s(np * nc);

            for (decltype(sat.size()) p = 0; p < np; ++p) {
                const auto& sat_p = sat[p];
                double*     sp    = & s[0*nc + p];

                for (decltype(sat[0].size()) c = 0;
                     c < nc; ++c, sp += np)
                {
                    *sp = sat_p[c];
                }
            }

            return s;
        }
    } // namespace Details


    /**
     * Compute initial state by an equilibration procedure.
     *
     * The following state fields are modified:
     *   pressure(),
     *   saturation(),
     *   surfacevol(),
     *   gasoilratio(),
     *   rv().
     *
     * \param[in] grid     Grid.
     * \param[in] props    Property object, pvt and capillary properties are used.
     * \param[in] deck     Simulation deck, used to obtain EQUIL and related data.
     * \param[in] gravity  Acceleration of gravity, assumed to be in Z direction.
     * \param[in] applySwatInit     Make it possible to not apply SWATINIT even if it
     *                              is present in the deck
     */
    template<class MaterialLawManager, class Grid>
    void initStateEquil(const Grid& grid,
                        std::shared_ptr<MaterialLawManager> materialLawManager,
                        const Opm::Deck& deck,
                        const Opm::EclipseState& eclipseState,
                        const double gravity,
                        BlackoilState& state,
                        bool applySwatinit = true)
    {

        typedef EQUIL::DeckDependent::InitialStateComputer ISC;

        PhaseUsage pu = phaseUsageFromDeck(deck);

        //Check for presence of kw SWATINIT
        std::vector<double> swat_init = {};
        if (eclipseState.get3DProperties().hasDeckDoubleGridProperty("SWATINIT") && applySwatinit) {
            const std::vector<double>& swat_init_ecl = eclipseState.
                    get3DProperties().getDoubleGridProperty("SWATINIT").getData();
            const int nc = UgGridHelpers::numCells(grid);
            swat_init.resize(nc);
            const int* gc = UgGridHelpers::globalCell(grid);
            for (int c = 0; c < nc; ++c) {
                const int deck_pos = (gc == NULL) ? c : gc[c];
                swat_init[c] = swat_init_ecl[deck_pos];
            }
        }

        ISC isc(materialLawManager, pu, deck, eclipseState, grid, gravity, swat_init);
        const int ref_phase = pu.phase_used[BlackoilPhases::Liquid]
            ? pu.phase_pos[BlackoilPhases::Liquid]
            : pu.phase_pos[BlackoilPhases::Aqua];
        state.pressure() = isc.press()[ref_phase];
        state.saturation() = Details::convertSats(isc.saturation());
        state.gasoilratio() = isc.rs();
        state.rv() = isc.rv();

        //initBlackoilSurfvolUsingRSorRV(UgGridHelpers::numCells(grid), props, state);
    }



} // namespace Opm

#endif // OPM_INITSTATEEQUIL_IMPL_HEADER_INCLUDED
