/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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
#include <opm/core/props/BlackoilPhases.hpp>

#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

namespace Opm
{
    namespace Details {
        template <class RHS>
        class RK4IVP : public std::binary_function<double,double,double> {
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
            template <class Density>
            class Water {
            public:
                Water(const Density& rho,
                      const int      np,
                      const int      ix,
                      const double   norm_grav)
                    : rho_(rho)
                    , svol_(np, 0)
                    , ix_(ix)
                    , g_(norm_grav)
                {
                    svol_[ix_] = 1.0;
                }

                double
                operator()(const double /* depth */,
                           const double press) const
                {
                    return this->density(press) * g_;
                }

            private:
                const Density&      rho_;
                std::vector<double> svol_;
                const int           ix_;
                const double        g_;

                double
                density(const double press) const
                {
                    const std::vector<double>& rho = rho_(press, svol_);

                    return rho[ix_];
                }
            };

            template <class Density, class RS>
            class Oil {
            public:
                Oil(const Density& rho,
                    const RS&      rs,
                    const int      np,
                    const int      oix,
                    const int      gix,
                    const double   norm_grav)
                    : rho_(rho)
                    , rs_(rs)
                    , svol_(np, 0)
                    , oix_(oix)
                    , gix_(gix)
                    , g_(norm_grav)
                {
                    svol_[oix_] = 1.0;
                }

                double
                operator()(const double depth,
                           const double press) const
                {
                    return this->density(depth, press) * g_;
                }

            private:
                const Density&              rho_;
                const RS&                   rs_;
                mutable std::vector<double> svol_;
                const int                   oix_;
                const int                   gix_;
                const double                g_;

                double
                density(const double depth,
                        const double press) const
                {
                    if (gix_ >= 0) {
                        svol_[gix_] = rs_(depth, press);
                    }

                    const std::vector<double>& rho = rho_(press, svol_);
                    return rho[oix_];
                }
            };

            template <class Density, class RV>
            class Gas {
            public:
                Gas(const Density& rho,
                    const RV&      rv,
                    const int      np,
                    const int      gix,
                    const int      oix,
                    const double   norm_grav)
                    : rho_(rho)
                    , rv_(rv)
                    , svol_(np, 0)
                    , gix_(gix)
                    , oix_(oix)
                    , g_(norm_grav)
                {
                    svol_[gix_] = 1.0;
                }

                double
                operator()(const double depth,
                           const double press) const
                {
                    return this->density(depth, press) * g_;
                }

            private:
                const Density&              rho_;
                const RV&                   rv_;
                mutable std::vector<double> svol_;
                const int                   gix_;
                const int                   oix_;
                const double                g_;

                double
                density(const double depth,
                        const double press) const
                {
                    if (oix_ >= 0) {
                        svol_[oix_] = rv_(depth, press);
                    }

                    const std::vector<double>& rho = rho_(press, svol_);
                    return rho[gix_];
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
            template <class PressFunction,
                      class CellRange>
            void
            assign(const UnstructuredGrid&             G    ,
                   const std::array<PressFunction, 2>& f    ,
                   const double                        split,
                   const CellRange&                    cells,
                   std::vector<double>&                p    )
            {
                const int nd = G.dimensions;

                enum { up = 0, down = 1 };

                std::vector<double>::size_type c = 0;
                for (typename CellRange::const_iterator
                         ci = cells.begin(), ce = cells.end();
                     ci != ce; ++ci, ++c)
                {
                    assert (c < p.size());

                    const double z = G.cell_centroids[(*ci)*nd + (nd - 1)];
                    p[c] = (z < split) ? f[up](z) : f[down](z);
                }
            }

            template <class Region,
                      class CellRange>
            void
            water(const UnstructuredGrid&     G     ,
                  const Region&               reg   ,
                  const std::array<double,2>& span  ,
                  const double                grav  ,
                  const double                po_woc,
                  const CellRange&            cells ,
                  std::vector<double>&        press )
            {
                using PhasePressODE::Water;
                typedef Water<typename Region::CalcDensity> ODE;

                const PhaseUsage& pu = reg.phaseUsage();

                const int wix = PhaseIndex::water(pu);
                ODE drho(reg.densityCalculator(), pu.num_phases, wix, grav);

                const double z0 = reg.zwoc();
                const double p0 = po_woc - reg.pcow_woc(); // Pcow = Po - Pw

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> WPress;
                std::array<WPress,2> wpress = {
                    {
                        WPress(drho, up  , p0, 100)
                        ,
                        WPress(drho, down, p0, 100)
                    }
                };

                assign(G, wpress, z0, cells, press);
            }

            template <class Region,
                      class CellRange>
            void
            oil(const UnstructuredGrid&     G     ,
                const Region&               reg   ,
                const std::array<double,2>& span  ,
                const double                grav  ,
                const CellRange&            cells ,
                std::vector<double>&        press ,
                double&                     po_woc,
                double&                     po_goc)
            {
                using PhasePressODE::Oil;
                typedef Oil<typename Region::CalcDensity,
                            typename Region::CalcDissolution> ODE;

                const PhaseUsage& pu = reg.phaseUsage();

                const int oix = PhaseIndex::oil(pu);
                const int gix = PhaseIndex::gas(pu);
                ODE drho(reg.densityCalculator(),
                         reg.dissolutionCalculator(),
                         pu.num_phases, oix, gix, grav);

                const double z0 = reg.datum();
                const double p0 = reg.pressure();

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> OPress;
                std::array<OPress,2> opress = {
                    {
                        OPress(drho, up  , p0, 100)
                        ,
                        OPress(drho, down, p0, 100)
                    }
                };

                assign(G, opress, z0, cells, press);

                po_woc = -std::numeric_limits<double>::max();
                po_goc = -std::numeric_limits<double>::max();

                const double woc = reg.zwoc();
                // Compute Oil pressure at WOC
                if      (z0 > woc) { po_woc = opress[0](woc); } // WOC above datum
                else if (z0 < woc) { po_woc = opress[1](woc); } // WOC below datum
                else               { po_woc = p0;             } // WOC *at*  datum

                const double goc = reg.zgoc();
                // Compute Oil pressure at GOC
                if      (z0 > goc) { po_goc = opress[0](goc); } // GOC above datum
                else if (z0 < goc) { po_goc = opress[1](goc); } // GOC below datum
                else               { po_goc = p0;             } // GOC *at*  datum
            }

            template <class Region,
                      class CellRange>
            void
            gas(const UnstructuredGrid&     G     ,
                const Region&               reg   ,
                const std::array<double,2>& span  ,
                const double                grav  ,
                const double                po_goc,
                const CellRange&            cells ,
                std::vector<double>&        press )
            {
                using PhasePressODE::Gas;
                typedef Gas<typename Region::CalcDensity,
                            typename Region::CalcEvaporation> ODE;

                const PhaseUsage& pu = reg.phaseUsage();

                const int gix = PhaseIndex::gas(pu);
                const int oix = PhaseIndex::oil(pu);

                ODE drho(reg.densityCalculator(),
                         reg.evaporationCalculator(),
                         pu.num_phases, gix, oix, grav);

                const double z0 = reg.zgoc();
                const double p0 = po_goc + reg.pcgo_goc(); // Pcog = Pg - Po

                std::array<double,2> up   = {{ z0, span[0] }};
                std::array<double,2> down = {{ z0, span[1] }};

                typedef Details::RK4IVP<ODE> GPress;
                std::array<GPress,2> gpress = {
                    {
                        GPress(drho, up  , p0, 100)
                        ,
                        GPress(drho, down, p0, 100)
                    }
                };

                assign(G, gpress, z0, cells, press);
            }
        } // namespace PhasePressure

        template <class Region,
                  class CellRange>
        void
        equilibrateOWG(const UnstructuredGrid&             G,
                       const Region&                       reg,
                       const double                        grav,
                       const std::array<double,2>&         span,
                       const CellRange&                    cells,
                       std::vector< std::vector<double> >& press)
        {
            const PhaseUsage& pu = reg.phaseUsage();

            double po_woc = -1, po_goc = -1;
            if (PhaseUsed::oil(pu)) {
                const int oix = PhaseIndex::oil(pu);
                PhasePressure::oil(G, reg, span, grav, cells,
                                   press[ oix ], po_woc, po_goc);
            }

            if (PhaseUsed::water(pu)) {
                const int wix = PhaseIndex::water(pu);
                PhasePressure::water(G, reg, span, grav, po_woc,
                                     cells, press[ wix ]);
            }

            if (PhaseUsed::gas(pu)) {
                const int gix = PhaseIndex::gas(pu);
                PhasePressure::gas(G, reg, span, grav, po_goc,
                                   cells, press[ gix ]);
            }
        }
    } // namespace Details

    namespace equil {
        template <class Region,
                  class CellRange>
        std::vector< std::vector<double> >
        phasePressures(const UnstructuredGrid& G,
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
                assert (G.dimensions == 3);

                const int nd = G.dimensions;

                // Define short-name aliases to reduce visual clutter.
                const double* const nc  = & G.node_coordinates[0];

                const int*    const cfp = & G.cell_facepos[0];
                const int*    const cf  = & G.cell_faces[0];

                const int*    const fnp = & G.face_nodepos[0];
                const int*    const fn  = & G.face_nodes[0];

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
                for (typename CellRange::const_iterator
                         ci = cells.begin(), ce = cells.end();
                     ci != ce; ++ci, ++ncell)
                {
                    for (const int
                             *fi = & cf[ cfp[*ci + 0] ],
                             *fe = & cf[ cfp[*ci + 1] ];
                         fi != fe; ++fi)
                    {
                        for (const int
                                 *i = & fn[ fnp[*fi + 0] ],
                                 *e = & fn[ fnp[*fi + 1] ];
                             i != e; ++i)
                        {
                            const double z = nc[(*i)*nd + (nd - 1)];

                            if (z < span[0]) { span[0] = z; }
                            if (z > span[1]) { span[1] = z; }
                        }
                    }
                }
            }

            const int np = reg.phaseUsage().num_phases;

            typedef std::vector<double> pval;
            std::vector<pval> press(np, pval(ncell, 0.0));

            const double z0   = reg.datum();
            const double zwoc = reg.zwoc ();
            const double zgoc = reg.zgoc ();

            if (! ((zgoc > z0) || (z0 > zwoc))) {
                // Datum depth in oil zone  (zgoc <= z0 <= zwoc)
                Details::equilibrateOWG(G, reg, grav, span, cells, press);
            }

            return press;
        }
    } // namespace equil
} // namespace Opm

#endif // OPM_INITSTATEEQUIL_IMPL_HEADER_INCLUDED
