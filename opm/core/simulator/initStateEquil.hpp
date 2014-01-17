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

#ifndef OPM_INITSTATEEQUIL_HEADER_INCLUDED
#define OPM_INITSTATEEQUIL_HEADER_INCLUDED

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/core/utility/Units.hpp>

#include <array>
#include <cassert>
#include <utility>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{
    namespace equil {
        template <class Props>
        class DensityCalculator;

        template <>
        class DensityCalculator< BlackoilPropertiesInterface > {
        public:
            DensityCalculator(const BlackoilPropertiesInterface& props,
                              const int                          c)
                : props_(props)
                , c_(1, c)
            {
            }

            std::vector<double>
            operator()(const double               p,
                       const std::vector<double>& z) const
            {
                const int np = props_.numPhases();
                std::vector<double> A(np * np, 0);

                assert (z.size() == std::vector<double>::size_type(np));

                double* dAdp = 0;
                props_.matrix(1, &p, &z[0], &c_[0], &A[0], dAdp);

                std::vector<double> rho(np, 0.0);
                props_.density(1, &A[0], &rho[0]);

                return rho;
            }

        private:
            const BlackoilPropertiesInterface& props_;
            const std::vector<int>             c_;
        };

        namespace miscibility {
            struct NoMixing {
                double
                operator()(const double /* depth */,
                           const double /* press */) const
                {
                    return 0.0;
                }
            };

            class RsVD {
            public:
                RsVD(const std::vector<double>& depth,
                     const std::vector<double>& rs)
                    : depth_(depth)
                    , rs_(rs)
                {
                }

                double
                operator()(const double depth,
                           const double /* press */) const
                {
                    return linearInterpolation(depth_, rs_, depth);
                }

            private:
                std::vector<double> depth_;
                std::vector<double> rs_;
            };
        } // namespace miscibility

        template < class Region = std::vector<int> >
        class RegionMapping {
        public:
            explicit
            RegionMapping(const Region& reg)
                : reg_(reg)
            {
                rev_.init(reg_);
            }

            typedef typename Region::value_type                  RegionId;
            typedef typename Region::size_type                   CellId;
            typedef typename std::vector<CellId>::const_iterator CellIter;

            class CellRange {
            public:
                CellRange(const CellIter b,
                          const CellIter e)
                    : b_(b), e_(e)
                {}

                typedef CellIter iterator;
                typedef CellIter const_iterator;

                iterator begin() const { return b_; }
                iterator end()   const { return e_; }

            private:
                iterator b_;
                iterator e_;
            };

            RegionId
            numRegions() const { return RegionId(rev_.p.size()) - 1; }

            RegionId
            region(const CellId c) const { return reg_[c]; }

            CellRange
            cells(const RegionId r) const {
                const RegionId i = r - rev_.low;
                return CellRange(rev_.c.begin() + rev_.p[i + 0],
                                 rev_.c.begin() + rev_.p[i + 1]);
            }

        private:
            Region reg_;

            struct {
                typedef typename std::vector<CellId>::size_type Pos;
                std::vector<Pos>    p;
                std::vector<CellId> c;
                RegionId            low;

                void
                init(const Region& reg)
                {
                    typedef typename Region::const_iterator CI;
                    const std::pair<CI,CI>
                        m = std::minmax_element(reg.begin(), reg.end());

                    low  = *m.first;

                    const typename Region::size_type
                        n = *m.second - low + 1;

                    p.resize(n + 1);  std::fill(p.begin(), p.end(), Pos(0));
                    for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                        p[ reg[i] - low + 1 ] += 1;
                    }

                    for (typename std::vector<Pos>::size_type
                             i = 1, sz = p.size(); i < sz; ++i) {
                        p[0] += p[i];
                        p[i]  = p[0] - p[i];
                    }

                    assert (p[0] ==
                            static_cast<typename Region::size_type>(reg.size()));

                    c.resize(reg.size());
                    for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                        c[ p[ reg[i] - low + 1 ] ++ ] = i;
                    }

                    p[0] = 0;
                }
            } rev_;
        };

        struct EquilRecord {
            struct {
                double depth;
                double press;
            } main, woc, goc;
        };

        template <class DensCalc,
                  class RS = miscibility::NoMixing,
                  class RV = miscibility::NoMixing>
        class EquilReg {
        public:
            EquilReg(const EquilRecord& rec,
                     const DensCalc&    density,
                     const RS&          rs,
                     const RV&          rv,
                     const PhaseUsage&  pu)
                : rec_    (rec)
                , density_(density)
                , rs_     (rs)
                , rv_     (rv)
                , pu_     (pu)
            {
            }

            typedef DensCalc CalcDensity;
            typedef RS       CalcDissolution;
            typedef RV       CalcEvaporation;

            double datum()    const { return this->rec_.main.depth; }
            double pressure() const { return this->rec_.main.press; }

            double zwoc()     const { return this->rec_.woc .depth; }
            double pcow_woc() const { return this->rec_.woc .press; }

            double zgoc()     const { return this->rec_.goc .depth; }
            double pcgo_goc() const { return this->rec_.goc .press; }

            const CalcDensity&
            densityCalculator() const { return this->density_; }

            const CalcDissolution&
            dissolutionCalculator() const { return this->rs_; }

            const CalcEvaporation&
            evaporationCalculator() const { return this->rv_; }

            const PhaseUsage&
            phaseUsage() const { return this->pu_; }

        private:
            EquilRecord rec_;
            DensCalc    density_;
            RS          rs_;
            RV          rv_;
            PhaseUsage  pu_;
        };

        template <class Region, class CellRange>
        std::vector< std::vector<double> >
        phasePressures(const UnstructuredGrid& G,
                       const Region&           reg,
                       const CellRange&        cells,
                       const double            grav = unit::gravity);
    } // namespace equil
} // namespace Opm

#include <opm/core/simulator/initStateEquil_impl.hpp>

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
