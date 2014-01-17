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
