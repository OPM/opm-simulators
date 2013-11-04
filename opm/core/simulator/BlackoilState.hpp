/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILSTATE_HEADER_INCLUDED
#define OPM_BLACKOILSTATE_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <vector>
#include <cmath>

namespace Opm
{

    /// Simulator state for a blackoil simulator.
    class BlackoilState
    {
    public:

        void init(const UnstructuredGrid& g, const int num_phases)
        {
            num_phases_ = num_phases;
            press_.resize(g.number_of_cells, 0.0);
            fpress_.resize(g.number_of_faces, 0.0);
            flux_.resize(g.number_of_faces, 0.0);
            sat_.resize(num_phases * g.number_of_cells, 0.0);
            for (int cell = 0; cell < g.number_of_cells; ++cell) {
                // Defaulting the second saturation to 1.0.
                // This will usually be oil in a water-oil case,
                // gas in an oil-gas case.
                // For proper initialization, one should not rely on this,
                // but use available phase information instead.
                sat_[num_phases*cell + 1] = 1.0;
            }
            gor_.resize(g.number_of_cells, 0.0);
        }

        enum ExtremalSat { MinSat, MaxSat };

        /// Set the first saturation to either its min or max value in
        /// the indicated cells. The second saturation value s2 is set
        /// to (1.0 - s1) for each cell. Any further saturation values
        /// are unchanged.
        void setFirstSat(const std::vector<int>& cells,
                         const Opm::BlackoilPropertiesInterface& props,
                         ExtremalSat es)
        {
            if (cells.empty()) {
                return;
            }
            const int n = cells.size();
            assert(n > 0);
            std::vector<double> smin(num_phases_*n);
            std::vector<double> smax(num_phases_*n);
            props.satRange(n, &cells[0], &smin[0], &smax[0]);
            const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
            for (int ci = 0; ci < n; ++ci) {
                const int cell = cells[ci];
                sat_[num_phases_*cell] = svals[num_phases_*ci];
                sat_[num_phases_*cell + 1] = 1.0 - sat_[num_phases_*cell];
            }
        }

        int numPhases() const
        {
            return num_phases_;
        }



        bool equals(const BlackoilState& other, double epsilon = 1e-8) const {
            bool equal = (num_phases_ == other.num_phases_);

            equal = equal && (vectorApproxEqual( pressure() , other.pressure() , epsilon));
            equal = equal && (vectorApproxEqual( facepressure() , other.facepressure() , epsilon));
            equal = equal && (vectorApproxEqual( faceflux() , other.faceflux() , epsilon));
            equal = equal && (vectorApproxEqual( surfacevol() , other.surfacevol() , epsilon));
            equal = equal && (vectorApproxEqual( saturation() , other.saturation() , epsilon));
            equal = equal && (vectorApproxEqual( gasoilratio() , other.gasoilratio() , epsilon));

            return equal;
        }


        std::vector<double>& pressure    () { return press_ ; }
        std::vector<double>& facepressure() { return fpress_; }
        std::vector<double>& faceflux    () { return flux_  ; }
        std::vector<double>& surfacevol  () { return surfvol_; }
        std::vector<double>& saturation  () { return sat_   ; }
        std::vector<double>& gasoilratio () { return gor_   ; }

        const std::vector<double>& pressure    () const { return press_ ; }
        const std::vector<double>& facepressure() const { return fpress_; }
        const std::vector<double>& faceflux    () const { return flux_  ; }
        const std::vector<double>& surfacevol  () const { return surfvol_; }
        const std::vector<double>& saturation  () const { return sat_   ; }
        const std::vector<double>& gasoilratio () const { return gor_   ; }

    private:
        int num_phases_;
        std::vector<double> press_ ;
        std::vector<double> fpress_;
        std::vector<double> flux_  ;
        std::vector<double> surfvol_;
        std::vector<double> sat_   ;
        std::vector<double> gor_   ;


        static bool vectorApproxEqual(const std::vector<double>& v1, const std::vector<double>& v2 , double epsilon) {
            if (v1.size() != v2.size())
                return false;
            
            for (size_t i = 0; i < v1.size(); i++)
                if (std::abs(v1[i] - v2[i]) > epsilon * (std::abs(v1[i]) + std::abs(v2[i])))
                    return false;
            
            return true;
        }

    };

} // namespace Opm


#endif // OPM_BLACKOILSTATE_HEADER_INCLUDED
