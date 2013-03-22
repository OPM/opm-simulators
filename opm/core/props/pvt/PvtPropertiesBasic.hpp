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

#ifndef OPM_PVTPROPERTIESBASIC_HEADER_INCLUDED
#define OPM_PVTPROPERTIESBASIC_HEADER_INCLUDED

#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm
{

    /// Class collecting simple pvt properties for 1-3 phases.
    /// All phases are incompressible and have constant viscosities.
    /// For all the methods, the following apply: p and z are unused.
    /// Output arrays shall be of size n*numPhases(), and must be valid
    /// before calling the method.
    /// NOTE: This class is intentionally similar to BlackoilPvtProperties.
    class PvtPropertiesBasic
    {
    public:
        /// Default constructor.
        PvtPropertiesBasic();

        /// Initialize from parameters.
        /// The following parameters are accepted (defaults):
        ///   - num_phases        (2)      --  Must be 1, 2 or 3.
        ///   - rho1, rho2, rho3  (1.0e3)  --  Density in kg/m^3
        ///   - mu1, mu2, mu3     (1.0)    --  Viscosity in cP
        void init(const parameter::ParameterGroup& param);

        /// Initialize from arguments.
        /// Basic multi phase fluid pvt properties.
        void init(const int num_phases,
                  const std::vector<double>& rho,
                  const std::vector<double>& visc);

        /// Number of active phases.
        int numPhases() const;

        /// Densities of stock components at surface conditions.
        /// \return  Array of size numPhases().
        const double* surfaceDensities() const;

        /// Viscosity as a function of p and z.
        void mu(const int n,
                const double* p,
                const double* z,
                double* output_mu) const;

        /// Formation volume factor as a function of p and z.
        void B(const int n,
               const double* p,
               const double* z,
               double* output_B) const;

        /// Formation volume factor and p-derivative as functions of p and z.
        void dBdp(const int n,
                  const double* p,
                  const double* z,
                  double* output_B,
                  double* output_dBdp) const;

        /// Solution factor as a function of p and z.
        void R(const int n,
               const double* p,
               const double* z,
               double* output_R) const;

        /// Solution factor and p-derivative as functions of p and z.
        void dRdp(const int n,
                  const double* p,
                  const double* z,
                  double* output_R,
                  double* output_dRdp) const;

    private:
        std::vector<double> density_;
        std::vector<double> viscosity_;
        std::vector<double> formation_volume_factor_;
    };

}



#endif // OPM_PVTPROPERTIESBASIC_HEADER_INCLUDED
