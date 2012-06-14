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

#ifndef OPM_INCOMPPROPERTIESDEFAULTPOLYMER_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESDEFAULTPOLYMER_HEADER_INCLUDED

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <vector>

namespace Opm
{

    class IncompPropertiesDefaultPolymer : public Opm::IncompPropertiesBasic
    {
    public:
        /// Construct from parameters.
	/// The following parameters are accepted (defaults):
	///    num_phases         (2)         Must be 1 or 2.
	///    relperm_func       ("Linear")  Must be "Constant", "Linear" or "Quadratic".
	///    rho1 [rho2, rho3]  (1.0e3)     Density in kg/m^3
	///    mu1 [mu2, mu3]     (1.0)       Viscosity in cP
	///    porosity           (1.0)       Porosity
	///    permeability       (100.0)     Permeability in mD
        IncompPropertiesDefaultPolymer(const Opm::parameter::ParameterGroup& param, int dim, int num_cells)
            : Opm::IncompPropertiesBasic(param, dim, num_cells)
        {
            ASSERT(numPhases() == 2);
            sw_.resize(3);
            sw_[0] = 0.2;
            sw_[1] = 0.7;
            sw_[2] = 1.0;
            krw_.resize(3);
            krw_[0] = 0.0;
            krw_[1] = 0.7;
            krw_[2] = 1.0;
            so_.resize(2);
            so_[0] = 0.3;
            so_[1] = 0.8;
            kro_.resize(2);
            kro_[0] = 0.0;
            kro_[1] = 1.0;
        }

        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)

        virtual void relperm(const int n,
                             const double* s,
                             const int* /*cells*/,
                             double* kr,
                             double* dkrds) const
        {
            // ASSERT(dkrds == 0);
            // We assume two phases flow
            for (int i = 0; i < n; ++i) {
                kr[2*i] = krw(s[2*i]);
                kr[2*i+1] = kro(s[2*i+1]);
                if (dkrds != 0) {
                    dkrds[4*i + 0] = krw_dsw(s[2*i]);
                    dkrds[4*i + 3] = kro_dso(s[2*i+1]);
                    dkrds[4*i + 1] = 0.0;
                    dkrds[4*i + 2] = 0.0;
                }
            }
        }

	/// Obtain the range of allowable saturation values.
	/// In cell cells[i], saturation of phase p is allowed to be
	/// in the interval [smin[i*P + p], smax[i*P + p]].
        /// \param[in]  n      Number of data points.
        /// \param[in]  cells  Array of n cell indices.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        virtual void satRange(const int n,
                              const int* /*cells*/,
                              double* smin,
                              double* smax) const
        {
            const int np = 2;
            for (int i = 0; i < n; ++i) {
                smin[np*i + 0] = sw_[0];
                smax[np*i + 0] = sw_.back();
                smin[np*i + 1] = 1.0 - sw_[0];
                smax[np*i + 1] = 1.0 - sw_.back();
            }
        }

    private:
        double krw(double s) const
        {
            return Opm::linearInterpolation(sw_, krw_, s);
        }

        double krw_dsw(double s) const
        {
            return Opm::linearInterpolationDerivative(sw_, krw_, s);
        }


        double kro(double s) const
        {
            return Opm::linearInterpolation(so_, kro_, s);
        }

        double kro_dso(double s) const
        {
            return Opm::linearInterpolationDerivative(so_, kro_, s);
        }

        std::vector<double> sw_;
        std::vector<double> krw_;
        std::vector<double> so_;
        std::vector<double> kro_;
    };

} // namespace Opm

#endif // OPM_INCOMPPROPERTIESDEFAULTPOLYMER_HEADER_INCLUDED
