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

#ifndef OPM_SATURATIONPROPSBASIC_HEADER_INCLUDED
#define OPM_SATURATIONPROPSBASIC_HEADER_INCLUDED

#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm
{


    /// Class encapsulating basic saturation function behaviour,
    /// by which we mean constant, linear or quadratic relative
    /// permeability functions for a maximum of two phases,
    /// and zero capillary pressure.
    ///
    /// TODO: This class can easily be extended to three phases,
    /// by adding three-phase relperm behaviour.
    class SaturationPropsBasic
    {
    public:
        /// Default constructor.
        SaturationPropsBasic();

        /// Initialize from parameters.
        /// The following parameters are accepted (defaults):
        ///    num_phases   (2)          Must be 1 or 2.
        ///    relperm_func ("Linear")   Must be "Constant", "Linear" or "Quadratic".
        void init(const parameter::ParameterGroup& param);

        enum RelPermFunc { Constant, Linear, Quadratic };

        /// Initialize from arguments a basic Saturation property.
        void init(const int num_phases,
                  const RelPermFunc& relperm_func)
        {
            num_phases_ = num_phases;
            relperm_func_ = relperm_func;
        }

        /// \return   P, the number of phases.
        int numPhases() const;

        /// Relative permeability.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        void relperm(const int n,
                     const double* s,
                     double* kr,
                     double* dkrds) const;

        /// Capillary pressure.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dpc_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        void capPress(const int n,
                      const double* s,
                      double* pc,
                      double* dpcds) const;

        /// Obtain the range of allowable saturation values.
        /// \param[in]  n      Number of data points.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        void satRange(const int n,
                      double* smin,
                      double* smax) const;


    private:
        int num_phases_;
        RelPermFunc relperm_func_;
    };



} // namespace Opm




#endif // OPM_SATURATIONPROPSBASIC_HEADER_INCLUDED
