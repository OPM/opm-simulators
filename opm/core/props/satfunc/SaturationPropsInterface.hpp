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

#ifndef OPM_SATURATIONPROPSINTERFACE_HEADER_INCLUDED
#define OPM_SATURATIONPROPSINTERFACE_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>


namespace Opm
{

    class SaturationPropsInterface : public BlackoilPhases
    {
    public:
        /// Virtual destructor.
        virtual ~SaturationPropsInterface() {};

        /// \return   P, the number of phases.
        virtual int numPhases() const = 0;

        /// Relative permeability.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        virtual void relperm(const int n,
                             const double* s,
                             const int* cells,
                             double* kr,
                             double* dkrds) const = 0;

        /// Capillary pressure.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dpc_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        virtual void capPress(const int n,
                              const double* s,
                              const int* cells,
                              double* pc,
                              double* dpcds) const = 0;

        /// Obtain the range of allowable saturation values.
        /// \param[in]  n      Number of data points.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        virtual void satRange(const int n,
                              const int* cells,
                              double* smin,
                              double* smax) const = 0;
                                           
        /// Update saturation state for the hysteresis tracking 
        /// \param[in]  n      Number of data points. 
        /// \param[in]  s      Array of nP saturation values.             
        virtual void updateSatHyst(const int n,
                                   const int* cells,
                                   const double* s) = 0;
 
        /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
        /// \param[in]     cell  Cell index.
        /// \param[in]     pcow  P_oil - P_water.
        /// \param[in/out] swat  Water saturation. / Possibly modified Water saturation.      
        virtual void swatInitScaling(const int cell, 
                                     const double pcow, 
                                     double & swat) = 0;

    };



} // namespace Opm



#endif // OPM_SATURATIONPROPSINTERFACE_HEADER_INCLUDED
