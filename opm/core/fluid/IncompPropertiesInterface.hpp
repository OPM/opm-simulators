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

#ifndef OPM_INCOMPPROPERTIESINTERFACE_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESINTERFACE_HEADER_INCLUDED

namespace Opm
{

    /// Abstract base class for incompressible fluid and reservoir properties.
    ///
    /// Supports variable number of spatial dimensions, called D.
    /// Supports variable number of phases, called P.
    /// In general, when arguments call for n values of some vector or
    /// matrix property, such as saturation, they shall always be
    /// ordered cellwise:
    ///   [s^1_0 s^2_0 s^3_0 s^1_1 s^2_2 ... ]
    /// in which s^i_j denotes saturation of phase i in cell j.
    class IncompPropertiesInterface
    {
    public:
        virtual ~IncompPropertiesInterface() {}

        // ---- Rock interface ----

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const = 0;

        /// \return   N, the number of cells.
        virtual int numCells() const = 0;

        /// \return   Array of N porosity values.
        virtual const double* porosity() const = 0;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const = 0;


        // ---- Fluid interface ----

        /// \return   P, the number of phases (also the number of components).
        virtual int numPhases() const = 0;

        /// \return Array of P viscosity values.
        virtual const double* viscosity() const = 0;

        /// Densities of fluid phases at surface conditions.
        /// \return Array of P density values.
        virtual const double* density() const = 0;

        /// Densities of fluid phases at surface conditions.
        /// Note: a reasonable question to ask is why there can be
        /// different densities at surface and reservoir conditions,
        /// when the phases are assumed incompressible. The answer is
        /// that even if we approximate the phases as being
        /// incompressible during simulation, the density difference
        /// between surface and reservoir may be larger. For accurate
        /// reporting and using data given in terms of surface values,
        /// we need to handle this difference.
        /// \return Array of P density values.
        virtual const double* surfaceDensity() const = 0;

        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
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


        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
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
	/// In cell cells[i], saturation of phase p is allowed to be
	/// in the interval [smin[i*P + p], smax[i*P + p]].
        /// \param[in]  n      Number of data points.
        /// \param[in]  cells  Array of n cell indices.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        virtual void satRange(const int n,
                              const int* cells,
                              double* smin,
                              double* smax) const = 0;
    };



} // namespace Opm


#endif // OPM_INCOMPPROPERTIESINTERFACE_HEADER_INCLUDED
