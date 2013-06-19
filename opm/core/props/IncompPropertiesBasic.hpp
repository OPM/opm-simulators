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

#ifndef OPM_INCOMPPROPERTIESBASIC_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESBASIC_HEADER_INCLUDED

#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/props/rock/RockBasic.hpp>
#include <opm/core/props/pvt/PvtPropertiesBasic.hpp>
#include <opm/core/props/satfunc/SaturationPropsBasic.hpp>

namespace Opm
{

    /// Concrete class implementing the incompressible property
    /// interface, reading all necessary input from parameters.
    ///
    /// Supports variable number of spatial dimensions, called D.
    /// Supports variable number of phases, called P.
    /// In general, when arguments call for n values of some vector or
    /// matrix property, such as saturation, they shall always be
    /// ordered cellwise:
    ///   \f[ [s^1_0, s^2_0, s^3_0, s^1_1, s^2_2, \ldots ] \f]
    /// in which \f$ s^i_j \f$ denotes saturation of phase i in cell j.
    class IncompPropertiesBasic : public IncompPropertiesInterface
    {
    public:
        /// Construct from parameters.
        /// Note that all values passed through param should be in convenient units,
        /// as documented below.
        /// The following parameters are accepted (defaults):
        ///   - \c num_phases            (2)        -- Must be 1 or 2.
        ///   - \c relperm_func          ("Linear") -- Must be "Constant", "Linear" or "Quadratic".
        ///   - \c rho1 \c rho2, \c rho3 (1.0e3)    -- Density in kg/m^3.
        ///   - \c mu1 \c mu2, \c mu3    (1.0)      -- Viscosity in cP.
        ///   - \c porosity              (1.0)      -- Porosity.
        ///   - \c permeability          (100.0)    -- Permeability in mD.
        IncompPropertiesBasic(const parameter::ParameterGroup& param,
                              const int dim,
                              const int num_cells);


        /// Construct properties from arguments.
        /// Note that all values should be given in SI units
        /// for this constructor.
        /// \param[in] num_phases     Must be 1 or 2.
        /// \param[in] rho            Phase densities in kg/m^3.
        /// \param[in] mu             Phase viscosities in Pa*s.
        /// \param[in] porosity       Must be in [0,1].
        /// \param[in] permeability   Permeability in m^2.
        /// \param[in] dim            Must be 2 or 3.
        /// \param[in] num_cells      The number of grid cells.
        IncompPropertiesBasic(const int num_phases,
                              const SaturationPropsBasic::RelPermFunc& relpermfunc,
                              const std::vector<double>& rho,
                              const std::vector<double>& mu,
                              const double porosity,
                              const double permeability,
                              const int dim,
                              const int num_cells);

        /// Destructor.
        virtual ~IncompPropertiesBasic();

        // ---- Rock interface ----

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const;

        /// \return   N, the number of cells.
        virtual int numCells() const;

        /// \return   Array of N porosity values.
        virtual const double* porosity() const;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const;


        // ---- Fluid interface ----

        /// \return   P, the number of phases (also the number of components).
        virtual int numPhases() const;

        /// \return Array of P viscosity values.
        virtual const double* viscosity() const;

        /// Densities of fluid phases at reservoir conditions.
        /// \return Array of P density values.
        virtual const double* density() const;

        /// Densities of fluid phases at surface conditions.
        /// \return Array of P density values.
        virtual const double* surfaceDensity() const;

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
                             const int* cells,
                             double* kr,
                             double* dkrds) const;


        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dpc_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)
        virtual void capPress(const int n,
                              const double* s,
                              const int* cells,
                              double* pc,
                              double* dpcds) const;


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
                              double* smax) const;
    private:
        RockBasic rock_;
        PvtPropertiesBasic pvt_;
        SaturationPropsBasic satprops_;
        std::vector<double> viscosity_;
    };



} // namespace Opm


#endif // OPM_INCOMPPROPERTIESBASIC_HEADER_INCLUDED
