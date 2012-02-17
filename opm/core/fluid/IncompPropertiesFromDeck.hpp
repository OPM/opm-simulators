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

#ifndef OPM_INCOMPPROPERTIESFROMDECK_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESFROMDECK_HEADER_INCLUDED

#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/fluid/RockFromDeck.hpp>
#include <opm/core/fluid/PvtPropertiesIncompFromDeck.hpp>
#include <opm/core/fluid/SaturationPropsFromDeck.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>

namespace Opm
{

    /// Concrete class implementing the incompressible property
    /// interface, reading all data and properties from eclipse deck
    /// input.
    ///
    /// Supports variable number of spatial dimensions, called D.
    /// Supports variable number of phases, called P.
    /// In general, when arguments call for n values of some vector or
    /// matrix property, such as saturation, they shall always be
    /// ordered cellwise:
    ///   [s^1_0 s^2_0 s^3_0 s^1_1 s^2_2 ... ]
    /// in which s^i_j denotes saturation of phase i in cell j.
    class IncompPropertiesFromDeck : public IncompPropertiesInterface
    {
    public:
        /// Construct from deck and cell mapping.
        /// \param  deck         eclipse input parser
        /// \param  global_cell  mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        IncompPropertiesFromDeck(const EclipseGridParser& deck,
				 const std::vector<int>& global_cell);

	/// Destructor.
        virtual ~IncompPropertiesFromDeck();

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

        /// \return Array of P density values.
        virtual const double* density() const;

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
        RockFromDeck rock_;
	PvtPropertiesIncompFromDeck pvt_;
        SaturationPropsFromDeck satprops_;
    };



} // namespace Opm


#endif // OPM_INCOMPPROPERTIESFROMDECK_HEADER_INCLUDED
