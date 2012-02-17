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


#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{

    IncompPropertiesFromDeck::IncompPropertiesFromDeck(const EclipseGridParser& deck,
						       const std::vector<int>& global_cell)
    {
        rock_.init(deck, global_cell);
	pvt_.init(deck);
        satprops_.init(deck);
	if (pvt_.numPhases() != satprops_.numPhases()) {
	    THROW("IncompPropertiesFromDeck::IncompPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
		  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
	}
    }

    IncompPropertiesFromDeck::~IncompPropertiesFromDeck()
    {
    }


    /// \return   D, the number of spatial dimensions.
    int IncompPropertiesFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int IncompPropertiesFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* IncompPropertiesFromDeck::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* IncompPropertiesFromDeck::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int IncompPropertiesFromDeck::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return Array of P viscosity values.
    const double* IncompPropertiesFromDeck::viscosity() const
    {
	return pvt_.viscosity();
    }

    /// \return Array of P density values.
    const double* IncompPropertiesFromDeck::density() const
    {
	return pvt_.surfaceDensities();
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
    void IncompPropertiesFromDeck::relperm(const int n,
					   const double* s,
					   const int* /*cells*/,
					   double* kr,
					   double* dkrds) const
    {
        satprops_.relperm(n, s, kr, dkrds);
    }


    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)
    void IncompPropertiesFromDeck::capPress(const int n,
					    const double* s,
					    const int* /*cells*/,
					    double* pc,
					    double* dpcds) const
    {
        satprops_.capPress(n, s, pc, dpcds);
    }


    /// Obtain the range of allowable saturation values.
    /// In cell cells[i], saturation of phase p is allowed to be
    /// in the interval [smin[i*P + p], smax[i*P + p]].
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void IncompPropertiesFromDeck::satRange(const int n,
					    const int* /*cells*/,
					    double* smin,
					    double* smax) const
    {
	satprops_.satRange(n, smin, smax);
    }

} // namespace Opm

