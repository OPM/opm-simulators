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

#include <opm/core/fluid/BlackoilPropertiesFromDeck.hpp>

namespace Opm
{

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(const Dune::EclipseGridParser& deck)
    {
        pvt_.init(deck);
        satprops_.init(deck);
    }

    BlackoilPropertiesFromDeck::~BlackoilPropertiesFromDeck()
    {
    }


    /// \return   D, the number of spatial dimensions.
    int BlackoilPropertiesFromDeck::numDimensions() const
    {
        return 3;
    }

    /// \return   N, the number of cells.
    int BlackoilPropertiesFromDeck::numCells() const
    {
        return -1;
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropertiesFromDeck::porosity() const
    {
        return 0;
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropertiesFromDeck::permeability() const
    {
        return 0;
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int BlackoilPropertiesFromDeck::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  p      Array of n pressure values.
    /// \param[in]  z      Array of nP surface volume values.
    /// \param[in]  cells  Array of n cell indices to be associated with the p and z values.
    /// \param[out] mu     Array of nP viscosity values, array must be valid before calling.
    /// \param[out] dmudp  If non-null: array of nP viscosity derivative values,
    ///                    array must be valid before calling.
    void BlackoilPropertiesFromDeck::viscosity(const int n,
                                               const double* p,
                                               const double* z,
                                               const int* /*cells*/,
                                               double* mu,
                                               double* dmudp) const
    {
        if (dmudp) {
            THROW("BlackoilPropertiesFromDeck::viscosity()  --  derivatives of viscosity not yet implemented.");
        } else {
            pvt_.mu(n, p, z, mu);
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  p      Array of n pressure values.
    /// \param[in]  z      Array of nP surface volume values.
    /// \param[in]  cells  Array of n cell indices to be associated with the p and z values.
    /// \param[out] A      Array of nP^2 values, array must be valid before calling.
    ///                    The P^2 values for a cell give the matrix A = RB^{-1} which
    ///                    relates z to u by z = Au. The matrices are output in Fortran order.
    /// \param[out] dAdp   If non-null: array of nP^2 matrix derivative values,
    ///                    array must be valid before calling. The matrices are output
    ///                    in Fortran order.
    void BlackoilPropertiesFromDeck::matrix(const int n,
                                            const double* p,
                                            const double* z,
                                            const int* cells,
                                            double* A,
                                            double* dAdp) const
    {
        THROW("BlackoilPropertiesFromDeck::matrix() not yet implemented.");
    }


    /// \param[in]  n      Number of data points.
    /// \param[in]  A      Array of nP^2 values, where the P^2 values for a cell give the
    ///                    matrix A = RB^{-1} which relates z to u by z = Au. The matrices
    ///                    are assumed to be in Fortran order, and are typically the result
    ///                    of a call to the method matrix().
    /// \param[out] rho    Array of nP density values, array must be valid before calling.
    void BlackoilPropertiesFromDeck::density(const int n,
                                             const double* A,
                                             double* rho) const
    {
        int np = numPhases();
        const double* sdens = pvt_.surfaceDensities();
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            for (int phase = 0; phase < np; ++phase) {
                rho[np*i + phase] = 0.0;
                for (int comp = 0; comp < np; ++comp) {
                    rho[np*i + phase] += A[n*np*np + np*phase + comp]*sdens[comp];
                }
            }
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void BlackoilPropertiesFromDeck::relperm(const int n,
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
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void BlackoilPropertiesFromDeck::capPress(const int n,
                                              const double* s,
                                              const int* /*cells*/,
                                              double* pc,
                                              double* dpcds) const
    {
        satprops_.relperm(n, s, pc, dpcds);
    }



} // namespace Opm

