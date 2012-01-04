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
        fluid_.init(deck);
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
        return 3;
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
                                               const int* cells,
                                               double* mu,
                                               double* dmudp) const
    {
        state_.phase_pressure.resize(n);
        state_.surface_volume_density.resize(n);
        int num_phases = numPhases();
        assert(num_phases == BlackoilFluid::numPhases);
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            state_.phase_pressure[i] = p[i];
            for (int phase = 0; phase < num_phases; ++phase) {
                state_.surface_volume_density[i][phase] = z[num_phases*i + phase];
            }
        }
        if (dmudp) {
            THROW("Sorry, derivatives of viscosity not yet done.");
        } else {
            fluid_.computePvtNoDerivs(state_); // Unnecessarily computes B and R
            const double* beg_mu = &(state_.viscosity[0][0]);
            std::copy(beg_mu, beg_mu + n*num_phases, mu);
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
        state_.phase_pressure.resize(n);
        state_.surface_volume_density.resize(n);
        int num_phases = numPhases();
        assert(num_phases == BlackoilFluid::numPhases);
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            state_.phase_pressure[i] = p[i];
            for (int phase = 0; phase < num_phases; ++phase) {
                state_.surface_volume_density[i][phase] = z[num_phases*i + phase];
            }
        }
        if (dAdp) {
            THROW("Sorry, derivatives of A matrix not yet done.");
        } else {
            fluid_.computeBAndR(state_);
            fluid_.computeStateMatrix(state_);
            const double* beg_A = &(state_.state_matrix[0][0][0]);
            std::copy(beg_A, beg_A + n*num_phases*num_phases, A);
        }
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
        int num_phases = numPhases();
        assert(num_phases == BlackoilFluid::numPhases);
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            BlackoilFluid::PhaseVec dens = fluid_.phaseDensities(A + n*num_phases*num_phases);
            for (int phase = 0; phase < num_phases; ++phase) {
                rho[num_phases*i + phase] = dens[phase];
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
                                             const int* cells,
                                             double* kr,
                                             double* dkrds) const
    {
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
                                              const int* cells,
                                              double* pv,
                                              double* dpcds) const
    {
    }



} // namespace Opm

