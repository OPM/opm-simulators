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

#include "config.h"
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm
{
    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr newParserDeck,
                                                           const UnstructuredGrid& grid,
                                                           bool init_rock)
    {
        init(newParserDeck, grid.number_of_cells, grid.global_cell, grid.cartdims,
             grid.cell_centroids, grid.dimensions, init_rock);
    }

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr newParserDeck,
                                                           const UnstructuredGrid& grid,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        init(newParserDeck, grid.number_of_cells, grid.global_cell, grid.cartdims, grid.cell_centroids, 
             grid.dimensions, param, init_rock);
    }

    BlackoilPropertiesFromDeck::~BlackoilPropertiesFromDeck()
    {
    }


    /// \return   D, the number of spatial dimensions.
    int BlackoilPropertiesFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int BlackoilPropertiesFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropertiesFromDeck::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropertiesFromDeck::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int BlackoilPropertiesFromDeck::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return   Object describing the active phases.
    PhaseUsage BlackoilPropertiesFromDeck::phaseUsage() const
    {
        return pvt_.phaseUsage();
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
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::viscosity()  --  derivatives of viscosity not yet implemented.");
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
                                            const int* /*cells*/,
                                            double* A,
                                            double* dAdp) const
    {
        const int np = numPhases();
        B_.resize(n*np);
        R_.resize(n*np);
        if (dAdp) {
            dB_.resize(n*np);
            dR_.resize(n*np);
            pvt_.dBdp(n, p, z, &B_[0], &dB_[0]);
            pvt_.dRdp(n, p, z, &R_[0], &dR_[0]);
        } else {
            pvt_.B(n, p, z, &B_[0]);
            pvt_.R(n, p, z, &R_[0]);
        }
        const int* phase_pos = pvt_.phasePosition();
        bool oil_and_gas = pvt_.phaseUsed()[BlackoilPhases::Liquid] &&
            pvt_.phaseUsed()[BlackoilPhases::Vapour];
        const int o = phase_pos[BlackoilPhases::Liquid];
        const int g = phase_pos[BlackoilPhases::Vapour];

        // Compute A matrix
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double* m = A + i*np*np;
            std::fill(m, m + np*np, 0.0);
            // Diagonal entries.
            for (int phase = 0; phase < np; ++phase) {
                m[phase + phase*np] = 1.0/B_[i*np + phase];
            }
            // Off-diagonal entries.
            if (oil_and_gas) {
                m[o + g*np] = R_[i*np + g]/B_[i*np + g];
                m[g + o*np] = R_[i*np + o]/B_[i*np + o];
            }
        }

        // Derivative of A matrix.
        // A     = R*inv(B) whence
        //
        // dA/dp = (dR/dp*inv(B) + R*d(inv(B))/dp)
        //       = (dR/dp*inv(B) - R*inv(B)*(dB/dp)*inv(B))
        //       = (dR/dp - A*(dB/dp)) * inv(B)
        //
        // The B matrix is diagonal and that fact is exploited in the
        // following implementation.
        if (dAdp) {
// #pragma omp parallel for
            // (1): dA/dp <- A
            std::copy(A, A + n*np*np, dAdp);

            for (int i = 0; i < n; ++i) {
                double*       m  = dAdp + i*np*np;

                // (2): dA/dp <- -dA/dp*(dB/dp) == -A*(dB/dp)
                const double* dB = & dB_[i * np];
                for (int col = 0; col < np; ++col) {
                    for (int row = 0; row < np; ++row) {
                        m[col*np + row] *= - dB[ col ]; // Note sign.
                    }
                }

                if (oil_and_gas) {
                    // (2b): dA/dp += dR/dp (== dR/dp - A*(dB/dp))
                    const double* dR = & dR_[i * np];

                    m[o*np + g] += dR[ o ];
                    m[g*np + o] += dR[ g ];
                }

                // (3): dA/dp *= inv(B) (== final result)
                const double* B = & B_[i * np];
                for (int col = 0; col < np; ++col) {
                    for (int row = 0; row < np; ++row) {
                        m[col*np + row] /= B[ col ];
                    }
                }
            }
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
        const int np = numPhases();
        const double* sdens = pvt_.surfaceDensities();
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            for (int phase = 0; phase < np; ++phase) {
                rho[np*i + phase] = 0.0;
                for (int comp = 0; comp < np; ++comp) {
                    rho[np*i + phase] += A[i*np*np + np*phase + comp]*sdens[comp];
                }
            }
        }
    }

    /// Densities of stock components at surface conditions.
    /// \return Array of P density values.
    const double* BlackoilPropertiesFromDeck::surfaceDensity() const
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
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void BlackoilPropertiesFromDeck::relperm(const int n,
                                             const double* s,
                                             const int* cells,
                                             double* kr,
                                             double* dkrds) const
    {
        satprops_->relperm(n, s, cells, kr, dkrds);
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
                                              double* pc,
                                              double* dpcds) const
    {
        satprops_->capPress(n, s, cells, pc, dpcds);
    }


    /// Obtain the range of allowable saturation values.
    /// In cell cells[i], saturation of phase p is allowed to be
    /// in the interval [smin[i*P + p], smax[i*P + p]].
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void BlackoilPropertiesFromDeck::satRange(const int n,
                                              const int* cells,
                                              double* smin,
                                              double* smax) const
    {
        satprops_->satRange(n, cells, smin, smax);
    }


} // namespace Opm

