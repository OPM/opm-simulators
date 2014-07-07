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
#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{

    BlackoilPropertiesBasic::BlackoilPropertiesBasic(const parameter::ParameterGroup& param,
                                                     const int dim,
                                                     const int num_cells)
    {
        double poro = param.getDefault("porosity", 1.0);
        using namespace Opm::unit;
        using namespace Opm::prefix;
        double perm = param.getDefault("permeability", 100.0)*milli*darcy;
        rock_.init(dim, num_cells, poro, perm);
        pvt_.init(param);
        satprops_.init(param);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesBasic::BlackoilPropertiesBasic() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
    }

    BlackoilPropertiesBasic::~BlackoilPropertiesBasic()
    {
    }


    /// \return   D, the number of spatial dimensions.
    int BlackoilPropertiesBasic::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int BlackoilPropertiesBasic::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropertiesBasic::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropertiesBasic::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int BlackoilPropertiesBasic::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return   Object describing the active phases.
    PhaseUsage BlackoilPropertiesBasic::phaseUsage() const
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
    void BlackoilPropertiesBasic::viscosity(const int n,
                                            const double* p,
                                            const double* z,
                                            const int* /*cells*/,
                                            double* mu,
                                            double* dmudp) const
    {
        if (dmudp) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesBasic::viscosity()  --  derivatives of viscosity not yet implemented.");
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
    void BlackoilPropertiesBasic::matrix(const int n,
                                         const double* /*p*/,
                                         const double* /*z*/,
                                         const int* /*cells*/,
                                         double* A,
                                         double* dAdp) const
    {
        const int np = numPhases();
        assert(np <= 2);
        double B[2]; // Must be enough since component classes do not handle more than 2.
        pvt_.B(1, 0, 0, B);
        // Compute A matrix
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double* m = A + i*np*np;
            std::fill(m, m + np*np, 0.0);
            // Diagonal entries only.
            for (int phase = 0; phase < np; ++phase) {
                m[phase + phase*np] = 1.0/B[phase];
            }
        }

        // Derivative of A matrix.
        if (dAdp) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                double* m = dAdp + i*np*np;
                std::fill(m, m + np*np, 0.0);
            }
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  A      Array of nP^2 values, where the P^2 values for a cell give the
    ///                    matrix A = RB^{-1} which relates z to u by z = Au. The matrices
    ///                    are assumed to be in Fortran order, and are typically the result
    ///                    of a call to the method matrix().
    /// \param[in]  cells  The index of the grid cell of each data point.
    /// \param[out] rho    Array of nP density values, array must be valid before calling.
    void BlackoilPropertiesBasic::density(const int n,
                                          const double* A,
                                          const int* /*cells*/,
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
    const double* BlackoilPropertiesBasic::surfaceDensity(int /*cellIdx*/) const
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
    void BlackoilPropertiesBasic::relperm(const int n,
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
    void BlackoilPropertiesBasic::capPress(const int n,
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
    void BlackoilPropertiesBasic::satRange(const int n,
                                           const int* /*cells*/,
                                           double* smin,
                                           double* smax) const
    {
        satprops_.satRange(n, smin, smax);
    }


    /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
    /// \param[in]     cell   Cell index. 
    /// \param[in]     pcow   P_oil - P_water.
    /// \param[in/out] swat   Water saturation. / Possibly modified Water saturation.      
    void BlackoilPropertiesBasic::swatInitScaling(const int cell, 
                                                  const double pcow, 
                                                  double & swat)
    {
        OPM_THROW(std::runtime_error, "BlackoilPropertiesBasic::swatInitScaling()  --  not implemented.");
    }

} // namespace Opm

