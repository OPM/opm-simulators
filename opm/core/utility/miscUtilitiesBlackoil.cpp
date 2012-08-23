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


#include <opm/core/utility/miscUtilitiesBlackoil.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid.h>
#include <opm/core/fluid/BlackoilPropertiesInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iterator>

namespace Opm
{

    /// @brief Computes injected and produced volumes of all phases.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// @param[in]  props     fluid and rock properties.
    /// @param[in]  p         pressure (one value per cell)
    /// @param[in]  z         surface-volume values (for all P phases)
    /// @param[in]  s         saturation values (for all P phases)
    /// @param[in]  src       if < 0: total outflow, if > 0: first phase inflow.
    /// @param[in]  dt        timestep used
    /// @param[out] injected  must point to a valid array with P elements,
    ///                       where P = s.size()/src.size().
    /// @param[out] produced  must also point to a valid array with P elements.
    void computeInjectedProduced(const BlackoilPropertiesInterface& props,
                                 const std::vector<double>& press,
                                 const std::vector<double>& z,
                                 const std::vector<double>& s,
                                 const std::vector<double>& src,
                                 const double dt,
                                 double* injected,
                                 double* produced)
    {
        const int num_cells = src.size();
        const int np = s.size()/src.size();
        if (int(s.size()) != num_cells*np) {
            THROW("Sizes of s and src vectors do not match.");
        }
        std::fill(injected, injected + np, 0.0);
        std::fill(produced, produced + np, 0.0);
        std::vector<double> visc(np);
        std::vector<double> mob(np);
        for (int c = 0; c < num_cells; ++c) {
            if (src[c] > 0.0) {
                injected[0] += src[c]*dt;
            } else if (src[c] < 0.0) {
                const double flux = -src[c]*dt;
                const double* sat = &s[np*c];
                props.relperm(1, sat, &c, &mob[0], 0);
                props.viscosity(1, &press[c], &z[np*c], &c, &visc[0], 0);
                double totmob = 0.0;
                for (int p = 0; p < np; ++p) {
                    mob[p] /= visc[p];
                    totmob += mob[p];
                }
                for (int p = 0; p < np; ++p) {
                    produced[p] += (mob[p]/totmob)*flux;
                }
            }
        }
    }



    /// @brief Computes total mobility for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  p         pressure (one value per cell)
    /// @param[in]  z         surface-volume values (for all P phases)
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobilities.
    void computeTotalMobility(const Opm::BlackoilPropertiesInterface& props,
                              const std::vector<int>& cells,
                              const std::vector<double>& press,
                              const std::vector<double>& z,
                              const std::vector<double>& s,
                              std::vector<double>& totmob)
    {
        std::vector<double> pmobc;

        computePhaseMobilities(props, cells, press, z, s, pmobc);

        const std::size_t                 np = props.numPhases();
        const std::vector<int>::size_type nc = cells.size();

        totmob.clear();
        totmob.resize(nc, 0.0);

        for (std::vector<int>::size_type c = 0; c < nc; ++c) {
            for (std::size_t p = 0; p < np; ++p) {
                totmob[ c ] += pmobc[c*np + p];
            }
        }
    }

    /*
    /// @brief Computes total mobility and omega for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  p         pressure (one value per cell)
    /// @param[in]  z         surface-volume values (for all P phases)
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobility
    /// @param[out] omega     fractional-flow weighted fluid densities.
    void computeTotalMobilityOmega(const Opm::BlackoilPropertiesInterface& props,
                                   const std::vector<int>& cells,
                                   const std::vector<double>& p,
                                   const std::vector<double>& z,
                                   const std::vector<double>& s,
                                   std::vector<double>& totmob,
                                   std::vector<double>& omega)
    {
        std::vector<double> pmobc;

        computePhaseMobilities(props, cells, p, z, s, pmobc);

        const std::size_t                 np = props.numPhases();
        const std::vector<int>::size_type nc = cells.size();

        totmob.clear();
        totmob.resize(nc, 0.0);
        omega.clear();
        omega.resize(nc, 0.0);

        const double* rho = props.density();
        for (std::vector<int>::size_type c = 0; c < nc; ++c) {
            for (std::size_t p = 0; p < np; ++p) {
                totmob[ c ] += pmobc[c*np + p];
                omega [ c ] += pmobc[c*np + p] * rho[ p ];
            }

            omega[ c ] /= totmob[ c ];
        }
    }
    */

    /// @brief Computes phase mobilities for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  p         pressure (one value per cell)
    /// @param[in]  z         surface-volume values (for all P phases)
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] pmobc     phase mobilities (for all phases).
    void computePhaseMobilities(const Opm::BlackoilPropertiesInterface& props,
                                const std::vector<int>&                 cells,
                                const std::vector<double>&              p,
                                const std::vector<double>&              z,
                                const std::vector<double>&              s,
                                std::vector<double>&                    pmobc)
    {
        const int nc = props.numCells();
        const int np = props.numPhases();

        ASSERT (int(s.size()) == nc * np);

        std::vector<double> mu(nc*np);
        props.viscosity(nc, &p[0], &z[0], &cells[0], &mu[0], 0);

        pmobc.clear();
        pmobc.resize(nc*np, 0.0);
        double* dpmobc = 0;
        props.relperm(nc, &s[0], &cells[0],
                      &pmobc[0], dpmobc);

        std::transform(pmobc.begin(), pmobc.end(),
                       mu.begin(),
                       pmobc.begin(),
                       std::divides<double>());
    }

    /// Computes the fractional flow for each cell in the cells argument
    /// @param[in]  props            rock and fluid properties
    /// @param[in]  cells            cells with which the saturation values are associated
    /// @param[in]  p                pressure (one value per cell)
    /// @param[in]  z                surface-volume values (for all P phases)
    /// @param[in]  s                saturation values (for all phases)
    /// @param[out] fractional_flow  the fractional flow for each phase for each cell.
    void computeFractionalFlow(const Opm::BlackoilPropertiesInterface& props,
                               const std::vector<int>& cells,
                               const std::vector<double>& p,
                               const std::vector<double>& z,
                               const std::vector<double>& s,
                               std::vector<double>& fractional_flows)
    {
        const int num_phases = props.numPhases();

        computePhaseMobilities(props, cells, p, z, s, fractional_flows);

        for (std::vector<int>::size_type i = 0; i < cells.size(); ++i) {
            double phase_sum = 0.0;
            for (int phase = 0; phase < num_phases; ++phase) {
                phase_sum += fractional_flows[i * num_phases + phase];
            }
            for (int phase = 0; phase < num_phases; ++phase) {
                fractional_flows[i * num_phases + phase] /= phase_sum;
            }
        }
    }

    /// Computes the surface volume densities from saturations by the formula
    ///     z = A s
    /// for a number of data points, where z is the surface volume density,
    /// s is the saturation (both as column vectors) and A is the
    /// phase-to-component relation matrix.
    /// @param[in]  n            number of data points
    /// @param[in]  np           number of phases, must be 2 or 3
    /// @param[in]  A            array containing n square matrices of size num_phases^2,
    ///                          in Fortran ordering, typically the output of a call
    ///                          to the matrix() method of a BlackoilProperties* class.
    /// @param[in]  saturation   concatenated saturation values (for all P phases)
    /// @param[out] surfacevol   concatenated surface-volume values (for all P phases)
    void computeSurfacevol(const int n,
                           const int np,
                           const double* A,
                           const double* saturation,
                           double* surfacevol)
    {
        // Note: since this is a simple matrix-vector product, it can
        // be done by a BLAS call, but then we have to reorder the A
        // matrix data.
        std::fill(surfacevol, surfacevol + n*np, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int row = 0; row < np; ++row) {
                for (int col = 0; col < np; ++col) {
                    surfacevol[i*np + row] += A[i*np*np + row*np + col] * saturation[i*np + col];
                }
            }
        }
    }

} // namespace Opm
