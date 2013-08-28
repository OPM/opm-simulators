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
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iterator>

namespace Opm
{

    /// @brief Computes injected and produced surface volumes of all phases.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// Note 3: Gives surface volume values, not reservoir volumes
    ///         (as the incompressible version of the function does).
    ///         Also, assumes that transport_src is given in surface volumes
    ///         for injector terms!
    /// @param[in]  props           fluid and rock properties.
    /// @param[in]  state           state variables (pressure, sat, surfvol)
    /// @param[in]  transport_src   if < 0: total resv outflow, if > 0: first phase surfv inflow
    /// @param[in]  dt              timestep used
    /// @param[out] injected        must point to a valid array with P elements,
    ///                             where P = s.size()/src.size().
    /// @param[out] produced        must also point to a valid array with P elements.
    void computeInjectedProduced(const BlackoilPropertiesInterface& props,
                                 const BlackoilState& state,
                                 const std::vector<double>& transport_src,
                                 const double dt,
                                 double* injected,
                                 double* produced)
    {
        const int num_cells = transport_src.size();
        if (props.numCells() != num_cells) {
            OPM_THROW(std::runtime_error, "Size of transport_src vector does not match number of cells in props.");
        }
        const int np = props.numPhases();
        if (int(state.saturation().size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of state vectors do not match number of cells.");
        }
        const std::vector<double>& press = state.pressure();
        const std::vector<double>& s = state.saturation();
        const std::vector<double>& z = state.surfacevol();
        std::fill(injected, injected + np, 0.0);
        std::fill(produced, produced + np, 0.0);
        std::vector<double> visc(np);
        std::vector<double> mob(np);
        std::vector<double> A(np*np);
        std::vector<double> prod_resv_phase(np);
        std::vector<double> prod_surfvol(np);
        for (int c = 0; c < num_cells; ++c) {
            if (transport_src[c] > 0.0) {
                // Inflowing transport source is a surface volume flux
                // for the first phase.
                injected[0] += transport_src[c]*dt;
            } else if (transport_src[c] < 0.0) {
                // Outflowing transport source is a total reservoir
                // volume flux.
                const double flux = -transport_src[c]*dt;
                const double* sat = &s[np*c];
                props.relperm(1, sat, &c, &mob[0], 0);
                props.viscosity(1, &press[c], &z[np*c], &c, &visc[0], 0);
                props.matrix(1, &press[c], &z[np*c], &c, &A[0], 0);
                double totmob = 0.0;
                for (int p = 0; p < np; ++p) {
                    mob[p] /= visc[p];
                    totmob += mob[p];
                }
                std::fill(prod_surfvol.begin(), prod_surfvol.end(), 0.0);
                for (int p = 0; p < np; ++p) {
                    prod_resv_phase[p] = (mob[p]/totmob)*flux;
                    for (int q = 0; q < np; ++q) {
                        prod_surfvol[q] += prod_resv_phase[p]*A[q + np*p];
                    }
                }
                for (int p = 0; p < np; ++p) {
                    produced[p] += prod_surfvol[p];
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
            for (int col = 0; col < np; ++col) {
                for (int row = 0; row < np; ++row) {
                    surfacevol[i*np + row] += A[i*np*np + row + col*np] * saturation[i*np + col];
                }
            }
        }
    }


    /// Compute two-phase transport source terms from well terms.
    /// Note: Unlike the incompressible version of this function,
    ///       this version computes surface volume injection rates,
    ///       production rates are still total reservoir volumes.
    /// \param[in]  props         Fluid and rock properties.
    /// \param[in]  wells         Wells data structure.
    /// \param[in]  well_state    Well pressures and fluxes.
    /// \param[out] transport_src The transport source terms. They are to be interpreted depending on sign:
    ///                           (+) positive  inflow of first (water) phase (surface volume),
    ///                           (-) negative  total outflow of both phases (reservoir volume).
    void computeTransportSource(const BlackoilPropertiesInterface& props,
                                const Wells* wells,
                                const WellState& well_state,
                                std::vector<double>& transport_src)
    {
        int nc = props.numCells();
        transport_src.clear();
        transport_src.resize(nc, 0.0);
        // Well contributions.
        if (wells) {
            const int nw = wells->number_of_wells;
            const int np = wells->number_of_phases;
            if (np != 2) {
                OPM_THROW(std::runtime_error, "computeTransportSource() requires a 2 phase case.");
            }
            std::vector<double> A(np*np);
            for (int w = 0; w < nw; ++w) {
                const double* comp_frac = wells->comp_frac + np*w;
                for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                    const int perf_cell = wells->well_cells[perf];
                    double perf_rate = well_state.perfRates()[perf];
                    if (perf_rate > 0.0) {
                        // perf_rate is a total inflow reservoir rate, we want a surface water rate.
                        if (wells->type[w] != INJECTOR) {
                            std::cout << "**** Warning: crossflow in well "
                                      << w << " perf " << perf - wells->well_connpos[w]
                                      << " ignored. Reservoir rate was "
                                      << perf_rate/Opm::unit::day << " m^3/day." << std::endl;
                            perf_rate = 0.0;
                        } else {
                            ASSERT(std::fabs(comp_frac[0] + comp_frac[1] - 1.0) < 1e-6);
                            perf_rate *= comp_frac[0]; // Water reservoir volume rate.
                            props.matrix(1, &well_state.perfPress()[perf], comp_frac, &perf_cell, &A[0], 0);
                            perf_rate *= A[0];         // Water surface volume rate.
                        }
                    }
                    transport_src[perf_cell] += perf_rate;
                }
            }
        }
    }


} // namespace Opm
