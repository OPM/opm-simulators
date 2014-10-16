/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
//#include <opm/autodiff/IncompPropsAdInterface.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/PolymerState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/polymer/fullyimplicit/utilities.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <iostream>
#include <iterator>


namespace Opm
{

    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;
    typedef Eigen::Array<double,
                         Eigen::Dynamic,
                         Eigen::Dynamic,
                         Eigen::RowMajor> DataBlock;
    /// Compute two-phase transport source terms from well terms.
    /// Note: Unlike the incompressible version of this function,
    ///       this version computes surface volume injection rates,
    ///       production rates are still total reservoir volumes.
    /// \param[in]  props         Fluid and rock properties.
    /// \param[in]  wells         Wells data structure.
    /// \param[in]  well_state    Well pressures and fluxes.
    /// \param[out] transport_src The transport source terms. They are to be interpreted depending on sign:
    ///                           (+) positive  inflow of first (water) phase (reservoir volume),
    ///                           (-) negative  total outflow of both phases (reservoir volume).
    void computeTransportSource(const BlackoilPropsAdInterface& props,
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
                            assert(std::fabs(comp_frac[0] + comp_frac[1] - 1.0) < 1e-6);
                            perf_rate *= comp_frac[0]; // Water reservoir volume rate.
                        }
                    }
                    transport_src[perf_cell] += perf_rate;
                }
            }
        }
    }

    /// @brief Computes injected and produced volumes of all phases,
    ///        and injected and produced polymer mass - in the compressible case.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// @param[in]  props     fluid and rock properties.
    /// @param[in]  polyprops polymer properties
    /// @param[in]  state     state variables (pressure, fluxes etc.)
    /// @param[in]  transport_src  if < 0: total reservoir volume outflow,
    ///                       if > 0: first phase *surface volume* inflow.
    /// @param[in]  inj_c     injected concentration by cell
    /// @param[in]  dt        timestep used
    /// @param[out] injected  must point to a valid array with P elements,
    ///                       where P = s.size()/transport_src.size().
    /// @param[out] produced  must also point to a valid array with P elements.
    /// @param[out] polyinj   injected mass of polymer
    /// @param[out] polyprod  produced mass of polymer
    // This function need a incompProps based on Ad.
    /*
    void computeInjectedProduced(const IncompPropsAdInterface& props,
                                 const Opm::PolymerPropsAd& polymer_props,
                                 const PolymerState& state,
                                 const std::vector<double>& transport_src,
                                 const std::vector<double>& inj_c,
                                 const double dt,
                                 double* injected,
                                 double* produced,
                                 double& polyinj,
                                 double& polyprod)
    {
        const int num_cells = transport_src.size();
        if (props.numCells() != num_cells) {
            OPM_THROW(std::runtime_error, "Size of transport_src vector does not match number of cells in props.");
        }
        const int np = props.numPhases();
        if (int(state.saturation().size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of state vectors do not match number of cells.");
        }
		std::vector<int> cells(num_cells);
		const V p = Eigen::Map<const V>(&state.pressure()[0], num_cells, 1);
        const DataBlock s = Eigen::Map<const DataBlock>(&state.saturation()[0], num_cells, np);
		const V sw = s.col(0);
		const V so = s.col(1);
		const V c = Eigen::Map<const V>(&state.concentration()[0], num_cells, 1);
		const V cmax = Eigen::Map<const V>(&state.maxconcentration()[0], num_cells, 1);
		const V trans_src = Eigen::Map<const V>(&transport_src[0], num_cells, 1);
		V src = V::Constant(num_cells, -1.0); // negative is injec, positive is producer.
		for (int cell = 0; cell < num_cells; ++cell) {
			cells[cell] = cell;
			if(transport_src[cell] > 0.0) {
				src[cell] = 1.0;
			}
		}
		const Selector<double> src_selector(src);
		const V one = V::Constant(num_cells, 1.0);
		const V zero = V::Zero(num_cells);
		const std::vector<V> kr = props.relperm(sw, so, cells);

        const V krw_eff = polymer_props.effectiveRelPerm(c, cmax, kr[0]);
        const double* mus = props.viscosity();
		const V inv_muw_eff = polymer_props.effectiveInvWaterVisc(c, mus);
		std::vector<V> mob(np);
		mob[0] = krw_eff * inv_muw_eff;
		mob[1] = kr[1] / mus[1];
		
		const V watmob_c = src_selector.select(mob[0], one);
		const V oilmob_c = src_selector.select(mob[1], zero);
		const V flux = trans_src * dt;
	    const V totmob_c = watmob_c + oilmob_c;
		const V wat_src = flux * (watmob_c / totmob_c);
		const V oil_src = flux * (oilmob_c / totmob_c);
		const V mc = polymer_props.polymerWaterVelocityRatio(c);
		
        polyinj = 0.0;
        polyprod = 0.0;
		std::fill(injected, injected + np , 0.0);
		std::fill(produced, produced + np , 0.0);
		for (int cell = 0; cell < num_cells; ++cell) {
			if (wat_src[cell] < 0) {
				injected[0] += wat_src[cell];
				polyinj += injected[0] * inj_c[cell];
			} else {
				produced[0] += wat_src[cell];
				produced[1] += oil_src[cell];
				polyprod += produced[0] * mc[cell];
			}
		}
    }

    */
    /// @brief Computes injected and produced volumes of all phases,
    ///        and injected and produced polymer mass - in the compressible case.
    /// Note 1: assumes that only the first phase is injected.
    /// Note 2: assumes that transport has been done with an
    ///         implicit method, i.e. that the current state
    ///         gives the mobilities used for the preceding timestep.
    /// @param[in]  props     fluid and rock properties.
    /// @param[in]  polyprops polymer properties
    /// @param[in]  state     state variables (pressure, fluxes etc.)
    /// @param[in]  transport_src  if < 0: total reservoir volume outflow,
    ///                       if > 0: first phase *surface volume* inflow.
    /// @param[in]  inj_c     injected concentration by cell
    /// @param[in]  dt        timestep used
    /// @param[out] injected  must point to a valid array with P elements,
    ///                       where P = s.size()/transport_src.size().
    /// @param[out] produced  must also point to a valid array with P elements.
    /// @param[out] polyinj   injected mass of polymer
    /// @param[out] polyprod  produced mass of polymer
    void computeInjectedProduced(const BlackoilPropsAdInterface& props,
                                 const Opm::PolymerPropsAd& polymer_props,
                                 const PolymerBlackoilState& state,
                                 const std::vector<double>& transport_src,
                                 const std::vector<double>& inj_c,
                                 const double dt,
                                 double* injected,
                                 double* produced,
                                 double& polyinj,
                                 double& polyprod)
    {
        const int num_cells = transport_src.size();
        if (props.numCells() != num_cells) {
            OPM_THROW(std::runtime_error, "Size of transport_src vector does not match number of cells in props.");
        }
        const int np = props.numPhases();
        if (int(state.saturation().size()) != num_cells*np) {
            OPM_THROW(std::runtime_error, "Sizes of state vectors do not match number of cells.");
        }
		std::vector<int> cells(num_cells);
		const V p = Eigen::Map<const V>(&state.pressure()[0], num_cells, 1);
        const DataBlock s = Eigen::Map<const DataBlock>(&state.saturation()[0], num_cells, np);
		const V sw = s.col(0);
		const V so = s.col(1);
		const V c = Eigen::Map<const V>(&state.concentration()[0], num_cells, 1);
		const V cmax = Eigen::Map<const V>(&state.maxconcentration()[0], num_cells, 1);
		const V trans_src = Eigen::Map<const V>(&transport_src[0], num_cells, 1);
		V src = V::Constant(num_cells, -1.0); // negative is injec, positive is producer.
		for (int cell = 0; cell < num_cells; ++cell) {
			cells[cell] = cell;
			if(transport_src[cell] > 0.0) {
				src[cell] = 1.0;
			}
		}
        //Add PhasePresence make muOil() happy.
        std::vector<PhasePresence> phaseCondition(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            phaseCondition[c] = PhasePresence();
            phaseCondition[c].setFreeWater();
            phaseCondition[c].setFreeOil();
        }
		const Selector<double> src_selector(src);
		const V one = V::Constant(num_cells, 1.0);
		const V zero = V::Zero(num_cells);
		const std::vector<V> kr = props.relperm(sw, so, zero, cells);
		const V muw = props.muWat(p, cells);
		const V muo = props.muOil(p, zero, phaseCondition, cells);
        const V krw_eff = polymer_props.effectiveRelPerm(c, cmax, kr[0]);
		const V inv_muw_eff = polymer_props.effectiveInvWaterVisc(c, muw.data());
		std::vector<V> mob(np);
		mob[0] = krw_eff * inv_muw_eff;
		mob[1] = kr[1] / muo;
		
		const V watmob_c = src_selector.select(mob[0], one);
		const V oilmob_c = src_selector.select(mob[1], zero);
		const V flux = trans_src * dt;
	    const V totmob_c = watmob_c + oilmob_c;
		const V wat_src = flux * (watmob_c / totmob_c);
		const V oil_src = flux * (oilmob_c / totmob_c);
		const V mc = polymer_props.polymerWaterVelocityRatio(c);
		
        polyinj = 0.0;
        polyprod = 0.0;
		std::fill(injected, injected + np , 0.0);
		std::fill(produced, produced + np , 0.0);
		for (int cell = 0; cell < num_cells; ++cell) {
			if (wat_src[cell] < 0) {
				injected[0] += wat_src[cell];
				polyinj += injected[0] * inj_c[cell];
			} else {
				produced[0] += wat_src[cell];
				produced[1] += oil_src[cell];
				polyprod += produced[0] * mc[cell];
			}
		}
    }



} //namespace Opm
