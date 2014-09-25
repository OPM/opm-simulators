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

#ifndef OPM_UTILITIES_HEADER_INCLUDED
#define OPM_UTILITIES_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/wells.h>
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
                                std::vector<double>& transport_src);

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
                                 double& polyprod);
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
                                 double& polyprod);

} //namespace Opm

#endif //OPM_UTILITIES_HEADER_INCLUDED
