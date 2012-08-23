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

#ifndef OPM_MISCUTILITIESBLACKOIL_HEADER_INCLUDED
#define OPM_MISCUTILITIESBLACKOIL_HEADER_INCLUDED

#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class BlackoilPropertiesInterface;

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
                                 const std::vector<double>& p,
                                 const std::vector<double>& z,
				 const std::vector<double>& s,
				 const std::vector<double>& src,
				 const double dt,
				 double* injected,
				 double* produced);

    /// @brief Computes total mobility for a set of saturation values.
    /// @param[in]  props     rock and fluid properties
    /// @param[in]  cells     cells with which the saturation values are associated
    /// @param[in]  p         pressure (one value per cell)
    /// @param[in]  z         surface-volume values (for all P phases)
    /// @param[in]  s         saturation values (for all phases)
    /// @param[out] totmob    total mobilities.
    void computeTotalMobility(const Opm::BlackoilPropertiesInterface& props,
			      const std::vector<int>& cells,
                              const std::vector<double>& p,
                              const std::vector<double>& z,
			      const std::vector<double>& s,
			      std::vector<double>& totmob);

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
				   std::vector<double>& omega);


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
                                std::vector<double>&                    pmobc);
    

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
                               std::vector<double>& fractional_flows);


    /// Computes the surface volume densities from saturations by the formula
    ///     z = A s
    /// for a number of data points, where z is the surface volume density,
    /// s is the saturation (both as column vectors) and A is the
    /// phase-to-component relation matrix.
    /// @param[in]  n            number of data points
    /// @param[in]  np           number of phases, must be 2 or 3
    /// @param[in]  A            array containing n square matrices of size num_phases,
    ///                          in Fortran ordering, typically the output of a call
    ///                          to the matrix() method of a BlackoilProperties* class.
    /// @param[in]  saturation   concatenated saturation values (for all P phases)
    /// @param[out] surfacevol   concatenated surface-volume values (for all P phases)
    void computeSurfacevol(const int n,
                           const int np,
                           const double* A,
                           const double* saturation,
                           double* surfacevol);

} // namespace Opm

#endif // OPM_MISCUTILITIESBLACKOIL_HEADER_INCLUDED
