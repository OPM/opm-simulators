/*
  Copyright 2021, Equinor

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_REGIONAVERAGECALCULATOR_HPP_HEADER_INCLUDED
#define OPM_REGIONAVERAGECALCULATOR_HPP_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/RegionAttributeHelpers.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * \file
 * Facility for converting component rates at surface conditions to
 * phase (voidage) rates at reservoir conditions.
 *
 * This uses the average hydrocarbon pressure to define fluid
 * properties.  The facility is intended to support Reservoir Voidage
 * rates only ('RESV').
 */

namespace Opm {
    namespace RegionAverageCalculator {

        /**
         * Computes hydrocarbon weighed average pressures over regions
         *
         * \tparam FluidSystem Fluid system class. Expected to be a BlackOilFluidSystem
         *
         * \tparam Region Type of a forward region mapping.  Expected
         * to provide indexed access through \code operator[]()
         * \endcode as well as inner types \c value_type, \c
         * size_type, and \c const_iterator.  Typically \code
         * std::vector<int> \endcode.
         */
        template <class FluidSystem, class Region>
        class AverageRegionalPressure {
        public:
            /**
             * Constructor.
             *
             * \param[in] region Forward region mapping.  Often
             * corresponds to the "FIPNUM" mapping of an ECLIPSE input
             * deck.
             */
            AverageRegionalPressure(const PhaseUsage& phaseUsage,
                                    const Region&   region)
                : phaseUsage_(phaseUsage)
                , rmap_ (region)
                , attr_ (rmap_, Attributes())
            {
            }


            /**
             * Compute pore volume averaged hydrocarbon state pressure
             */
            template <typename ElementContext, class Simulator>
            void defineState(const Simulator& simulator)
            {
                int numRegions = 0;
                const auto& gridView = simulator.gridView();
                const auto& comm = gridView.comm();
                for (const auto& reg : rmap_.activeRegions()) {
                    numRegions = std::max(numRegions, reg);
                }
                numRegions = comm.max(numRegions);
                for (int reg = 1; reg <= numRegions ; ++ reg) {
                    if (!attr_.has(reg))
                        attr_.insert(reg, Attributes());
                }
                // create map from cell to region
                // and set all attributes to zero
                for (int reg = 1; reg <= numRegions ; ++ reg) {
                    auto& ra = attr_.attributes(reg);
                    ra.pressure = 0.0;
                    ra.pv = 0.0;

                }

                // quantities for pore volume average
                std::unordered_map<RegionId, Attributes> attributes_pv;

                // quantities for hydrocarbon volume average
                std::unordered_map<RegionId, Attributes> attributes_hpv;

                for (int reg = 1; reg <= numRegions ; ++ reg) {
                    attributes_pv.insert({reg, Attributes()});
                    attributes_hpv.insert({reg, Attributes()});
                }

                ElementContext elemCtx( simulator );

                OPM_BEGIN_PARALLEL_TRY_CATCH();
                for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                    const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                    const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                    const auto& fs = intQuants.fluidState();
                    // use pore volume weighted averages.
                    const double pv_cell =
                            simulator.model().dofTotalVolume(cellIdx)
                            * intQuants.porosity().value();

                    // only count oil and gas filled parts of the domain
                    double hydrocarbon = 1.0;
                    const auto& pu = phaseUsage_;
                    if (RegionAttributeHelpers::PhaseUsed::water(pu)) {
                        hydrocarbon -= fs.saturation(FluidSystem::waterPhaseIdx).value();
                    }

                    const int reg = rmap_.region(cellIdx);
                    assert(reg >= 0);

                    // sum p, rs, rv, and T.
                    const double hydrocarbonPV = pv_cell*hydrocarbon;
                    if (hydrocarbonPV > 0.) {
                        auto& attr = attributes_hpv[reg];
                        attr.pv += hydrocarbonPV;
                        if (RegionAttributeHelpers::PhaseUsed::oil(pu)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * hydrocarbonPV;
                        } else {
                            assert(RegionAttributeHelpers::PhaseUsed::gas(pu));
                            attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * hydrocarbonPV;
                        }
                    }

                    if (pv_cell > 0.) {
                        auto& attr = attributes_pv[reg];
                        attr.pv += pv_cell;
                        if (RegionAttributeHelpers::PhaseUsed::oil(pu)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * pv_cell;
                        } else if (RegionAttributeHelpers::PhaseUsed::gas(pu)) {
                             attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * pv_cell;
                        } else {
                            assert(RegionAttributeHelpers::PhaseUsed::water(pu));
                            attr.pressure += fs.pressure(FluidSystem::waterPhaseIdx).value() * pv_cell;
                        }
                    }
                }
                OPM_END_PARALLEL_TRY_CATCH("AverageRegionalPressure::defineState(): ", simulator.vanguard().grid().comm());

                for (int reg = 1; reg <= numRegions ; ++ reg) {
                      auto& ra = attr_.attributes(reg);
                      const double hpv_sum = comm.sum(attributes_hpv[reg].pv);
                      // TODO: should we have some epsilon here instead of zero?
                      if (hpv_sum > 0.) {
                          const auto& attri_hpv = attributes_hpv[reg];
                          const double p_hpv_sum = comm.sum(attri_hpv.pressure);
                          ra.pressure = p_hpv_sum / hpv_sum;
                      } else {
                          // using the pore volume to do the averaging
                          const auto& attri_pv = attributes_pv[reg];
                          const double pv_sum = comm.sum(attri_pv.pv);
                          // pore volums can be zero if a fipnum region is empty
                          if (pv_sum > 0) {
                            const double p_pv_sum = comm.sum(attri_pv.pressure);
                            ra.pressure = p_pv_sum / pv_sum;
                          }
                      }
                }
            }

            /**
             * Region identifier.
             *
             * Integral type.
             */
            typedef typename RegionMapping<Region>::RegionId RegionId;

            /**
             * Average pressure
             *
             */
            double
            pressure(const RegionId r) const
            {
                if (r == 0 ) // region 0 is the whole field
                {
                    double pressure = 0.0;
                    int num_active_regions = 0;
                    for (const auto& attr :  attr_.attributes()) {
                        const auto& value = *attr.second;
                        const auto& ra = value.attr_;
                        pressure += ra.pressure;
                        num_active_regions ++;
                    }
                    return pressure / num_active_regions;
                }

                const auto& ra = attr_.attributes(r);
                return ra.pressure;
            }


        private:
            /**
             * Fluid property object.
             */
            const PhaseUsage phaseUsage_;

            /**
             * "Fluid-in-place" region mapping (forward and reverse).
             */
            const RegionMapping<Region> rmap_;

            /**
             * Derived property attributes for each active region.
             */
            struct Attributes {
                Attributes()
                    : pressure(0.0)
                    , pv(0.0)

                {}

                double pressure;
                double pv;

            };

            RegionAttributeHelpers::RegionAttributes<RegionId, Attributes> attr_;

        };
    } // namespace RegionAverageCalculator
} // namespace Opm

#endif  /* OPM_REGIONAVERAGECALCULATOR_HPP_HEADER_INCLUDED */
