/*
  Copyright 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2017, IRIS

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

#ifndef OPM_RATECONVERTER_HPP_HEADER_INCLUDED
#define OPM_RATECONVERTER_HPP_HEADER_INCLUDED

#include <opm/grid/utility/RegionMapping.hpp>

#include <opm/simulators/wells/RegionAttributeHelpers.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <array>
#include <cassert>
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
    namespace RateConverter {

        /**
         * Convert component rates at surface conditions to phase
         * (voidage) rates at reservoir conditions.
         *
         * The conversion uses fluid properties evaluated at average
         * hydrocarbon pressure in regions or field.
         *
         * \tparam FluidSystem Fluid system class. Expected to be a
         *    BlackOilFluidSystem
         *
         * \tparam Region Type of a forward region mapping.  Expected to
         *    provide indexed access through \code operator[]() \endcode as
         *    well as inner types \c value_type, \c size_type, and \c
         *    const_iterator.  Typically \code std::vector<int> \endcode.
         */
        template <class FluidSystem, class Region>
        class SurfaceToReservoirVoidage {
        public:
            using Scalar = typename FluidSystem::Scalar;
            /**
             * Constructor.
             *
             * \param[in] region Forward region mapping.  Often corresponds
             * to the "FIPNUM" mapping of an ECLIPSE input deck.
             */
            SurfaceToReservoirVoidage(const Region&     region)
                : rmap_      (region)
                , attr_      (rmap_, Attributes())
            {}

            /**
             * Compute pore volume averaged hydrocarbon state pressure, rs and rv.
             *
             * Fluid properties are evaluated at average hydrocarbon
             * state for purpose of conversion from surface rate to
             * reservoir voidage rate.
             *
             */
            template <typename ElementContext, class Simulator>
            void defineState(const Simulator& simulator)
            {
                // create map from cell to region and set all attributes to
                // zero
                for (const auto& reg : rmap_.activeRegions()) {
                    auto& ra = attr_.attributes(reg);
                    ra.pressure = 0.0;
                    ra.temperature = 0.0;
                    ra.rs = 0.0;
                    ra.rv = 0.0;
                    ra.pv = 0.0;
                    ra.saltConcentration = 0.0;
                    ra.rsw = 0.0;
                    ra.rvw = 0.0;
                }

                // quantities for pore volume average
                std::unordered_map<RegionId, Attributes> attributes_pv;

                // quantities for hydrocarbon volume average
                std::unordered_map<RegionId, Attributes> attributes_hpv;

                for (const auto& reg : rmap_.activeRegions()) {
                    attributes_pv.insert({reg, Attributes()});
                    attributes_hpv.insert({reg, Attributes()});
                }

                ElementContext elemCtx( simulator );
                const auto& gridView = simulator.gridView();
                const auto& comm = gridView.comm();

                OPM_BEGIN_PARALLEL_TRY_CATCH();
                for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                    const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                    const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                    const auto& fs = intQuants.fluidState();
                    // use pore volume weighted averages.
                    const Scalar pv_cell =
                            simulator.model().dofTotalVolume(simulator.vanguard().gridEquilIdxToGridIdx(cellIdx))
                            * intQuants.porosity().value();

                    // only count oil and gas filled parts of the domain
                    Scalar hydrocarbon = 1.0;
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        hydrocarbon -= fs.saturation(FluidSystem::waterPhaseIdx).value();
                    }

                    const int reg = rmap_.region(cellIdx);
                    assert(reg >= 0);

                    // sum p, rs, rv, and T.
                    const Scalar hydrocarbonPV = pv_cell*hydrocarbon;
                    if (hydrocarbonPV > 0.) {
                        auto& attr = attributes_hpv[reg];
                        attr.pv += hydrocarbonPV;
                        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                            attr.rs += fs.Rs().value() * hydrocarbonPV;
                            attr.rv += fs.Rv().value() * hydrocarbonPV;
                        }
                        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * hydrocarbonPV;
                            attr.temperature += fs.temperature(FluidSystem::oilPhaseIdx).value() * hydrocarbonPV;
                        } else {
                            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
                            attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * hydrocarbonPV;
                            attr.temperature += fs.temperature(FluidSystem::gasPhaseIdx).value() * hydrocarbonPV;
                        }
                        attr.saltConcentration += fs.saltConcentration().value() * hydrocarbonPV;
                        if (FluidSystem::enableDissolvedGasInWater()) {
                            attr.rsw += fs.Rsw().value() * hydrocarbonPV; // scale with total volume?
                        }
                        if (FluidSystem::enableVaporizedWater()) {
                            attr.rvw += fs.Rvw().value() * hydrocarbonPV; // scale with total volume?
                        }
                    }

                    if (pv_cell > 0.) {
                        auto& attr = attributes_pv[reg];
                        attr.pv += pv_cell;
                        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                            attr.rs += fs.Rs().value() * pv_cell;
                            attr.rv += fs.Rv().value() * pv_cell;
                        }
                        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * pv_cell;
                            attr.temperature += fs.temperature(FluidSystem::oilPhaseIdx).value() * pv_cell;
                        } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                             attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * pv_cell;
                             attr.temperature += fs.temperature(FluidSystem::gasPhaseIdx).value() * pv_cell;
                        } else {
                            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
                            attr.pressure += fs.pressure(FluidSystem::waterPhaseIdx).value() * pv_cell;
                            attr.temperature += fs.temperature(FluidSystem::waterPhaseIdx).value() * pv_cell;
                        }
                        attr.saltConcentration += fs.saltConcentration().value() * pv_cell;
                        if (FluidSystem::enableDissolvedGasInWater()) {
                            attr.rsw += fs.Rsw().value() * pv_cell;
                        }
                        if (FluidSystem::enableVaporizedWater()) {
                            attr.rvw += fs.Rvw().value() * pv_cell;
                        }
                    }
                }

                OPM_END_PARALLEL_TRY_CATCH("SurfaceToReservoirVoidage::defineState() failed: ", simulator.vanguard().grid().comm());

                this->sumRates(attributes_hpv,
                               attributes_pv,
                               comm);
            }

            /**
             * Region identifier.
             *
             * Integral type.
             */
            using RegionId = typename RegionMapping<Region>::RegionId;

            /**
             * Compute coefficients for surface-to-reservoir voidage
             * conversion.
             *
             * \tparam Input Type representing contiguous collection
             * of component rates at surface conditions.  Must support
             * direct indexing through \code operator[]()\endcode.
             *
             * \tparam Coeff Type representing contiguous collection
             * of surface-to-reservoir conversion coefficients.  Must
             * support direct indexing through \code operator[]()
             * \endcode.
             *
             *
             * \param[in] r Fluid-in-place region of the well
             * \param[in] pvtRegionIdx PVT region of the well
             *
             *
             * \param[out] coeff Surface-to-reservoir conversion
             * coefficients that can be used to compute total reservoir
             * volumes from surface volumes with the formula
             *               q_{rT} = \sum_p coeff[p] q_{sp}.
             * However, individual phase reservoir volumes cannot be calculated from
             * these coefficients (i.e. q_{rp} is not equal to coeff[p] q_{sp})
             * since they can depend on more than one surface volume rate when
             * we have dissolved gas or vaporized oil.
             */
            template <class Coeff>
            void
            calcCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const;

            template <class Coeff , class Rates>
            void
            calcCoeff(const RegionId r, const int pvtRegionIdx, const Rates& surface_rates, Coeff& coeff) const;

            template <class Coeff>
            void
            calcCoeff(  const int     pvtRegionIdx,
                        const Scalar  p,
                        const Scalar  rs,
                        const Scalar  rv,
                        const Scalar  rsw,
                        const Scalar  rvw,
                        const Scalar  T,
                        const Scalar  saltConcentration,
                        Coeff&        coeff) const;

            template <class Coeff>
            void
            calcInjCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const;

            /**
             * Convert surface volume flow rates to reservoir voidage flow
             * rates.
             *
             * State dependent version.  Client must call \code
             * defineState() \endcode prior to invoking this member
             * function.
             *
             * \tparam Rates Type representing contiguous collection of
             *    surface flow rates.  Must support direct indexing through
             *    \code operator[]() \endcode.
             *
             * \param[in] r Zero based fluid-in-place region index.
             *
             * \param[in] pvtRegionIdx Zero based PVT region index.
             *
             * \param[in] surface_rates surface volume flow rates for all
             *    active phases.
             *
             * \param[out] voidage_rates reservoir volume flow rates for all
             *    active phases.
             */
            template <class Rates>
            void calcReservoirVoidageRates(const RegionId r,
                                           const int      pvtRegionIdx,
                                           const Rates&   surface_rates,
                                           Rates&         voidage_rates) const;

            /**
             * Convert surface volume flow rates to reservoir voidage flow
             * rates.
             *
             * State independent version.
             *
             * \tparam Rates Type representing contiguous collection of
             *    surface flow rates.  Must support direct indexing through
             *    \code operator[]() \endcode.
             *
             * \param[in] pvtRegionIdx PVT region.
             *
             * \param[in] p Fluid pressure.
             *
             * \param[in] rs Dissolved gas/oil ratio.
             *
             * \param[in] rv Vaporised oil/gas ratio.
             *
             * \param[in] rsw Dissolved gas/water ratio.
             *
             * \param[in] rwv Vaporised water/gas ratio.
             *
             * \param[in] T Temperature.  Unused in non-thermal simulation
             *    runs.
             *
             * \param[in] saltConcentration Salt concentration.  Unused in
             *    simulation runs without salt precipitation.
             *
             * \param[in] surface_rates Surface volume flow rates for all
             *    active phases.
             *
             * \param[out] voidage_rates Reservoir volume flow rates for all
             *    active phases.
             */
            template <typename SurfaceRates, typename VoidageRates>
            void calcReservoirVoidageRates(const int           pvtRegionIdx,
                                           const Scalar        p,
                                           const Scalar        rs,
                                           const Scalar        rv,
                                           const Scalar        rsw,
                                           const Scalar        rvw,
                                           const Scalar        T,
                                           const Scalar        saltConcentration,
                                           const SurfaceRates& surface_rates,
                                           VoidageRates&       voidage_rates) const;

            template <class Rates>
            std::pair<Scalar, Scalar>
            inferDissolvedVaporisedRatio(const Scalar rsMax,
                                         const Scalar rvMax,
                                         const Rates& surface_rates) const;

            /**
             * Compute coefficients for surface-to-reservoir voidage
             * conversion for solvent.
             *
             *
             * \param[in] r Fluid-in-place region of the well
             * \param[in] pvtRegionIdx PVT region of the well
             *
             *
             * \param[out] double Surface-to-reservoir conversion
             * coefficients for solvent.
             */
            template <class SolventModule>
            void
            calcCoeffSolvent(const RegionId r, const int pvtRegionIdx, Scalar& coeff) const
            {
                const auto& ra = attr_.attributes(r);
                const Scalar p = ra.pressure;
                const Scalar T = ra.temperature;
                const auto& solventPvt = SolventModule::solventPvt();
                const Scalar bs = solventPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p);
                coeff = 1.0 / bs;
            }

        private:
            /**
             * "Fluid-in-place" region mapping (forward and reverse).
             */
            const RegionMapping<Region> rmap_;

            /**
             * Derived property attributes for each active region.
             */
            struct Attributes {
                Attributes()
                    : data{0.0}
                    , pressure(data[0])
                    , temperature(data[1])
                    , rs(data[2])
                    , rv(data[3])
                    , rsw(data[4])
                    , rvw(data[5])
                    , pv(data[6])
                    , saltConcentration(data[7])
                {}

                Attributes(const Attributes& rhs)
                    : Attributes()
                {
                    this->data = rhs.data;
                }

                Attributes& operator=(const Attributes& rhs)
                {
                    this->data = rhs.data;
                    return *this;
                }

                std::array<Scalar,8> data;
                Scalar& pressure;
                Scalar& temperature;
                Scalar& rs;
                Scalar& rv;
                Scalar& rsw;
                Scalar& rvw;
                Scalar& pv;
                Scalar& saltConcentration;
            };

            void sumRates(std::unordered_map<RegionId,Attributes>& attributes_hpv,
                          std::unordered_map<RegionId,Attributes>& attributes_pv,
                          Parallel::Communication comm);

            RegionAttributeHelpers::RegionAttributes<RegionId, Attributes> attr_;
        };

    } // namespace RateConverter
} // namespace Opm

#endif  /* OPM_RATECONVERTER_HPP_HEADER_INCLUDED */
