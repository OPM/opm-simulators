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

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/grid/utility/RegionMapping.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>

#include <dune/grid/common/gridenums.hh>
#include <algorithm>
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
    namespace RateConverter {
        /**
         * Convenience tools for implementing the rate conversion
         * facility.
         */
        namespace Details {
            namespace Select {
                template <class RegionID, bool>
                struct RegionIDParameter
                {
                    using type =
                        typename std::remove_reference<RegionID>::type &;
                };

                template <class RegionID>
                struct RegionIDParameter<RegionID, true>
                {
                    using type = RegionID;
                };
            } // Select

            /**
             * \brief Computes the temperature, pressure, and counter increment.
             *
             * In a parallel run only cells owned contribute to the cell average.
             * \tparam is_parallel Whether this is a parallel run.
             */
             template<bool is_parallel>
            struct AverageIncrementCalculator
            {
                /**
                 * \brief Computes the temperature, pressure, and counter increment.
                 * \param pressure    The pressure.
                 * \param temperature The temperature.
                 * \param rs          The rs.
                 * \param rv          The rv.
                 * \param cell        The current cell index.
                 * \param ownership   A vector indicating whether a cell is owned
                 *                    by this process (value 1), or not (value 0).
                 * \param cell        The cell index.
                 */
                std::tuple<double, double, double, double, int>
                operator()(const std::vector<double>& pressure,
                           const std::vector<double>& temperature,
                           const std::vector<double>& rs,
                           const std::vector<double>& rv,
                           const std::vector<double>& ownership,
                           std::size_t cell){
                    if ( ownership[cell] )
                    {
                        return std::make_tuple(pressure[cell],
                                               temperature[cell],
                                               rs[cell],
                                               rv[cell],
                                               1);
                    }
                    else
                    {
                        return std::make_tuple(0, 0, 0, 0, 0);
                    }
                }
            };
            template<>
            struct AverageIncrementCalculator<false>
            {
                std::tuple<double, double, double, double, int>
                operator()(const std::vector<double>& pressure,
                           const std::vector<double>& temperature,
                           const std::vector<double>& rs,
                           const std::vector<double>& rv,
                           const std::vector<double>&,
                           std::size_t cell){
                    return std::make_tuple(pressure[cell],
                                           temperature[cell],
                                           rs[cell],
                                           rv[cell],
                                           1);
                }
            };
            /**
             * Provide mapping from Region IDs to user-specified collection
             * of per-region attributes.
             *
             * \tparam RegionId Region identifier type.  Must be hashable by
             * \code std::hash<> \endcode.  Typically a built-in integer
             * type--e.g., \c int.
             *
             * \tparam Attributes User-defined type that represents
             * collection of attributes that have meaning in a per-region
             * aggregate sense.  Must be copy-constructible.
             */
            template <typename RegionId, class Attributes>
            class RegionAttributes
            {
            public:
                /**
                 * Expose \c RegionId as a vocabulary type for use in query
                 * methods.
                 */
                using RegionID =
                    typename Select::RegionIDParameter
                    <RegionId, std::is_integral<RegionId>::value>::type;

                /**
                 * Constructor.
                 *
                 * \tparam RMap Class type that implements the RegionMapping
                 * protocol.  Typically an instantiation of \code
                 * Opm::RegionMapping<> \endcode.
                 *
                 * \param[in] rmap Specific region mapping that provides
                 * reverse lookup from regions to cells.
                 *
                 * \param[in] attr Pre-constructed initialiser for \c
                 * Attributes.
                 */
                template <class RMap>
                RegionAttributes(const RMap&       rmap,
                                 const Attributes& attr)
                {
                    using VT = typename AttributeMap::value_type;

                    for (const auto& r : rmap.activeRegions()) {
                        auto v = std::make_unique<Value>(attr);

                        const auto stat = attr_.insert(VT(r, std::move(v)));

                        if (stat.second) {
                            // New value inserted.
                            const auto& cells = rmap.cells(r);

                            assert (! cells.empty());

                            // Region's representative cell.
                            stat.first->second->cell_ = cells[0];
                        }
                    }
                }

                /**
                 * Retrieve representative cell in region.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Representative cell in region \p reg.
                 */
                int cell(const RegionID reg) const
                {
                    return this->find(reg).cell_;
                }

                /**
                 * Request read-only access to region's attributes.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Read-only access to region \p reg's per-region
                 * attributes.
                 */
                const Attributes& attributes(const RegionID reg) const
                {
                    return this->find(reg).attr_;
                }

                /**
                 * Request modifiable access to region's attributes.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Read-write access to region \p reg's per-region
                 * attributes.
                 */
                Attributes& attributes(const RegionID reg)
                {
                    return this->find(reg).attr_;
                }

            private:
                /**
                 * Aggregate per-region attributes along with region's
                 * representative cell.
                 */
                struct Value {
                    Value(const Attributes& attr)
                        : attr_(attr)
                        , cell_(-1)
                    {}

                    Attributes attr_;
                    int        cell_;
                };

                using ID =
                    typename std::remove_reference<RegionId>::type;

                using AttributeMap =
                    std::unordered_map<ID, std::unique_ptr<Value>>;

                AttributeMap attr_;

                /**
                 * Read-only access to region's properties.
                 */
                const Value& find(const RegionID reg) const
                {
                    const auto& i = attr_.find(reg);

                    if (i == attr_.end()) {
                        throw std::invalid_argument("Unknown region ID");
                    }

                    return *i->second;
                }

                /**
                 * Read-write access to region's properties.
                 */
                Value& find(const RegionID reg)
                {
                    const auto& i = attr_.find(reg);

                    if (i == attr_.end()) {
                        throw std::invalid_argument("Unknown region ID");
                    }

                    return *i->second;
                }
            };

            /**
             * Convenience functions for querying presence/absence of
             * active phases.
             */
            namespace PhaseUsed {
                /**
                 * Active water predicate.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Whether or not water is an active phase.
                 */
                inline bool
                water(const PhaseUsage& pu)
                {
                    return pu.phase_used[ BlackoilPhases::Aqua ] != 0;
                }

                /**
                 * Active oil predicate.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Whether or not oil is an active phase.
                 */
                inline bool
                oil(const PhaseUsage& pu)
                {
                    return pu.phase_used[ BlackoilPhases::Liquid ] != 0;
                }

                /**
                 * Active gas predicate.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Whether or not gas is an active phase.
                 */
                inline bool
                gas(const PhaseUsage& pu)
                {
                    return pu.phase_used[ BlackoilPhases::Vapour ] != 0;
                }
            } // namespace PhaseUsed

            /**
             * Convenience functions for querying numerical IDs
             * ("positions") of active phases.
             */
            namespace PhasePos {
                /**
                 * Numerical ID of active water phase.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Non-negative index/position of water if
                 * active, -1 if not.
                 */
                inline int
                water(const PhaseUsage& pu)
                {
                    int p = -1;

                    if (PhaseUsed::water(pu)) {
                        p = pu.phase_pos[ BlackoilPhases::Aqua ];
                    }

                    return p;
                }

                /**
                 * Numerical ID of active oil phase.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Non-negative index/position of oil if
                 * active, -1 if not.
                 */
                inline int
                oil(const PhaseUsage& pu)
                {
                    int p = -1;

                    if (PhaseUsed::oil(pu)) {
                        p = pu.phase_pos[ BlackoilPhases::Liquid ];
                    }

                    return p;
                }

                /**
                 * Numerical ID of active gas phase.
                 *
                 * \param[in] pu Active phase object.
                 *
                 * \return Non-negative index/position of gas if
                 * active, -1 if not.
                 */
                inline int
                gas(const PhaseUsage& pu)
                {
                    int p = -1;

                    if (PhaseUsed::gas(pu)) {
                        p = pu.phase_pos[ BlackoilPhases::Vapour ];
                    }

                    return p;
                }
            } // namespace PhasePos
        } // namespace Details

        /**
         * Convert component rates at surface conditions to phase
         * (voidage) rates at reservoir conditions.
         *
         * The conversion uses fluid properties evaluated at average
         * hydrocarbon pressure in regions or field.
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
        class SurfaceToReservoirVoidage {
        public:
            /**
             * Constructor.
             *
             * \param[in] region Forward region mapping.  Often
             * corresponds to the "FIPNUM" mapping of an ECLIPSE input
             * deck.
             */
            SurfaceToReservoirVoidage(const PhaseUsage& phaseUsage,
                                      const Region&   region)
                : phaseUsage_(phaseUsage)
                , rmap_ (region)
                , attr_ (rmap_, Attributes())
            {
                active = Details::PhaseUsed::gas(phaseUsage_) || Details::PhaseUsed::oil(phaseUsage_);
            }


            /**
             * Compute pore volume averaged hydrocarbon state pressure, rs and rv.
             *
             * Fluid properties are evaluated at average hydrocarbon
             * state for purpose of conversion from surface rate to
             * reservoir voidage rate.
             *
             */
            template <typename ElementContext, class EbosSimulator>
            void defineState(const EbosSimulator& simulator)
            {
                if (!active)
                    return;

                // create map from cell to region
                // and set all attributes to zero
                for (const auto& reg : rmap_.activeRegions()) {
                    auto& ra = attr_.attributes(reg);
                    ra.pressure = 0.0;
                    ra.temperature = 0.0;
                    ra.rs = 0.0;
                    ra.rv = 0.0;
                    ra.pv = 0.0;
                    ra.saltConcentration = 0.0;

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

                const auto& elemEndIt = gridView.template end</*codim=*/0>();
                for (auto elemIt = gridView.template begin</*codim=*/0>();
                     elemIt != elemEndIt;
                     ++elemIt)
                {

                    const auto& elem = *elemIt;
                    if (elem.partitionType() != Dune::InteriorEntity)
                        continue;

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
                    if (Details::PhaseUsed::water(pu)) {
                        hydrocarbon -= fs.saturation(FluidSystem::waterPhaseIdx).value();
                    }

                    const int reg = rmap_.region(cellIdx);
                    assert(reg >= 0);

                    // sum p, rs, rv, and T.
                    const double hydrocarbonPV = pv_cell*hydrocarbon;
                    if (hydrocarbonPV > 0.) {
                        auto& attr = attributes_hpv[reg];
                        attr.pv += hydrocarbonPV;
                        if (Details::PhaseUsed::oil(pu) && Details::PhaseUsed::gas(pu)) {
                            attr.rs += fs.Rs().value() * hydrocarbonPV;
                            attr.rv += fs.Rv().value() * hydrocarbonPV;
                        }
                        if (Details::PhaseUsed::oil(pu)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * hydrocarbonPV;
                            attr.temperature += fs.temperature(FluidSystem::oilPhaseIdx).value() * hydrocarbonPV;
                        } else {
                            assert(Details::PhaseUsed::gas(pu));
                            attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * hydrocarbonPV;
                            attr.temperature += fs.temperature(FluidSystem::gasPhaseIdx).value() * hydrocarbonPV;
                        }
                        attr.saltConcentration += fs.saltConcentration().value() * hydrocarbonPV;
                    }

                    if (pv_cell > 0.) {
                        auto& attr = attributes_pv[reg];
                        attr.pv += pv_cell;
                        if (Details::PhaseUsed::oil(pu) && Details::PhaseUsed::gas(pu)) {
                            attr.rs += fs.Rs().value() * pv_cell;
                            attr.rv += fs.Rv().value() * pv_cell;
                        }
                        if (Details::PhaseUsed::oil(pu)) {
                            attr.pressure += fs.pressure(FluidSystem::oilPhaseIdx).value() * pv_cell;
                            attr.temperature += fs.temperature(FluidSystem::oilPhaseIdx).value() * pv_cell;
                        } else {
                             assert(Details::PhaseUsed::gas(pu));
                             attr.pressure += fs.pressure(FluidSystem::gasPhaseIdx).value() * pv_cell;
                             attr.temperature += fs.temperature(FluidSystem::gasPhaseIdx).value() * pv_cell;
                        }
                        attr.saltConcentration += fs.saltConcentration().value() * pv_cell;
                    }
                }

                for (const auto& reg : rmap_.activeRegions()) {
                      auto& ra = attr_.attributes(reg);
                      const double hpv_sum = comm.sum(attributes_hpv[reg].pv);
                      // TODO: should we have some epsilon here instead of zero?
                      if (hpv_sum > 0.) {
                          const auto& attri_hpv = attributes_hpv[reg];
                          const double p_hpv_sum = comm.sum(attri_hpv.pressure);
                          const double T_hpv_sum = comm.sum(attri_hpv.temperature);
                          const double rs_hpv_sum = comm.sum(attri_hpv.rs);
                          const double rv_hpv_sum = comm.sum(attri_hpv.rv);
                          const double sc_hpv_sum = comm.sum(attri_hpv.saltConcentration);

                          ra.pressure = p_hpv_sum / hpv_sum;
                          ra.temperature = T_hpv_sum / hpv_sum;
                          ra.rs = rs_hpv_sum / hpv_sum;
                          ra.rv = rv_hpv_sum / hpv_sum;
                          ra.pv = hpv_sum;
                          ra.saltConcentration = sc_hpv_sum / hpv_sum;
                      } else {
                          // using the pore volume to do the averaging
                          const auto& attri_pv = attributes_pv[reg];
                          const double pv_sum = comm.sum(attri_pv.pv);
                          assert(pv_sum > 0.);
                          const double p_pv_sum = comm.sum(attri_pv.pressure);
                          const double T_pv_sum = comm.sum(attri_pv.temperature);
                          const double rs_pv_sum = comm.sum(attri_pv.rs);
                          const double rv_pv_sum = comm.sum(attri_pv.rv);
                          const double sc_pv_sum = comm.sum(attri_pv.saltConcentration);

                          ra.pressure = p_pv_sum / pv_sum;
                          ra.temperature = T_pv_sum / pv_sum;
                          ra.rs = rs_pv_sum / pv_sum;
                          ra.rv = rv_pv_sum / pv_sum;
                          ra.pv = pv_sum;
                          ra.saltConcentration = sc_pv_sum / pv_sum;
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
            calcCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const
            {
                if (!active)
                    return;

                const auto& pu = phaseUsage_;
                const auto& ra = attr_.attributes(r);
                const double p = ra.pressure;
                const double T = ra.temperature;
                const double saltConcentration = ra.saltConcentration;

                const int   iw = Details::PhasePos::water(pu);
                const int   io = Details::PhasePos::oil  (pu);
                const int   ig = Details::PhasePos::gas  (pu);

                std::fill(& coeff[0], & coeff[0] + phaseUsage_.num_phases, 0.0);

                if (Details::PhaseUsed::water(pu)) {
                    // q[w]_r = q[w]_s / bw

                    const double bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, saltConcentration);

                    coeff[iw] = 1.0 / bw;
                }

                // Actual Rs and Rv:
                double Rs = ra.rs;
                double Rv = ra.rv;

                // Determinant of 'R' matrix
                const double detR = 1.0 - (Rs * Rv);

                if (Details::PhaseUsed::oil(pu)) {
                    // q[o]_r = 1/(bo * (1 - rs*rv)) * (q[o]_s - rv*q[g]_s)

                    const double bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rs);
                    const double den = bo * detR;

                    coeff[io] += 1.0 / den;

                    if (Details::PhaseUsed::gas(pu)) {
                        coeff[ig] -= ra.rv / den;
                    }
                }

                if (Details::PhaseUsed::gas(pu)) {
                    // q[g]_r = 1/(bg * (1 - rs*rv)) * (q[g]_s - rs*q[o]_s)

                    const double bg  = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rv);
                    const double den = bg * detR;

                    coeff[ig] += 1.0 / den;

                    if (Details::PhaseUsed::oil(pu)) {
                        coeff[io] -= ra.rs / den;
                    }
                }
            }

            template <class Coeff>
            void
            calcInjCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const
            {
                if (!active)
                    return;

                const auto& pu = phaseUsage_;
                const auto& ra = attr_.attributes(r);
                const double p = ra.pressure;
                const double T = ra.temperature;
                const double saltConcentration = ra.saltConcentration;

                const int   iw = Details::PhasePos::water(pu);
                const int   io = Details::PhasePos::oil  (pu);
                const int   ig = Details::PhasePos::gas  (pu);

                std::fill(& coeff[0], & coeff[0] + phaseUsage_.num_phases, 0.0);

                if (Details::PhaseUsed::water(pu)) {
                    // q[w]_r = q[w]_s / bw

                    const double bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, saltConcentration);

                    coeff[iw] = 1.0 / bw;
                }

                if (Details::PhaseUsed::oil(pu)) {
                    const double bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, 0.0);
                    coeff[io] += 1.0 / bo;
                }

                if (Details::PhaseUsed::gas(pu)) {
                    const double bg = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, 0.0);
                    coeff[ig] += 1.0 / bg;
                }
            }


            /**
             * Converting surface volume rates to reservoir voidage rates
             *
             * \tparam Rates Type representing contiguous collection
             * of surface-to-reservoir conversion coefficients.  Must
             * support direct indexing through \code operator[]()
             * \endcode.
             *
             *
             * \param[in] r Fluid-in-place region of the well
             * \param[in] pvtRegionIdx PVT region of the well
             * \param[in] surface_rates surface voluem rates for
             * all active phases
             *
             * \param[out] voidage_rates reservoir volume rates for
             * all active phases
             */
            template <class Rates >
            void
            calcReservoirVoidageRates(const RegionId r, const int pvtRegionIdx, const Rates& surface_rates,
                                      Rates& voidage_rates) const
            {
                if (!active)
                    return;

                assert(voidage_rates.size() == surface_rates.size());

                std::fill(voidage_rates.begin(), voidage_rates.end(), 0.0);

                const auto& pu = phaseUsage_;
                const auto& ra = attr_.attributes(r);
                const double p = ra.pressure;
                const double T = ra.temperature;
                const double saltConcentration = ra.saltConcentration;

                const int   iw = Details::PhasePos::water(pu);
                const int   io = Details::PhasePos::oil  (pu);
                const int   ig = Details::PhasePos::gas  (pu);

                if (Details::PhaseUsed::water(pu)) {
                    // q[w]_r = q[w]_s / bw

                    const double bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, saltConcentration);

                    voidage_rates[iw] = surface_rates[iw] / bw;
                }

                // Use average Rs and Rv:
                auto a = ra.rs;
                auto b = a;
                if (io >= 0 && ig >= 0) {
                    b = surface_rates[ig]/(surface_rates[io]+1.0e-15);
                }

                double Rs = std::min(a, b);

                a = ra.rv;
                b = a;
                if (io >= 0 && ig >= 0) {
                    b = surface_rates[io]/(surface_rates[ig]+1.0e-15);
                }

                double Rv = std::min(a, b);

                // Determinant of 'R' matrix
                const double detR = 1.0 - (Rs * Rv);

                if (Details::PhaseUsed::oil(pu)) {
                    // q[o]_r = 1/(bo * (1 - rs*rv)) * (q[o]_s - rv*q[g]_s)

                    const double bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rs);
                    const double den = bo * detR;

                    voidage_rates[io] = surface_rates[io];

                    if (Details::PhaseUsed::gas(pu)) {
                        voidage_rates[io] -= Rv * surface_rates[ig];
                    }

                    voidage_rates[io] /= den;
                }

                if (Details::PhaseUsed::gas(pu)) {
                    // q[g]_r = 1/(bg * (1 - rs*rv)) * (q[g]_s - rs*q[o]_s)

                    const double bg  = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rv);
                    const double den = bg * detR;

                    voidage_rates[ig] = surface_rates[ig];

                    if (Details::PhaseUsed::oil(pu)) {
                        voidage_rates[ig] -= Rs * surface_rates[io];
                    }

                    voidage_rates[ig] /= den;
                }
            }


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
            calcCoeffSolvent(const RegionId r, const int pvtRegionIdx, double& coeff) const
            {
                if (!active)
                    return;

                const auto& ra = attr_.attributes(r);
                const double p = ra.pressure;
                const double T = ra.temperature;
                const auto& solventPvt = SolventModule::solventPvt();
                const double bs = solventPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p);
                coeff = 1.0 / bs;
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

            bool active; //!< True if converter is active

            /**
             * Derived property attributes for each active region.
             */
            struct Attributes {
                Attributes()
                    : pressure   (0.0)
                    , temperature(0.0)
                    , rs(0.0)
                    , rv(0.0)
                    , pv(0.0)
                    , saltConcentration(0.0)
                {}

                double pressure;
                double temperature;
                double rs;
                double rv;
                double pv;
                double saltConcentration;
            };

            Details::RegionAttributes<RegionId, Attributes> attr_;

        };
    } // namespace RateConverter
} // namespace Opm

#endif  /* OPM_RATECONVERTER_HPP_HEADER_INCLUDED */
