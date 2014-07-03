/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Statoil ASA.

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

#include <opm/autodiff/BlackoilPropsAdInterface.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/RegionMapping.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
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
            /**
             * Count number of cells in all regions.
             *
             * This value is needed to compute the average (arithmetic
             * mean) hydrocarbon pressure in each region.
             *
             * \tparam RMap Region mapping.  Typically an instance of
             * class Opm::RegionMapping<> from module opm-core.
             *
             * \param[in] m Specific region mapping.
             *
             * \return Array containing number of cells in each region
             * defined by the region mapping.
             */
            template <class RMap>
            Eigen::ArrayXd
            countCells(const RMap& m)
            {
                // Note: Floating point type (double) to represent
                // cell counts is intentional.  The count will be
                // subsequently used to compute average (pressure)
                // values only, and that operation is safer if we
                // guarantee a floating point type here.
                Eigen::ArrayXd n(m.numRegions());

                for (typename RMap::RegionId
                         r = 0, nr = m.numRegions(); r < nr; ++r)
                {
                    n(r) = double(m.cells(r).size());
                }

                return n;
            }

            /**
             * Extract representative cell in each region.
             *
             * These are the cells for which fluid properties will be
             * computed.
             *
             * \tparam Cells Type of cell container.  Must be
             * commensurable with the properties object's input
             * requirements and support a single (integer) argument
             * constructor that specifies the number of regions.
             * Typically \code std::vector<int> \endcode.
             *
             * \tparam RMap Region mapping.  Typically an instance of
             * class Opm::RegionMapping<> from module opm-core.
             *
             * \param[in] m Specific region mapping.
             *
             * \return Array of representative cells, one cell in each
             * region defined by @c m.
             */
            template <class Cells, class RMap>
            Cells
            representative(const RMap& m)
            {
                Cells c(m.numRegions());

                for (typename RMap::RegionId
                         r = 0, nr = m.numRegions(); r < nr; ++r)
                {
                    c[r] = *m.cells(r).begin();
                }

                return c;
            }

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
         * \tparam Property Fluid property object.  Expected to
         * feature the formation volume factor functions of the
         * BlackoilPropsAdInterface.
         *
         * \tparam Region Type of a forward region mapping.  Expected
         * to provide indexed access through \code operator[]()
         * \endcode as well as inner types \c value_type, \c
         * size_type, and \c const_iterator.  Typically \code
         * std::vector<int> \endcode.
         */
        template <class Property, class Region>
        class SurfaceToReservoirVoidage {
        public:
            /**
             * Constructor.
             *
             * \param[in] props Fluid property object.
             *
             * \param[in] region Forward region mapping.  Often
             * corresponds to the "FIPNUM" mapping of an ECLIPSE input
             * deck.
             */
            SurfaceToReservoirVoidage(const Property& props,
                                      const Region&   region)
                : props_   (props)
                , rmap_    (region)
                , repcells_(Details::representative<typename Property::Cells>(rmap_))
                , ncells_  (Details::countCells(rmap_))
                , p_avg_   (rmap_.numRegions())
                , Rmax_    (rmap_.numRegions(), props.numPhases())
            {}

            /**
             * Compute average hydrocarbon pressure and maximum
             * dissolution and evaporation at average hydrocarbon
             * pressure in all regions in field.
             *
             * Fluid properties are evaluated at average hydrocarbon
             * pressure for purpose of conversion from surface rate to
             * reservoir voidage rate.
             *
             * \param[in] state Dynamic reservoir state.
             */
            void
            defineState(const BlackoilState& state)
            {
                averagePressure(state);
                calcRmax();
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
             * \param[in] in Single tuple of active component rates at
             * surface conditions.
             *
             * \param[in] r Fluid-in-place region to which the
             * component rates correspond.
             *
             * \param[out] coeff Surface-to-reservoir conversion
             * coefficients for all active phases, corresponding to
             * input rates \c in in region \c r.
             */
            template <class Input,
                      class Coeff>
            void
            calcCoeff(const Input& in, const RegionId r, Coeff& coeff)
            {
                typedef typename Property::V V;

                const PhaseUsage&               pu = props_.phaseUsage();
                const V&                        p  = getRegPress(r);
                const typename Property::Cells& c  = getRegCell (r);

                const int iw = Details::PhasePos::water(pu);
                const int io = Details::PhasePos::oil  (pu);
                const int ig = Details::PhasePos::gas  (pu);

                std::fill(& coeff[0], & coeff[0] + props_.numPhases(), 0.0);

                if (Details::PhaseUsed::water(pu)) {
                    // q[w]_r = q[w]_s / bw

                    const V& bw = props_.bWat(p, c);

                    coeff[iw] = 1.0 / bw(0);
                }

                const Miscibility& m = calcMiscibility(in, r);

                // Determinant of 'R' matrix
                const double detR = 1.0 - (m.rs(0) * m.rv(0));

                if (Details::PhaseUsed::oil(pu)) {
                    // q[o]_r = 1/(bo * (1 - rs*rv)) * (q[o]_s - rv*q[g]_s)

                    const V&     bo  = props_.bOil(p, m.rs, m.cond, c);
                    const double den = bo(0) * detR;

                    coeff[io] += 1.0 / den;

                    if (Details::PhaseUsed::gas(pu)) {
                        coeff[ig] -= m.rv(0) / den;
                    }
                }

                if (Details::PhaseUsed::gas(pu)) {
                    // q[g]_r = 1/(bg * (1 - rs*rv)) * (q[g]_s - rs*q[o]_s)

                    const V&     bg  = props_.bGas(p, m.rv, m.cond, c);
                    const double den = bg(0) * detR;

                    coeff[ig] += 1.0 / den;

                    if (Details::PhaseUsed::oil(pu)) {
                        coeff[io] -= m.rs(0) / den;
                    }
                }
            }

        private:
            /**
             * Fluid property object.
             */
            const Property&  props_;

            /**
             * "Fluid-in-place" region mapping (forward and reverse).
             */
            const RegionMapping<Region> rmap_;

            /**
             * Representative cells in each FIP region.
             */
            const typename Property::Cells repcells_;

            /**
             * Number of cells in each region.
             *
             * Floating-point type (double) for purpose of average
             * pressure calculation.
             */
            const Eigen::ArrayXd ncells_;

            /**
             * Average hydrocarbon pressure in each FIP region.
             */
            Eigen::ArrayXd p_avg_;

            /**
             * Maximum dissolution and evaporation ratios at average
             * hydrocarbon pressure.
             *
             * Size (number of regions)-by-(number of fluid phases).
             * Water value is, strictly speaking, wasted if present.
             */
            Eigen::ArrayXXd Rmax_;

            /**
             * Aggregate structure defining fluid miscibility
             * conditions in single region with particular input
             * surface rates.
             */
            struct Miscibility {
                Miscibility()
                    : rs  (1)
                    , rv  (1)
                    , cond(1)
                {
                    rs << 0.0;
                    rv << 0.0;
                }

                /**
                 * Dissolved gas-oil ratio at particular component oil
                 * and gas rates at surface conditions.
                 *
                 * Limited by "RSmax" at average hydrocarbon pressure
                 * in region.
                 */
                typename Property::V rs;

                /**
                 * Evaporated oil-gas ratio at particular component oil
                 * and gas rates at surface conditions.
                 *
                 * Limited by "RVmax" at average hydrocarbon pressure
                 * in region.
                 */
                typename Property::V rv;

                /**
                 * Fluid condition in representative region cell.
                 *
                 * Needed for purpose of FVF evaluation.
                 */
                std::vector<PhasePresence> cond;
            };

            /**
             * Compute average hydrocarbon pressure in all regions.
             *
             * \param[in] state Dynamic reservoir state.
             */
            void
            averagePressure(const BlackoilState& state)
            {
                p_avg_.setZero();

                const std::vector<double>& p = state.pressure();
                for (std::vector<double>::size_type
                         i = 0, n = p.size(); i < n; ++i)
                {
                    p_avg_(rmap_.region(i)) += p[i];
                }

                p_avg_ /= ncells_;
            }

            /**
             * Compute maximum dissolution and evaporation ratios at
             * average hydrocarbon pressure.
             *
             * Uses the pressure value computed by averagePressure()
             * and must therefore be called *after* that method.
             */
            void
            calcRmax()
            {
                Rmax_.setZero();

                const PhaseUsage& pu = props_.phaseUsage();

                if (Details::PhaseUsed::oil(pu) &&
                    Details::PhaseUsed::gas(pu))
                {
                    const Eigen::ArrayXXd::Index
                        io = Details::PhasePos::oil(pu),
                        ig = Details::PhasePos::gas(pu);

                    // Note: Intentionally does not take capillary
                    // pressure into account.  This facility uses the
                    // average *hydrocarbon* pressure rather than
                    // average phase pressure.
                    Rmax_.col(io) = props_.rsSat(p_avg_, repcells_);
                    Rmax_.col(ig) = props_.rvSat(p_avg_, repcells_);
                }
            }

            /**
             * Compute fluid conditions in particular region for a
             * given set of component rates at surface conditions.
             *
             * \tparam Input Type representing collection of (active)
             * component rates at surface conditions.  Must support
             * direct indexing through \code operator[]()\endcode.
             *
             * \param[in] in Single tuple of active component rates at
             * surface conditions.
             *
             * \param[in] r Fluid-in-place region to which the
             * component rates correspond.
             *
             * \return Fluid conditions in region \c r corresponding
             * to surface component rates \c in.
             */
            template <class Input>
            Miscibility
            calcMiscibility(const Input& in, const RegionId r) const
            {
                const PhaseUsage& pu = props_.phaseUsage();

                const int io = Details::PhasePos::oil(pu);
                const int ig = Details::PhasePos::gas(pu);

                Miscibility m;
                PhasePresence& cond = m.cond[0];

                if (Details::PhaseUsed::water(pu)) {
                    cond.setFreeWater();
                }

                if (Details::PhaseUsed::oil(pu)) {
                    cond.setFreeOil();

                    if (Details::PhaseUsed::gas(pu)) {
                        const double rsmax = Rmax_(r, io);
                        const double rs =
                            (0.0 < std::abs(in[io]))
                            ? in[ig] / in[io]
                            : (0.0 < std::abs(in[ig])) ? rsmax : 0.0;

                        if (rsmax < rs) {
                            cond.setFreeGas();
                        }

                        m.rs(0) = std::min(rs, rsmax);
                    }
                }

                if (Details::PhaseUsed::gas(pu)) {
                    if (! Details::PhaseUsed::oil(pu)) {
                        // Oil *NOT* active -- not really supported.
                        cond.setFreeGas();
                    }

                    if (Details::PhaseUsed::oil(pu)) {
                        const double rvmax = Rmax_(r, ig);
                        const double rv =
                            (0.0 < std::abs(in[ig]))
                            ? (in[io] / in[ig])
                            : (0.0 < std::abs(in[io])) ? rvmax : 0.0;

                        m.rv(0) = std::min(rv, rvmax);
                    }
                }

                return m;
            }

            /**
             * Retrieve average hydrocarbon pressure in region.
             *
             * \param[in] r Particular region.
             *
             * \return Average hydrocarbon pressure in region \c r.
             */
            typename Property::V
            getRegPress(const RegionId r) const
            {
                typename Property::V p(1);
                p << p_avg_(r);

                return p;
            }

            /**
             * Retrieve representative cell of region
             *
             * \param[in] r Particular region.
             *
             * \return Representative cell of region \c r.
             */
            typename Property::Cells
            getRegCell(const RegionId r) const
            {
                typename Property::Cells c(1, repcells_[r]);

                return c;
            }
        };
    } // namespace RateConverter
} // namespace Opm

#endif  /* OPM_RATECONVERTER_HPP_HEADER_INCLUDED */
