/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_INITSTATEEQUIL_HEADER_INCLUDED
#define OPM_INITSTATEEQUIL_HEADER_INCLUDED

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/core/utility/Units.hpp>

#include <array>
#include <cassert>
#include <utility>
#include <vector>

/**
 * \file
 * Facilities for an ECLIPSE-style equilibration-based
 * initialisation scheme (keyword 'EQUIL').
 */
struct UnstructuredGrid;

namespace Opm
{
    /**
     * Types and routines that collectively implement a basic
     * ECLIPSE-style equilibration-based initialisation scheme.
     *
     * This namespace is intentionally nested to avoid name clashes
     * with other parts of OPM.
     */
    namespace equil {
        template <class Props>
        class DensityCalculator;

        /**
         * Facility for calculating phase densities based on the
         * BlackoilPropertiesInterface.
         *
         * Implements the crucial <CODE>operator()(p,svol)</CODE>
         * function that is expected by class EquilReg.
         */
        template <>
        class DensityCalculator< BlackoilPropertiesInterface > {
        public:
            /**
             * Constructor.
             *
             * \param[in] props Implementation of the
             * BlackoilPropertiesInterface.
             *
             * \param[in] c Single cell used as a representative cell
             * in a PVT region.
             */
            DensityCalculator(const BlackoilPropertiesInterface& props,
                              const int                          c)
                : props_(props)
                , c_(1, c)
            {
            }

            /**
             * Compute phase densities of all phases at phase point
             * given by (pressure, surface volume) tuple.
             *
             * \param[in] p Fluid pressure.
             *
             * \param[in] z Surface volumes of all phases.
             *
             * \return Phase densities at phase point.
             */
            std::vector<double>
            operator()(const double               p,
                       const std::vector<double>& z) const
            {
                const int np = props_.numPhases();
                std::vector<double> A(np * np, 0);

                assert (z.size() == std::vector<double>::size_type(np));

                double* dAdp = 0;
                props_.matrix(1, &p, &z[0], &c_[0], &A[0], dAdp);

                std::vector<double> rho(np, 0.0);
                props_.density(1, &A[0], &rho[0]);

                return rho;
            }

        private:
            const BlackoilPropertiesInterface& props_;
            const std::vector<int>             c_;
        };

        /**
         * Types and routines relating to phase mixing in
         * equilibration calculations.
         */
        namespace miscibility {
            /**
             * Type that implements "no phase mixing" policy.
             */
            struct NoMixing {
                /**
                 * Function call.
                 *
                 * \param[in] depth Depth at which to calculate RS
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.  In "no mixing
                 * policy", this is identically zero.
                 */
                double
                operator()(const double /* depth */,
                           const double /* press */) const
                {
                    return 0.0;
                }
            };

            /**
             * Type that implements "dissolved gas-oil ratio"
             * tabulated as a function of depth policy.  Data
             * typically taken from keyword 'RSVD'.
             */
            class RsVD {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] depth Depth nodes.
                 * \param[in] rs Dissolved gas-oil ratio at @c depth.
                 */
                RsVD(const std::vector<double>& depth,
                     const std::vector<double>& rs)
                    : depth_(depth)
                    , rs_(rs)
                {
                }

                /**
                 * Function call.
                 *
                 * \param[in] depth Depth at which to calculate RS
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.
                 */
                double
                operator()(const double depth,
                           const double /* press */) const
                {
                    return linearInterpolation(depth_, rs_, depth);
                }

            private:
                std::vector<double> depth_; /**< Depth nodes */
                std::vector<double> rs_;    /**< Dissolved gas-oil ratio */
            };
        } // namespace miscibility

        /**
         * Forward and reverse mappings between cells and
         * regions/partitions (e.g., the ECLIPSE-style 'SATNUM',
         * 'PVTNUM', or 'EQUILNUM' arrays).
         *
         * \tparam Region Type of a forward region mapping.  Expected
         *                to provide indexed access through
         *                operator[]() as well as inner types
         *                'value_type', 'size_type', and
         *                'const_iterator'.
         */
        template < class Region = std::vector<int> >
        class RegionMapping {
        public:
            /**
             * Constructor.
             *
             * \param[in] reg Forward region mapping, restricted to
             *                active cells only.
             */
            explicit
            RegionMapping(const Region& reg)
                : reg_(reg)
            {
                rev_.init(reg_);
            }

            /**
             * Type of forward (cell-to-region) mapping result.
             * Expected to be an integer.
             */
            typedef typename Region::value_type RegionId;

            /**
             * Type of reverse (region-to-cell) mapping (element)
             * result.
             */
            typedef typename Region::size_type CellId;

            /**
             * Type of reverse region-to-cell range bounds and
             * iterators.
             */
            typedef typename std::vector<CellId>::const_iterator CellIter;

            /**
             * Range of cells.  Result from reverse (region-to-cell)
             * mapping.
             */
            class CellRange {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] b Beginning of range.
                 * \param[in] e One past end of range.
                 */
                CellRange(const CellIter b,
                          const CellIter e)
                    : b_(b), e_(e)
                {}

                /**
                 * Read-only iterator on cell ranges.
                 */
                typedef CellIter const_iterator;

                /**
                 * Beginning of cell range.
                 */
                const_iterator begin() const { return b_; }

                /**
                 * One past end of cell range.
                 */
                const_iterator end()   const { return e_; }

            private:
                const_iterator b_;
                const_iterator e_;
            };

            /**
             * Number of declared regions in cell-to-region mapping.
             */
            RegionId
            numRegions() const { return RegionId(rev_.p.size()) - 1; }

            /**
             * Compute region number of given active cell.
             *
             * \param[in] c Active cell
             * \return Region to which @c c belongs.
             */
            RegionId
            region(const CellId c) const { return reg_[c]; }

            /**
             * Extract active cells in particular region.
             *
             * \param[in] r Region number
             * \returns Range of active cells in region @c r.
             */
            CellRange
            cells(const RegionId r) const {
                const RegionId i = r - rev_.low;
                return CellRange(rev_.c.begin() + rev_.p[i + 0],
                                 rev_.c.begin() + rev_.p[i + 1]);
            }

        private:
            /**
             * Copy of forward region mapping (cell-to-region).
             */
            Region reg_;

            /**
             * Reverse mapping (region-to-cell).
             */
            struct {
                typedef typename std::vector<CellId>::size_type Pos;
                std::vector<Pos>    p;   /**< Region start pointers */
                std::vector<CellId> c;   /**< Region cells */
                RegionId            low; /**< Smallest region number */

                /**
                 * Compute reverse mapping.  Standard linear insertion
                 * sort algorithm.
                 */
                void
                init(const Region& reg)
                {
                    typedef typename Region::const_iterator CI;
                    const std::pair<CI,CI>
                        m = std::minmax_element(reg.begin(), reg.end());

                    low  = *m.first;

                    const typename Region::size_type
                        n = *m.second - low + 1;

                    p.resize(n + 1);  std::fill(p.begin(), p.end(), Pos(0));
                    for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                        p[ reg[i] - low + 1 ] += 1;
                    }

                    for (typename std::vector<Pos>::size_type
                             i = 1, sz = p.size(); i < sz; ++i) {
                        p[0] += p[i];
                        p[i]  = p[0] - p[i];
                    }

                    assert (p[0] ==
                            static_cast<typename Region::size_type>(reg.size()));

                    c.resize(reg.size());
                    for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                        c[ p[ reg[i] - low + 1 ] ++ ] = i;
                    }

                    p[0] = 0;
                }
            } rev_; /**< Reverse mapping instance */
        };

        /**
         * Equilibration record.
         *
         * Layout and contents inspired by first six items of
         * ECLIPSE's 'EQUIL' records.  This is the minimum amount of
         * input data needed to define phase pressures in an
         * equilibration region.
         *
         * Data consists of three pairs of depth and pressure values:
         *   1. main
         *     - @c depth Main datum depth.
         *     - @c press Pressure at datum depth.
         *
         *   2. woc
         *     - @c depth Depth of water-oil contact
         *     - @c press water-oil capillary pressure at water-oil contact.
         *       Capillary pressure defined as "P_oil - P_water".
         *
         *   3. goc
         *     - @c depth Depth of gas-oil contact
         *     - @c press Gas-oil capillary pressure at gas-oil contact.
         *       Capillary pressure defined as "P_gas - P_oil".
         */
        struct EquilRecord {
            struct {
                double depth;
                double press;
            } main, woc, goc;
        };

        /**
         * Aggregate information base of an equilibration region.
         *
         * Provides inquiry methods for retrieving depths of contacs
         * and pressure values as well as a means of calculating fluid
         * densities, dissolved gas-oil ratio and vapourised oil-gas
         * ratios.
         *
         * \tparam DensCalc Type that provides access to a phase
         * density calculation facility.  Must implement an operator()
         * declared as
         * <CODE>
         * std::vector<double>
         * operator()(const double               press,
         *            const std::vector<double>& svol )
         * </CODE>
         * that calculates the phase densities of all phases in @c
         * svol at fluid pressure @c press.
         *
         * \tparam RS Type that provides access to a calculator for
         * (initial) dissolved gas-oil ratios as a function of depth
         * and (oil) pressure.  Must implement an operator() declared
         * as
         * <CODE>
         * double
         * operator()(const double depth,
         *            const double press)
         * </CODE>
         * that calculates the dissolved gas-oil ratio at depth @c
         * depth and (oil) pressure @c press.
         *
         * \tparam RV Type that provides access to a calculator for
         * (initial) vapourised oil-gas ratios as a function of depth
         * and (gas) pressure.  Must implement an operator() declared
         * as
         * <CODE>
         * double
         * operator()(const double depth,
         *            const double press)
         * </CODE>
         * that calculates the vapourised oil-gas ratio at depth @c
         * depth and (gas) pressure @c press.
         */
        template <class DensCalc,
                  class RS = miscibility::NoMixing,
                  class RV = miscibility::NoMixing>
        class EquilReg {
        public:
            /**
             * Constructor.
             *
             * \param[in] rec     Equilibration data of current region.
             * \param[in] density Density calculator of current region.
             * \param[in] rs      Calculator of dissolved gas-oil ratio.
             * \param[in] rv      Calculator of vapourised oil-gas ratio.
             * \param[in] pu      Summary of current active phases.
             */
            EquilReg(const EquilRecord& rec,
                     const DensCalc&    density,
                     const RS&          rs,
                     const RV&          rv,
                     const PhaseUsage&  pu)
                : rec_    (rec)
                , density_(density)
                , rs_     (rs)
                , rv_     (rv)
                , pu_     (pu)
            {
            }

            /**
             * Type of density calculator.
             */
            typedef DensCalc CalcDensity;

            /**
             * Type of dissolved gas-oil ratio calculator.
             */
            typedef RS CalcDissolution;

            /**
             * Type of vapourised oil-gas ratio calculator.
             */
            typedef RV CalcEvaporation;

            /**
             * Datum depth in current region
             */
            double datum()    const { return this->rec_.main.depth; }

            /**
             * Pressure at datum depth in current region.
             */
            double pressure() const { return this->rec_.main.press; }

            /**
             * Depth of water-oil contact.
             */
            double zwoc()     const { return this->rec_.woc .depth; }

            /**
             * water-oil capillary pressure at water-oil contact.
             *
             * \return P_o - P_w at WOC.
             */
            double pcow_woc() const { return this->rec_.woc .press; }

            /**
             * Depth of gas-oil contact.
             */
            double zgoc()     const { return this->rec_.goc .depth; }

            /**
             * Gas-oil capillary pressure at gas-oil contact.
             *
             * \return P_g - P_o at GOC.
             */
            double pcgo_goc() const { return this->rec_.goc .press; }

            /**
             * Retrieve phase density calculator of current region.
             */
            const CalcDensity&
            densityCalculator() const { return this->density_; }

            /**
             * Retrieve dissolved gas-oil ratio calculator of current
             * region.
             */
            const CalcDissolution&
            dissolutionCalculator() const { return this->rs_; }

            /**
             * Retrieve vapourised oil-gas ratio calculator of current
             * region.
             */
            const CalcEvaporation&
            evaporationCalculator() const { return this->rv_; }

            /**
             * Retrieve active fluid phase summary.
             */
            const PhaseUsage&
            phaseUsage() const { return this->pu_; }

        private:
            EquilRecord rec_;     /**< Equilibration data */
            DensCalc    density_; /**< Density calculator */
            RS          rs_;      /**< RS calculator */
            RV          rv_;      /**< RV calculator */
            PhaseUsage  pu_;      /**< Active phase summary */
        };

        /**
         * Compute initial phase pressures by means of equilibration.
         *
         * This function uses the information contained in an
         * equilibration record (i.e., depths and pressurs) as well as
         * a density calculator and related data to vertically
         * integrate the phase pressure ODE
         * \f[
         * \frac{\mathrm{d}p_{\alpha}}{\mathrm{d}z} =
         * \rho_{\alpha}(z,p_{\alpha})\cdot g
         * \f]
         * in which \f$\rho_{\alpha}$ denotes the fluid density of
         * fluid phase \f$\alpha\f$, \f$p_{\alpha}\f$ is the
         * corresponding phase pressure, \f$z\f$ is the depth and
         * \f$g\f$ is the acceleration due to gravity (assumed
         * directed downwords, in the positive \f$z\f$ direction).
         *
         * \tparam Region Type of an equilibration region information
         *                base.  Typically an instance of the EquilReg
         *                class template.
         *
         * \tparam CellRange Type of cell range that demarcates the
         *                cells pertaining to the current
         *                equilibration region.  Must implement
         *                methods begin() and end() to bound the range
         *                as well as provide an inner type,
         *                const_iterator, to traverse the range.
         *
         * \param[in] G     Grid.
         * \param[in] reg   Current equilibration region.
         * \param[in] cells Range that spans the cells of the current
         *                  equilibration region.
         * \param[in] grav  Acceleration of gravity.
         *
         * \return Phase pressures, one vector for each active phase,
         * of pressure values in each cell in the current
         * equilibration region.
         */
        template <class Region, class CellRange>
        std::vector< std::vector<double> >
        phasePressures(const UnstructuredGrid& G,
                       const Region&           reg,
                       const CellRange&        cells,
                       const double            grav = unit::gravity);

        namespace DeckDependent {
            inline
            std::vector<EquilRecord>
            getEquil(const EclipseGridParser& deck)
            {
                if (deck.hasField("EQUIL")) {
                    const EQUIL& eql = deck.getEQUIL();

                    typedef std::vector<EquilLine>::size_type sz_t;
                    const sz_t nrec = eql.equil.size();

                    std::vector<EquilRecord> ret;
                    ret.reserve(nrec);
                    for (sz_t r = 0; r < nrec; ++r) {
                        const EquilLine& rec = eql.equil[r];

                        EquilRecord record =
                            {
                                { rec.datum_depth_             ,
                                  rec.datum_depth_pressure_    }
                                ,
                                { rec.water_oil_contact_depth_ ,
                                  rec.oil_water_cap_pressure_  }
                                ,
                                { rec.gas_oil_contact_depth_   ,
                                  rec.gas_oil_cap_pressure_    }
                            };

                        ret.push_back(record);
                    }

                    return ret;
                }
                else {
                    OPM_THROW(std::domain_error,
                              "Deck does not provide equilibration data.");
                }
            }

            inline
            std::vector<int>
            equilnum(const EclipseGridParser& deck,
                     const UnstructuredGrid&  G   )
            {
                std::vector<int> eqlnum;
                if (deck.hasField("EQLNUM")) {
                    eqlnum = deck.getIntegerValue("EQLNUM");
                }
                else {
                    // No explicit equilibration region.
                    // All cells in region zero.
                    eqlnum.assign(G.number_of_cells, 0);
                }

                return eqlnum;
            }

            template <class InputDeck>
            class PhasePressureComputer;

            template <>
            class PhasePressureComputer<Opm::EclipseGridParser> {
            public:
                PhasePressureComputer(const BlackoilPropertiesInterface& props,
                                      const EclipseGridParser&           deck ,
                                      const UnstructuredGrid&            G    )
                    : pp_(props.numPhases(),
                          std::vector<double>(G.number_of_cells))
                {
                    const std::vector<EquilRecord> rec = getEquil(deck);
                    const RegionMapping<> eqlmap(equilnum(deck, G));

                    calcII(eqlmap, rec, props, G);
                }

                typedef std::vector<double> PVal;
                typedef std::vector<PVal>   PPress;

                const PPress& press() const { return pp_; }

            private:
                typedef DensityCalculator<BlackoilPropertiesInterface> RhoCalc;
                typedef EquilReg<RhoCalc> EqReg;

                PPress pp_;

                template <class RMap>
                void
                calcII(const RMap&                             reg  ,
                       const std::vector< EquilRecord >&       rec  ,
                       const Opm::BlackoilPropertiesInterface& props,
                       const UnstructuredGrid&                 G    )
                {
                    typedef miscibility::NoMixing NoMix;

                    for (typename RMap::RegionId
                             r = 0, nr = reg.numRegions();
                         r < nr; ++r)
                    {
                        const typename RMap::CellRange cells = reg.cells(r);

                        const int repcell = *cells.begin();
                        const RhoCalc calc(props, repcell);

                        const EqReg eqreg(rec[r], calc, NoMix(), NoMix(),
                                          props.phaseUsage());

                        const PPress& res = phasePressures(G, eqreg, cells);

                        for (int p = 0, np = props.numPhases(); p < np; ++p) {
                            PVal&                d = pp_[p];
                            PVal::const_iterator s = res[p].begin();
                            for (typename RMap::CellRange::const_iterator
                                     c = cells.begin(),
                                     e = cells.end();
                                 c != e; ++c, ++s)
                            {
                                d[*c] = *s;
                            }
                        }
                    }
                }
            };
        } // namespace DeckDependent
    } // namespace equil
} // namespace Opm

#include <opm/core/simulator/initStateEquil_impl.hpp>

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
