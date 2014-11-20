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

#ifndef OPM_EQUILIBRATIONHELPERS_HEADER_INCLUDED
#define OPM_EQUILIBRATIONHELPERS_HEADER_INCLUDED

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/core/utility/RegionMapping.hpp>
#include <opm/core/utility/RootFinders.hpp>

#include <memory>


/*
---- synopsis of EquilibrationHelpers.hpp ----

namespace Opm
{
    namespace Equil {

        template <class Props>
        class DensityCalculator;

        template <>
        class DensityCalculator< BlackoilPropertiesInterface >;

        namespace Miscibility {
            class RsFunction;
            class NoMixing;
            class RsVD;
            class RsSatAtContact;
        }

        struct EquilRecord;

        template <class DensCalc>
        class EquilReg;


        struct PcEq;

        inline double satFromPc(const BlackoilPropertiesInterface& props,
                                const int phase,
                                const int cell,
                                const double target_pc,
                                const bool increasing = false);
        struct PcEqSum
        inline double satFromSumOfPcs(const BlackoilPropertiesInterface& props,
                                      const int phase1,
                                      const int phase2,
                                      const int cell,
                                      const double target_pc);
    } // namespace Equil
} // namespace Opm

---- end of synopsis of EquilibrationHelpers.hpp ----
*/


namespace Opm
{
    /**
     * Types and routines that collectively implement a basic
     * ECLIPSE-style equilibration-based initialisation scheme.
     *
     * This namespace is intentionally nested to avoid name clashes
     * with other parts of OPM.
     */
    namespace Equil {


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
             * \param[in] T Temperature.
             *
             * \param[in] z Surface volumes of all phases.
             *
             * \return Phase densities at phase point.
             */
            std::vector<double>
            operator()(const double               p,
                       const double               T,
                       const std::vector<double>& z) const
            {
                const int np = props_.numPhases();
                std::vector<double> A(np * np, 0);

                assert (z.size() == std::vector<double>::size_type(np));

                double* dAdp = 0;
                props_.matrix(1, &p, &T, &z[0], &c_[0], &A[0], dAdp);

                std::vector<double> rho(np, 0.0);
                props_.density(1, &A[0], &c_[0], &rho[0]);

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
        namespace Miscibility {

            /**
             * Base class for phase mixing functions.
             */
            class RsFunction
            {
            public:
                /**
                 * Function call operator.
                 *
                 * \param[in] depth Depth at which to calculate RS
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RS
                 * value.
                 *
                 * \param[in] temp Temperature at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.
                 */
                virtual double operator()(const double depth,
                                          const double press,
                                          const double temp,
                                          const double sat = 0.0) const = 0;
            };


            /**
             * Type that implements "no phase mixing" policy.
             */
            class NoMixing : public RsFunction {
            public:
                /**
                 * Function call.
                 *
                 * \param[in] depth Depth at which to calculate RS
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RS
                 * value.
                 *
                 * \param[in] temp Temperature at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.  In "no mixing
                 * policy", this is identically zero.
                 */
                double
                operator()(const double /* depth */,
                           const double /* press */,
                           const double /* temp */,
                           const double /* sat */ = 0.0) const
                {
                    return 0.0;
                }
            };


            /**
             * Type that implements "dissolved gas-oil ratio"
             * tabulated as a function of depth policy.  Data
             * typically taken from keyword 'RSVD'.
             */
            class RsVD : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] props      property object
                 * \param[in] cell       any cell in the pvt region
                 * \param[in] depth Depth nodes.
                 * \param[in] rs Dissolved gas-oil ratio at @c depth.
                 */
                RsVD(const BlackoilPropertiesInterface& props,
                     const int cell,
                     const std::vector<double>& depth,
                     const std::vector<double>& rs)
                    : props_(props) 
                    , cell_(cell)
                    , depth_(depth)
                    , rs_(rs)
                {
                    auto pu = props_.phaseUsage();
                    std::fill(z_, z_ + BlackoilPhases::MaxNumPhases, 0.0);
                    z_[pu.phase_pos[BlackoilPhases::Vapour]] = 1e100;
                    z_[pu.phase_pos[BlackoilPhases::Liquid]] = 1.0;
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
                 * \param[in] temp Temperature at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.
                 */
                double
                operator()(const double depth,
                           const double press,
                           const double temp,
                           const double sat_gas = 0.0) const
                {
                    if (sat_gas > 0.0) {
                        return satRs(press, temp);
                    } else {
                        return std::min(satRs(press, temp), linearInterpolation(depth_, rs_, depth));
                    }
                }

            private:
                const BlackoilPropertiesInterface& props_;
                const int cell_;
                std::vector<double> depth_; /**< Depth nodes */
                std::vector<double> rs_;    /**< Dissolved gas-oil ratio */
                double z_[BlackoilPhases::MaxNumPhases];
                mutable double A_[BlackoilPhases::MaxNumPhases * BlackoilPhases::MaxNumPhases];

                double satRs(const double press, const double temp) const
                {
                    props_.matrix(1, &press, &temp, z_, &cell_, A_, 0);
                    // Rs/Bo is in the gas row and oil column of A_.
                    // 1/Bo is in the oil row and column.
                    // Recall also that it is stored in column-major order.
                    const int opos = props_.phaseUsage().phase_pos[BlackoilPhases::Liquid];
                    const int gpos = props_.phaseUsage().phase_pos[BlackoilPhases::Vapour];
                    const int np = props_.numPhases();
                    return A_[np*opos + gpos] / A_[np*opos + opos];
                }
            };


            /**
             * Type that implements "vaporized oil-gas ratio"
             * tabulated as a function of depth policy.  Data
             * typically taken from keyword 'RVVD'.
             */
            class RvVD : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] props      property object
                 * \param[in] cell       any cell in the pvt region
                 * \param[in] depth Depth nodes.
                 * \param[in] rv Dissolved gas-oil ratio at @c depth.
                 */
                RvVD(const BlackoilPropertiesInterface& props,
                     const int cell,
                     const std::vector<double>& depth,
                     const std::vector<double>& rv)
                    : props_(props) 
                    , cell_(cell)
                    , depth_(depth)
                    , rv_(rv)
                {
                    auto pu = props_.phaseUsage();
                    std::fill(z_, z_ + BlackoilPhases::MaxNumPhases, 0.0);
                    z_[pu.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    z_[pu.phase_pos[BlackoilPhases::Liquid]] = 1e100;
                }

                /**
                 * Function call.
                 *
                 * \param[in] depth Depth at which to calculate RV
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RV
                 * value.
                 *
                 * \param[in] temp Temperature at which to calculate RV
                 * value.
                 *
                 * \return Vaporized oil-gas ratio (RV) at depth @c
                 * depth and pressure @c press.
                 */
                double
                operator()(const double depth,
                           const double press,
                           const double temp,
                           const double sat_oil = 0.0 ) const
                {
                    if (sat_oil > 0.0) {
                        return satRv(press, temp);
                    } else {
                        return std::min(satRv(press, temp), linearInterpolation(depth_, rv_, depth));
                    }
                }

            private:
                const BlackoilPropertiesInterface& props_;
                const int cell_;
                std::vector<double> depth_; /**< Depth nodes */
                std::vector<double> rv_;    /**< Vaporized oil-gas ratio */
                double z_[BlackoilPhases::MaxNumPhases];
                mutable double A_[BlackoilPhases::MaxNumPhases * BlackoilPhases::MaxNumPhases];

                double satRv(const double press, const double temp) const
                {
                    props_.matrix(1, &press, &temp, z_, &cell_, A_, 0);
                    // Rv/Bg is in the oil row and gas column of A_.
                    // 1/Bg is in the gas row and column.
                    // Recall also that it is stored in column-major order.
                    const int opos = props_.phaseUsage().phase_pos[BlackoilPhases::Liquid];
                    const int gpos = props_.phaseUsage().phase_pos[BlackoilPhases::Vapour];
                    const int np = props_.numPhases();
                    return A_[np*gpos + opos] / A_[np*gpos + gpos];
                }
            };


            /**
             * Class that implements "dissolved gas-oil ratio" (Rs)
             * as function of depth and pressure as follows:
             *
             *   1. The Rs at the gas-oil contact is equal to the
             *      saturated Rs value, Rs_sat_contact.
             *
             *   2. The Rs elsewhere is equal to Rs_sat_contact, but
             *      constrained to the saturated value as given by the
             *      local pressure.
             *
             * This should yield Rs-values that are constant below the
             * contact, and decreasing above the contact.
             */
            class RsSatAtContact : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] props      property object
                 * \param[in] cell       any cell in the pvt region
                 * \param[in] p_contact  oil pressure at the contact
                 * \param[in] T_contact  temperature at the contact
                 */
                RsSatAtContact(const BlackoilPropertiesInterface& props, const int cell, const double p_contact,  const double T_contact)
                    : props_(props), cell_(cell)
                {
                    auto pu = props_.phaseUsage();
                    std::fill(z_, z_ + BlackoilPhases::MaxNumPhases, 0.0);
                    z_[pu.phase_pos[BlackoilPhases::Vapour]] = 1e100;
                    z_[pu.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    rs_sat_contact_ = satRs(p_contact, T_contact);
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
                 * \param[in] temp Temperature at which to calculate RS
                 * value.
                 *
                 * \return Dissolved gas-oil ratio (RS) at depth @c
                 * depth and pressure @c press.
                 */
                double
                operator()(const double /* depth */,
                           const double press,
                           const double temp,
                           const double sat_gas = 0.0) const
                {                     
                    if (sat_gas > 0.0) {
                        return satRs(press, temp);
                    } else {
                        return std::min(satRs(press, temp), rs_sat_contact_);
                    }
                }

            private:
                const BlackoilPropertiesInterface& props_;
                const int cell_;
                double z_[BlackoilPhases::MaxNumPhases];
                double rs_sat_contact_;
                mutable double A_[BlackoilPhases::MaxNumPhases * BlackoilPhases::MaxNumPhases];

                double satRs(const double press, const double temp) const
                {
                    props_.matrix(1, &press, &temp, z_, &cell_, A_, 0);
                    // Rs/Bo is in the gas row and oil column of A_.
                    // 1/Bo is in the oil row and column.
                    // Recall also that it is stored in column-major order.
                    const int opos = props_.phaseUsage().phase_pos[BlackoilPhases::Liquid];
                    const int gpos = props_.phaseUsage().phase_pos[BlackoilPhases::Vapour];
                    const int np = props_.numPhases();
                    return A_[np*opos + gpos] / A_[np*opos + opos];
                }
            };


            /**
             * Class that implements "vaporized oil-gas ratio" (Rv)
             * as function of depth and pressure as follows:
             *
             *   1. The Rv at the gas-oil contact is equal to the
             *      saturated Rv value, Rv_sat_contact.
             *
             *   2. The Rv elsewhere is equal to Rv_sat_contact, but
             *      constrained to the saturated value as given by the
             *      local pressure.
             *
             * This should yield Rv-values that are constant below the
             * contact, and decreasing above the contact.
             */
            class RvSatAtContact : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] props      property object
                 * \param[in] cell       any cell in the pvt region
                 * \param[in] p_contact  oil pressure at the contact
                 * \param[in] T_contact  temperature at the contact
                 */
                RvSatAtContact(const BlackoilPropertiesInterface& props, const int cell, const double p_contact, const double T_contact)
                    : props_(props), cell_(cell)
                {
                    auto pu = props_.phaseUsage();
                    std::fill(z_, z_ + BlackoilPhases::MaxNumPhases, 0.0);
                    z_[pu.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    z_[pu.phase_pos[BlackoilPhases::Liquid]] = 1e100;
                    rv_sat_contact_ = satRv(p_contact, T_contact);
                }

                /**
                 * Function call.
                 *
                 * \param[in] depth Depth at which to calculate RV
                 * value.
                 *
                 * \param[in] press Pressure at which to calculate RV
                 * value.
                 *
                 * \param[in] temp Temperature at which to calculate RV
                 * value.
                 *
                 * \return Dissolved oil-gas ratio (RV) at depth @c
                 * depth and pressure @c press.
                 */
                double
                operator()(const double /*depth*/,
                           const double press,
                           const double temp,
                           const double sat_oil = 0.0) const
                {
                    if (sat_oil > 0.0) {
                        return satRv(press, temp);
                    } else {
                        return std::min(satRv(press, temp), rv_sat_contact_);
                    }
                }

            private:
                const BlackoilPropertiesInterface& props_;
                const int cell_;
                double z_[BlackoilPhases::MaxNumPhases];
                double rv_sat_contact_;
                mutable double A_[BlackoilPhases::MaxNumPhases * BlackoilPhases::MaxNumPhases];

                double satRv(const double press, const double temp) const
                {
                    props_.matrix(1, &press, &temp, z_, &cell_, A_, 0);
                    // Rv/Bg is in the oil row and gas column of A_.
                    // 1/Bg is in the gas row and column.
                    // Recall also that it is stored in column-major order.
                    const int opos = props_.phaseUsage().phase_pos[BlackoilPhases::Liquid];
                    const int gpos = props_.phaseUsage().phase_pos[BlackoilPhases::Vapour];
                    const int np = props_.numPhases();
                    return A_[np*gpos + opos] / A_[np*gpos + gpos];
                }
            };

        } // namespace Miscibility



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
         *
         * For the time being, items 7-9 of ECLIPSE's 'EQUIL' records are also 
         * stored here, but might (should?) eventually be moved elsewhere.
         * 
         *   - @c live_oil_table_index Indicates type of initialisation for live oil.
         *      Positive value points to corresponding Rs vs. depth table.
         *   - @c wet_gas_table_index Indicates type of initialisation for wet gas.
         *      Positive value points to corresponding Rv vs. depth table.
         *   - @c N Defines accuracy of initialisation computations.  Currently
         *      only @c N=0 is supported.
         *      
         */
        struct EquilRecord {
            struct {
                double depth;
                double press;
            } main, woc, goc;
            int live_oil_table_index;
            int wet_gas_table_index;
            int N;
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
         */
        template <class DensCalc>
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
                     std::shared_ptr<Miscibility::RsFunction> rs,
                     std::shared_ptr<Miscibility::RsFunction> rv,
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
            typedef Miscibility::RsFunction CalcDissolution;

            /**
             * Type of vapourised oil-gas ratio calculator.
             */
            typedef Miscibility::RsFunction CalcEvaporation;

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
            dissolutionCalculator() const { return *this->rs_; }

            /**
             * Retrieve vapourised oil-gas ratio calculator of current
             * region.
             */
            const CalcEvaporation&
            evaporationCalculator() const { return *this->rv_; }

            /**
             * Retrieve active fluid phase summary.
             */
            const PhaseUsage&
            phaseUsage() const { return this->pu_; }

        private:
            EquilRecord rec_;     /**< Equilibration data */
            DensCalc    density_; /**< Density calculator */
            std::shared_ptr<Miscibility::RsFunction> rs_;      /**< RS calculator */
            std::shared_ptr<Miscibility::RsFunction> rv_;      /**< RV calculator */
            PhaseUsage  pu_;      /**< Active phase summary */
        };



        /// Functor for inverting capillary pressure function.
        /// Function represented is
        ///   f(s) = pc(s) - target_pc
        struct PcEq
        {
            PcEq(const BlackoilPropertiesInterface& props,
                 const int phase,
                 const int cell,
                 const double target_pc)
                : props_(props),
                  phase_(phase),
                  cell_(cell),
                  target_pc_(target_pc)
            {
                std::fill(s_, s_ + BlackoilPhases::MaxNumPhases, 0.0);
                std::fill(pc_, pc_ + BlackoilPhases::MaxNumPhases, 0.0);
            }
            double operator()(double s) const
            {
                s_[phase_] = s;
                props_.capPress(1, s_, &cell_, pc_, 0);
                return pc_[phase_] - target_pc_;
            }
        private:
            const BlackoilPropertiesInterface& props_;
            const int phase_;
            const int cell_;
            const int target_pc_;
            mutable double s_[BlackoilPhases::MaxNumPhases];
            mutable double pc_[BlackoilPhases::MaxNumPhases];
        };



        /// Compute saturation of some phase corresponding to a given
        /// capillary pressure.
        inline double satFromPc(const BlackoilPropertiesInterface& props,
                                const int phase,
                                const int cell,
                                const double target_pc,
                                const bool increasing = false)
        {
            // Find minimum and maximum saturations.
            double sminarr[BlackoilPhases::MaxNumPhases];
            double smaxarr[BlackoilPhases::MaxNumPhases];
            props.satRange(1, &cell, sminarr, smaxarr);
            const double s0 = increasing ? smaxarr[phase] : sminarr[phase];
            const double s1 = increasing ? sminarr[phase] : smaxarr[phase];

            // Create the equation f(s) = pc(s) - target_pc
            const PcEq f(props, phase, cell, target_pc);
            const double f0 = f(s0);
            const double f1 = f(s1);

            if (f0 <= 0.0) {
                return s0;
            } else if (f1 > 0.0) {
                return s1;
            } else {
                const int max_iter = 30;
                const double tol = 1e-6;
                int iter_used = -1;
                typedef RegulaFalsi<ThrowOnError> ScalarSolver;
                const double sol = ScalarSolver::solve(f, std::min(s0, s1), std::max(s0, s1), max_iter, tol, iter_used);
                return sol;
            }
        }


        /// Functor for inverting a sum of capillary pressure functions.
        /// Function represented is
        ///   f(s) = pc1(s) + pc2(1 - s) - target_pc
        struct PcEqSum
        {
            PcEqSum(const BlackoilPropertiesInterface& props,
                    const int phase1,
                    const int phase2,
                    const int cell,
                    const double target_pc)
                : props_(props),
                  phase1_(phase1),
                  phase2_(phase2),
                  cell_(cell),
                  target_pc_(target_pc)
            {
                std::fill(s_, s_ + BlackoilPhases::MaxNumPhases, 0.0);
                std::fill(pc_, pc_ + BlackoilPhases::MaxNumPhases, 0.0);
            }
            double operator()(double s) const
            {
                s_[phase1_] = s;
                s_[phase2_] = 1.0 - s;
                props_.capPress(1, s_, &cell_, pc_, 0);
                return pc_[phase1_] + pc_[phase2_] - target_pc_;
            }
        private:
            const BlackoilPropertiesInterface& props_;
            const int phase1_;
            const int phase2_;
            const int cell_;
            const int target_pc_;
            mutable double s_[BlackoilPhases::MaxNumPhases];
            mutable double pc_[BlackoilPhases::MaxNumPhases];
        };




        /// Compute saturation of some phase corresponding to a given
        /// capillary pressure, where the capillary pressure function
        /// is given as a sum of two other functions.
        inline double satFromSumOfPcs(const BlackoilPropertiesInterface& props,
                                      const int phase1,
                                      const int phase2,
                                      const int cell,
                                      const double target_pc)
        {
            // Find minimum and maximum saturations.
            double sminarr[BlackoilPhases::MaxNumPhases];
            double smaxarr[BlackoilPhases::MaxNumPhases];
            props.satRange(1, &cell, sminarr, smaxarr);
            const double smin = sminarr[phase1];
            const double smax = smaxarr[phase1];

            // Create the equation f(s) = pc1(s) + pc2(1-s) - target_pc
            const PcEqSum f(props, phase1, phase2, cell, target_pc);
            const double f0 = f(smin);
            const double f1 = f(smax);
            if (f0 <= 0.0) {
                return smin;
            } else if (f1 > 0.0) {
                return smax;
            } else {
                const int max_iter = 30;
                const double tol = 1e-6;
                int iter_used = -1;
                typedef RegulaFalsi<ThrowOnError> ScalarSolver;
                const double sol = ScalarSolver::solve(f, smin, smax, max_iter, tol, iter_used);
                return sol;
            }
        }

        /// Compute saturation from depth. Used for constant capillary pressure function
        inline double satFromDepth(const BlackoilPropertiesInterface& props,
                                   const double cellDepth,
                                   const double contactDepth,
                                   const int phase,
                                   const int cell,
                                   const bool increasing = false)
        {
            // Find minimum and maximum saturations.
            double sminarr[BlackoilPhases::MaxNumPhases];
            double smaxarr[BlackoilPhases::MaxNumPhases];
            props.satRange(1, &cell, sminarr, smaxarr);
            const double s0 = increasing ? smaxarr[phase] : sminarr[phase];
            const double s1 = increasing ? sminarr[phase] : smaxarr[phase];

            if (cellDepth < contactDepth){
                return s0;
            } else {
                return s1;
            }

        }

        /// Return true if capillary pressure function is constant
        bool isConstPc(const BlackoilPropertiesInterface& props,
                       const int phase,
                       const int cell)
        {
            // Find minimum and maximum saturations.
            double sminarr[BlackoilPhases::MaxNumPhases];
            double smaxarr[BlackoilPhases::MaxNumPhases];
            props.satRange(1, &cell, sminarr, smaxarr);

            // Create the equation f(s) = pc(s);
            const PcEq f(props, phase, cell, 0);
            const double f0 = f(sminarr[phase]);
            const double f1 = f(smaxarr[phase]);
            return std::abs(f0 - f1) < std::numeric_limits<double>::epsilon();
        }

    } // namespace Equil
} // namespace Opm


#endif // OPM_EQUILIBRATIONHELPERS_HEADER_INCLUDED
