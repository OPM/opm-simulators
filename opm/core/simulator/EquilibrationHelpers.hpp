/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2017 IRIS

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

#include <opm/common/utility/numeric/linearInterpolation.hpp>
#include <opm/core/utility/RegionMapping.hpp>
#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <memory>


/*
---- synopsis of EquilibrationHelpers.hpp ----

namespace Opm
{
    namespace EQUIL {

        namespace Miscibility {
            class RsFunction;
            class NoMixing;
            template <class FluidSystem>
            class RsVD;
            template <class FluidSystem>
            class RsSatAtContact;
        }

        class EquilReg;


        template <class FluidSystem,  class MaterialLaw, class MaterialLawManager>
        struct PcEq;

        template <class FluidSystem, class MaterialLaw, class MaterialLawManager >
        inline double satFromPc(const MaterialLawManager& materialLawManager,
                                const int phase,
                                const int cell,
                                const double target_pc,
                                const bool increasing = false)
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
        inline double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                                      const int phase1,
                                      const int phase2,
                                      const int cell,
                                      const double target_pc)
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
    namespace EQUIL {


    typedef Opm::FluidSystems::BlackOil<double> FluidSystemSimple;

    // Adjust oil pressure according to gas saturation and cap pressure
    typedef Opm::SimpleModularFluidState<double,
    /*numPhases=*/3,
    /*numComponents=*/3,
    FluidSystemSimple,
    /*storePressure=*/false,
    /*storeTemperature=*/false,
    /*storeComposition=*/false,
    /*storeFugacity=*/false,
    /*storeSaturation=*/true,
    /*storeDensity=*/false,
    /*storeViscosity=*/false,
    /*storeEnthalpy=*/false> SatOnlyFluidState;

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
            template <class FluidSystem>
            class RsVD : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] pvtRegionIdx The pvt region index
                 * \param[in] depth Depth nodes.
                 * \param[in] rs Dissolved gas-oil ratio at @c depth.
                 */
                RsVD(const int pvtRegionIdx,
                     const std::vector<double>& depth,
                     const std::vector<double>& rs)
                    : pvtRegionIdx_(pvtRegionIdx)
                    , depth_(depth)
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
                        return std::min(satRs(press, temp), linearInterpolationNoExtrapolation(depth_, rs_, depth));
                    }
                }

            private:
                const int pvtRegionIdx_;
                std::vector<double> depth_; /**< Depth nodes */
                std::vector<double> rs_;    /**< Dissolved gas-oil ratio */

                double satRs(const double press, const double temp) const
                {
                    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
                }
            };


            /**
             * Type that implements "vaporized oil-gas ratio"
             * tabulated as a function of depth policy.  Data
             * typically taken from keyword 'RVVD'.
             */
            template <class FluidSystem>
            class RvVD : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] pvtRegionIdx The pvt region index
                 * \param[in] depth Depth nodes.
                 * \param[in] rv Dissolved gas-oil ratio at @c depth.
                 */
                RvVD(const int pvtRegionIdx,
                     const std::vector<double>& depth,
                     const std::vector<double>& rv)
                    : pvtRegionIdx_(pvtRegionIdx)
                    , depth_(depth)
                    , rv_(rv)
                {
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
                    if (std::abs(sat_oil) > 1e-16) {
                        return satRv(press, temp);
                    } else {
                        return std::min(satRv(press, temp), linearInterpolationNoExtrapolation(depth_, rv_, depth));
                    }
                }

            private:
                const int pvtRegionIdx_;
                std::vector<double> depth_; /**< Depth nodes */
                std::vector<double> rv_;    /**< Vaporized oil-gas ratio */

                double satRv(const double press, const double temp) const
                {
                    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
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
            template <class FluidSystem>
            class RsSatAtContact : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] pvtRegionIdx The pvt region index
                 * \param[in] p_contact  oil pressure at the contact
                 * \param[in] T_contact  temperature at the contact
                 */
                RsSatAtContact(const int pvtRegionIdx, const double p_contact,  const double T_contact)
                    : pvtRegionIdx_(pvtRegionIdx)
                {
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
                const int pvtRegionIdx_;
                double rs_sat_contact_;

                double satRs(const double press, const double temp) const
                {
                    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
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
            template <class FluidSystem>
            class RvSatAtContact : public RsFunction {
            public:
                /**
                 * Constructor.
                 *
                 * \param[in] pvtRegionIdx The pvt region index
                 * \param[in] p_contact  oil pressure at the contact
                 * \param[in] T_contact  temperature at the contact
                 */
                RvSatAtContact(const int pvtRegionIdx, const double p_contact, const double T_contact)
                    :pvtRegionIdx_(pvtRegionIdx)
                {
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
                const int pvtRegionIdx_;
                double rv_sat_contact_;

                double satRv(const double press, const double temp) const
                {
                    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);;
                }
            };

        } // namespace Miscibility

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
        class EquilReg {
        public:
            /**
             * Constructor.
             *
             * \param[in] rec     Equilibration data of current region.
             * \param[in] rs      Calculator of dissolved gas-oil ratio.
             * \param[in] rv      Calculator of vapourised oil-gas ratio.
             * \param[in] pvtRegionIdx The pvt region index
             */
            EquilReg(const EquilRecord& rec,
                     std::shared_ptr<Miscibility::RsFunction> rs,
                     std::shared_ptr<Miscibility::RsFunction> rv,
                     const int pvtIdx)
                : rec_    (rec)
                , rs_     (rs)
                , rv_     (rv)
                , pvtIdx_ (pvtIdx)
            {
            }

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
            double datum()    const { return this->rec_.datumDepth(); }

            /**
             * Pressure at datum depth in current region.
             */
            double pressure() const { return this->rec_.datumDepthPressure(); }

            /**
             * Depth of water-oil contact.
             */
            double zwoc()     const { return this->rec_.waterOilContactDepth(); }

            /**
             * water-oil capillary pressure at water-oil contact.
             *
             * \return P_o - P_w at WOC.
             */
            double pcow_woc() const { return this->rec_.waterOilContactCapillaryPressure(); }

            /**
             * Depth of gas-oil contact.
             */
            double zgoc()     const { return this->rec_.gasOilContactDepth(); }

            /**
             * Gas-oil capillary pressure at gas-oil contact.
             *
             * \return P_g - P_o at GOC.
             */
            double pcgo_goc() const { return this->rec_.gasOilContactCapillaryPressure(); }


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
             * Retrieve pvtIdx of the region.
             */
            int pvtIdx() const { return this->pvtIdx_; }


        private:
            EquilRecord rec_;     /**< Equilibration data */
            std::shared_ptr<Miscibility::RsFunction> rs_;      /**< RS calculator */
            std::shared_ptr<Miscibility::RsFunction> rv_;      /**< RV calculator */
            const int pvtIdx_;
        };



        /// Functor for inverting capillary pressure function.
        /// Function represented is
        ///   f(s) = pc(s) - target_pc
        template <class FluidSystem,  class MaterialLaw, class MaterialLawManager>
        struct PcEq
        {
            PcEq(const MaterialLawManager& materialLawManager,
                 const int phase,
                 const int cell,
                 const double target_pc)
                : materialLawManager_(materialLawManager),
                  phase_(phase),
                  cell_(cell),
                  target_pc_(target_pc)
            {

            }
            double operator()(double s) const
            {
                const auto& matParams = materialLawManager_.materialLawParams(cell_);
                SatOnlyFluidState fluidState;
                fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
                fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
                fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
                fluidState.setSaturation(phase_, s);

                double pc[FluidSystem::numPhases];
                std::fill(pc, pc + FluidSystem::numPhases, 0.0);
                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                double sign = (phase_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
                double pcPhase = pc[FluidSystem::oilPhaseIdx] + sign *  pc[phase_];
                return pcPhase - target_pc_;
            }
        private:
            const MaterialLawManager& materialLawManager_;
            const int phase_;
            const int cell_;
            const double target_pc_;
        };

        template <class FluidSystem, class MaterialLawManager>
        double minSaturations(const MaterialLawManager& materialLawManager, const int phase, const int cell) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

            // Find minimum and maximum saturations.
            switch(phase) {
            case FluidSystem::waterPhaseIdx :
            {
                 return scaledDrainageInfo.Swl;
                 break;
            }
            case FluidSystem::gasPhaseIdx :
            {
                 return scaledDrainageInfo.Sgl;
                 break;
            }
            case FluidSystem::oilPhaseIdx :
            {
                 OPM_THROW(std::runtime_error, "Min saturation not implemented for oil phase.");
                 break;
            }
            default:  OPM_THROW(std::runtime_error, "Unknown phaseIdx .");
            }
            return -1.0;
        }

        template <class FluidSystem, class MaterialLawManager>
        double maxSaturations(const MaterialLawManager& materialLawManager, const int phase, const int cell) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

            // Find minimum and maximum saturations.
            switch(phase) {
            case FluidSystem::waterPhaseIdx :
            {
                 return scaledDrainageInfo.Swu;
                 break;
            }
            case FluidSystem::gasPhaseIdx :
            {
                 return scaledDrainageInfo.Sgu;
                 break;
            }
            case FluidSystem::oilPhaseIdx :
            {
                 OPM_THROW(std::runtime_error, "Max saturation not implemented for oil phase.");
                 break;
            }
            default:  OPM_THROW(std::runtime_error, "Unknown phaseIdx .");
            }
            return -1.0;
        }


        /// Compute saturation of some phase corresponding to a given
        /// capillary pressure.
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager >
        inline double satFromPc(const MaterialLawManager& materialLawManager,
                                const int phase,
                                const int cell,
                                const double target_pc,
                                const bool increasing = false)
        {
            // Find minimum and maximum saturations.           
            const double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
            const double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

            // Create the equation f(s) = pc(s) - target_pc
            const PcEq<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase, cell, target_pc);
            const double f0 = f(s0);
            const double f1 = f(s1);

            if (f0 <= 0.0) {
                return s0;
            } else if (f1 > 0.0) {
                return s1;
            } else {
                const int max_iter = 60;
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
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
        struct PcEqSum
        {
            PcEqSum(const MaterialLawManager& materialLawManager,
                    const int phase1,
                    const int phase2,
                    const int cell,
                    const double target_pc)
                : materialLawManager_(materialLawManager),
                  phase1_(phase1),
                  phase2_(phase2),
                  cell_(cell),
                  target_pc_(target_pc)
            {
            }
            double operator()(double s) const
            {
                const auto& matParams = materialLawManager_.materialLawParams(cell_);
                SatOnlyFluidState fluidState;
                fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
                fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
                fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
                fluidState.setSaturation(phase1_, s);
                fluidState.setSaturation(phase2_, 1.0 - s);

                double pc[FluidSystem::numPhases];
                std::fill(pc, pc + FluidSystem::numPhases, 0.0);

                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                double sign1 = (phase1_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
                double pc1 = pc[FluidSystem::oilPhaseIdx] + sign1 *  pc[phase1_];
                double sign2 = (phase2_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
                double pc2 = pc[FluidSystem::oilPhaseIdx] + sign2 *  pc[phase2_];
                return pc1 + pc2 - target_pc_;
            }
        private:
            const MaterialLawManager& materialLawManager_;
            const int phase1_;
            const int phase2_;
            const int cell_;
            const double target_pc_;
        };




        /// Compute saturation of some phase corresponding to a given
        /// capillary pressure, where the capillary pressure function
        /// is given as a sum of two other functions.
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
        inline double satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                                      const int phase1,
                                      const int phase2,
                                      const int cell,
                                      const double target_pc)
        {
            // Find minimum and maximum saturations.
            const double smin = minSaturations<FluidSystem>(materialLawManager, phase1, cell);
            const double smax = maxSaturations<FluidSystem>(materialLawManager, phase1, cell);

            // Create the equation f(s) = pc1(s) + pc2(1-s) - target_pc
            const PcEqSum<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase1, phase2, cell, target_pc);
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
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
        inline double satFromDepth(const MaterialLawManager& materialLawManager,
                                   const double cellDepth,
                                   const double contactDepth,
                                   const int phase,
                                   const int cell,
                                   const bool increasing = false)
        {
            const double s0 = increasing ? maxSaturations<FluidSystem>(materialLawManager, phase, cell) : minSaturations<FluidSystem>(materialLawManager, phase, cell);
            const double s1 = increasing ? minSaturations<FluidSystem>(materialLawManager, phase, cell) : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

            if (cellDepth < contactDepth){
                return s0;
            } else {
                return s1;
            }

        }

        /// Return true if capillary pressure function is constant
        template <class FluidSystem, class MaterialLaw, class MaterialLawManager>
        inline bool isConstPc(const MaterialLawManager& materialLawManager,
                              const int                          phase,
                              const int                          cell)
        {
            // Create the equation f(s) = pc(s);
            const PcEq<FluidSystem, MaterialLaw, MaterialLawManager> f(materialLawManager, phase, cell, 0);
            const double f0 = f(minSaturations<FluidSystem>(materialLawManager, phase, cell));
            const double f1 = f(maxSaturations<FluidSystem>(materialLawManager, phase, cell));
            return std::abs(f0 - f1) < std::numeric_limits<double>::epsilon();
        }

    } // namespace Equil
} // namespace Opm


#endif // OPM_EQUILIBRATIONHELPERS_HEADER_INCLUDED
