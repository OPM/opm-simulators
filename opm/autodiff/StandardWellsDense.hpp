/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.
  Copyright 2016 IRIS AS

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


#ifndef OPM_STANDARDWELLSDENSE_HEADER_INCLUDED
#define OPM_STANDARDWELLSDENSE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>
#include <opm/core/wells/WellCollection.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoilDense.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/WellSwitchingLogger.hpp>

namespace Opm {

enum WellVariablePositions {
    XvarWell = 0,
    WFrac = 1,
    GFrac = 2
};


        /// Class for handling the standard well model.
        template<typename FluidSystem, typename BlackoilIndices>
        class StandardWellsDense {
        public:

            // ---------      Types      ---------
            typedef WellStateFullyImplicitBlackoilDense WellState;
            typedef BlackoilModelParameters ModelParameters;

            typedef double Scalar;
            static const int blocksize = 3;
            typedef Dune::FieldVector<Scalar, blocksize    > VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, blocksize, blocksize > MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType> Mat;
            typedef Dune::BlockVector<VectorBlockType> BVector;
            typedef DenseAd::Evaluation<double, /*size=*/blocksize*2> EvalWell;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<BlackoilPropsAdFromDeck::FluidSystem, std::vector<int> >;

            // ---------  Public methods  ---------
            StandardWellsDense(const Wells* wells_arg,
                               WellCollection* well_collection,
                               const ModelParameters& param,
                               const bool terminal_output);

            void init(const PhaseUsage phase_usage_arg,
                      const std::vector<bool>& active_arg,
                      const VFPProperties*  vfp_properties_arg,
                      const double gravity_arg,
                      const std::vector<double>& depth_arg,
                      const std::vector<double>& pv_arg,
                      const RateConverterType* rate_converter);


            template <typename Simulator>
            SimulatorReport assemble(Simulator& ebosSimulator,
                                     const int iterationIdx,
                                     const double dt,
                                     WellState& well_state);

            template <typename Simulator>
            void assembleWellEq(Simulator& ebosSimulator,
                                const double dt,
                                WellState& well_state,
                                bool only_wells);

            template <typename Simulator>
            bool allow_cross_flow(const int w, Simulator& ebosSimulator) const;

            void localInvert(Mat& istlA) const;

            void print(Mat& istlA) const;

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax);

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax);

            // xw = inv(D)*(rw - C*x)
            void recoverVariable(const BVector& x, BVector& xw) const;

            int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

            int flowToEbosPvIdx( const int flowPv ) const;

            int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const;

            int ebosCompToFlowPhaseIdx( const int compIdx ) const;

            std::vector<double>
            extractPerfData(const std::vector<double>& in) const;

            int numPhases() const;

            int numCells() const;

            void resetWellControlFromState(WellState xw);

            const Wells& wells() const;

            const Wells* wellsPointer() const;

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            int numWellVars() const;

            /// Density of each well perforation
            const std::vector<double>& wellPerforationDensities() const;

            /// Diff to bhp for each well perforation.
            const std::vector<double>& wellPerforationPressureDiffs() const;

            typedef DenseAd::Evaluation<double, /*size=*/blocksize> Eval;

            EvalWell extendEval(Eval in) const;

            void setWellVariables(const WellState& xw);

            void print(EvalWell in) const;

            void computeAccumWells();

            template<typename intensiveQuants>
            void
            computeWellFlux(const int& w, const double& Tw, const intensiveQuants& intQuants,
                            const EvalWell& bhp, const double& cdp, const bool& allow_cf, std::vector<EvalWell>& cq_s)  const;

            template <typename Simulator>
            SimulatorReport solveWellEq(Simulator& ebosSimulator,
                                        const double dt,
                                        WellState& well_state);

            void printIf(const int c, const double x, const double y, const double eps, const std::string type) const;

            std::vector<double> residual() const;

            template <typename Simulator>
            bool getWellConvergence(Simulator& ebosSimulator,
                                    const int iteration) const;

            template<typename Simulator>
            void
            computeWellConnectionPressures(const Simulator& ebosSimulator,
                                           const WellState& xw);

            template<typename Simulator>
            void
            computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                        const WellState& xw,
                                                        std::vector<double>& b_perf,
                                                        std::vector<double>& rsmax_perf,
                                                        std::vector<double>& rvmax_perf,
                                                        std::vector<double>& surf_dens_perf) const;

            void updateWellState(const BVector& dwells,
                                 WellState& well_state) const;



            void updateWellControls(WellState& xw) const;

            /// upate the dynamic lists related to economic limits
            void
            updateListEconLimited(const Schedule& schedule,
                                  const int current_step,
                                  const Wells* wells_struct,
                                  const WellState& well_state,
                                  DynamicListEconLimited& list_econ_limited) const;

            void computeWellConnectionDensitesPressures(const WellState& xw,
                                                        const std::vector<double>& b_perf,
                                                        const std::vector<double>& rsmax_perf,
                                                        const std::vector<double>& rvmax_perf,
                                                        const std::vector<double>& surf_dens_perf,
                                                        const std::vector<double>& depth_perf,
                                                        const double grav);


            // TODO: Later we might want to change the function to only handle one well,
            // the requirement for well potential calculation can be based on individual wells.
            // getBhp() will be refactored to reduce the duplication of the code calculating the bhp from THP.
            template<typename Simulator>
            void
            computeWellPotentials(const Simulator& ebosSimulator,
                                  WellState& well_state)  const;


            WellCollection* wellCollection() const
            {
                return well_collection_;
            }





            const std::vector<double>&
            wellPerfEfficiencyFactors() const
            {
                return well_perforation_efficiency_factors_;
            }





            void calculateEfficiencyFactors()
            {
                if ( !localWellsActive() ) {
                    return;
                }

                const int nw = wells().number_of_wells;

                for (int w = 0; w < nw; ++w) {
                    const std::string well_name = wells().name[w];
                    const WellNode& well_node = wellCollection()->findWellNode(well_name);

                    const double well_efficiency_factor = well_node.getAccumulativeEfficiencyFactor();

                    // assign the efficiency factor to each perforation related.
                    for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w + 1]; ++perf) {
                        well_perforation_efficiency_factors_[perf] = well_efficiency_factor;
                    }
                }
            }




            void computeWellVoidageRates(const WellState& well_state,
                                         std::vector<double>& well_voidage_rates,
                                         std::vector<double>& voidage_conversion_coeffs) const
            {
                if ( !localWellsActive() ) {
                    return;
                }
                // TODO: for now, we store the voidage rates for all the production wells.
                // For injection wells, the rates are stored as zero.
                // We only store the conversion coefficients for all the injection wells.
                // Later, more delicate model will be implemented here.
                // And for the moment, group control can only work for serial running.
                const int nw = well_state.numWells();
                const int np = well_state.numPhases();

                // we calculate the voidage rate for each well, that means the sum of all the phases.
                well_voidage_rates.resize(nw, 0);
                // store the conversion coefficients, while only for the use of injection wells.
                voidage_conversion_coeffs.resize(nw * np, 1.0);

                std::vector<double> well_rates(np, 0.0);
                std::vector<double> convert_coeff(np, 1.0);

                for (int w = 0; w < nw; ++w) {
                    const bool is_producer = wells().type[w] == PRODUCER;

                    // not sure necessary to change all the value to be positive
                    if (is_producer) {
                        std::transform(well_state.wellRates().begin() + np * w,
                                       well_state.wellRates().begin() + np * (w + 1),
                                       well_rates.begin(), std::negate<double>());

                        // the average hydrocarbon conditions of the whole field will be used
                        const int fipreg = 0; // Not considering FIP for the moment.

                        rate_converter_->calcCoeff(well_rates, fipreg, convert_coeff);
                        well_voidage_rates[w] = std::inner_product(well_rates.begin(), well_rates.end(),
                                                                   convert_coeff.begin(), 0.0);
                    } else {
                        // TODO: Not sure whether will encounter situation with all zero rates
                        // and whether it will cause problem here.
                        std::copy(well_state.wellRates().begin() + np * w,
                                  well_state.wellRates().begin() + np * (w + 1),
                                  well_rates.begin());
                        // the average hydrocarbon conditions of the whole field will be used
                        const int fipreg = 0; // Not considering FIP for the moment.
                        rate_converter_->calcCoeff(well_rates, fipreg, convert_coeff);
                        std::copy(convert_coeff.begin(), convert_coeff.end(),
                                  voidage_conversion_coeffs.begin() + np * w);
                    }
                }
            }





            void applyVREPGroupControl(WellState& well_state) const
            {
                if ( wellCollection()->havingVREPGroups() ) {
                    std::vector<double> well_voidage_rates;
                    std::vector<double> voidage_conversion_coeffs;
                    computeWellVoidageRates(well_state, well_voidage_rates, voidage_conversion_coeffs);
                    wellCollection()->applyVREPGroupControls(well_voidage_rates, voidage_conversion_coeffs);

                    // for the wells under group control, update the currentControls for the well_state
                    for (const WellNode* well_node : wellCollection()->getLeafNodes()) {
                        if (well_node->isInjector() && !well_node->individualControl()) {
                            const int well_index = well_node->selfIndex();
                            well_state.currentControls()[well_index] = well_node->groupControlIndex();
                        }
                    }
                }
            }



        protected:
            bool wells_active_;
            const Wells*   wells_;

            // Well collection is used to enforce the group control
            WellCollection* well_collection_;

            ModelParameters param_;
            bool terminal_output_;

            PhaseUsage phase_usage_;
            std::vector<bool>  active_;
            const VFPProperties* vfp_properties_;
            double gravity_;
            const RateConverterType* rate_converter_;

            // The efficiency factor for each connection. It is specified based on wells and groups,
            // We calculate the factor for each connection for the computation of contributions to the mass balance equations.
            // By default, they should all be one.
            std::vector<double> well_perforation_efficiency_factors_;
            // the depth of the all the cell centers
            // for standard Wells, it the same with the perforation depth
            std::vector<double> cell_depths_;
            std::vector<double> pv_;

            std::vector<double> well_perforation_densities_;
            std::vector<double> well_perforation_pressure_diffs_;

            std::vector<EvalWell> wellVariables_;
            std::vector<double> F0_;

            Mat duneB_;
            Mat duneC_;
            Mat invDuneD_;

            BVector resWell_;

            mutable BVector Cx_;
            mutable BVector invDrw_;
            mutable BVector scaleAddRes_;

            double dbhpMaxRel() const {return param_.dbhp_max_rel_; }
            double dWellFractionMax() const {return param_.dwell_fraction_max_; }

            // protected methods
            EvalWell getBhp(const int wellIdx) const {
                const WellControls* wc = wells().ctrls[wellIdx];
                if (well_controls_get_current_type(wc) == BHP) {
                    EvalWell bhp = 0.0;
                    const double target_rate = well_controls_get_current_target(wc);
                    bhp.setValue(target_rate);
                    return bhp;
                } else if (well_controls_get_current_type(wc) == THP) {
                    const int control = well_controls_get_current(wc);
                    const double thp = well_controls_get_current_target(wc);
                    const double alq = well_controls_iget_alq(wc, control);
                    const int table_id = well_controls_iget_vfp(wc, control);
                    EvalWell aqua = 0.0;
                    EvalWell liquid = 0.0;
                    EvalWell vapour = 0.0;
                    EvalWell bhp = 0.0;
                    double vfp_ref_depth = 0.0;

                    const Opm::PhaseUsage& pu = phase_usage_;

                    if (active_[ Water ]) {
                        aqua = getQs(wellIdx, pu.phase_pos[ Water]);
                    }
                    if (active_[ Oil ]) {
                        liquid = getQs(wellIdx, pu.phase_pos[ Oil ]);
                    }
                    if (active_[ Gas ]) {
                        vapour = getQs(wellIdx, pu.phase_pos[ Gas ]);
                    }
                    if (wells().type[wellIdx] == INJECTOR) {
                        bhp = vfp_properties_->getInj()->bhp(table_id, aqua, liquid, vapour, thp);
                        vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();
                    } else {
                        bhp = vfp_properties_->getProd()->bhp(table_id, aqua, liquid, vapour, thp, alq);
                        vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();
                    }

                    // pick the density in the top layer
                    const int perf = wells().well_connpos[wellIdx];
                    const double rho = well_perforation_densities_[perf];
                    const double dp = wellhelpers::computeHydrostaticCorrection(wells(), wellIdx, vfp_ref_depth, rho, gravity_);
                    bhp -= dp;
                    return bhp;

                }
                const int nw = wells().number_of_wells;
                return wellVariables_[nw*XvarWell + wellIdx];
            }

            EvalWell getQs(const int wellIdx, const int phaseIdx) const {
                EvalWell qs = 0.0;
                const WellControls* wc = wells().ctrls[wellIdx];
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                const double target_rate = well_controls_get_current_target(wc);

                if (wells().type[wellIdx] == INJECTOR) {
                    const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
                    if (comp_frac == 0.0)
                        return qs;

                    if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                        return wellVariables_[nw*XvarWell + wellIdx];
                    }
                    qs.setValue(target_rate);
                    return qs;
                }

                // Producers
                if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP ) {
                    return wellVariables_[nw*XvarWell + wellIdx] * wellVolumeFractionScaled(wellIdx,phaseIdx);
                }

                if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                    // checking how many phases are included in the rate control
                    // to decide wheter it is a single phase rate control or not
                    const double* distr = well_controls_get_current_distr(wc);
                    int num_phases_under_rate_control = 0;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            num_phases_under_rate_control += 1;
                        }
                    }

                    // there should be at least one phase involved
                    assert(num_phases_under_rate_control > 0);

                    // when it is a single phase rate limit
                    if (num_phases_under_rate_control == 1) {
                        if (distr[phaseIdx] == 1.0) {
                            qs.setValue(target_rate);
                            return qs;
                        }

                        int currentControlIdx = 0;
                        for (int i = 0; i < np; ++i) {
                            currentControlIdx += wells().comp_frac[np*wellIdx + i] * i;
                        }

                        const double eps = 1e-6;
                        if (wellVolumeFractionScaled(wellIdx,currentControlIdx) < eps) {
                            return qs;
                        }
                        return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / wellVolumeFractionScaled(wellIdx,currentControlIdx));
                    }

                    // when it is a combined two phase rate limit, such like LRAT
                    // we neec to calculate the rate for the certain phase
                    if (num_phases_under_rate_control == 2) {
                        EvalWell combined_volume_fraction = 0.;
                        for (int p = 0; p < np; ++p) {
                            if (distr[p] == 1.0) {
                                combined_volume_fraction += wellVolumeFractionScaled(wellIdx, p);
                            }
                        }
                        return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / combined_volume_fraction);
                    }

                    // suppose three phase combined limit is the same with RESV
                    // not tested yet.
                }
                // ReservoirRate
                return target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx);
            }

            EvalWell wellVolumeFraction(const int wellIdx, const int phaseIdx) const {
                const int nw = wells().number_of_wells;
                if (phaseIdx == Water) {
                   return wellVariables_[WFrac * nw + wellIdx];
                }

                if (phaseIdx == Gas) {
                   return wellVariables_[GFrac * nw + wellIdx];
                }

                // Oil fraction
                EvalWell well_fraction = 1.0;
                if (active_[Water]) {
                    well_fraction -= wellVariables_[WFrac * nw + wellIdx];
                }

                if (active_[Gas]) {
                    well_fraction -= wellVariables_[GFrac * nw + wellIdx];
                }
                return well_fraction;
            }

            EvalWell wellVolumeFractionScaled(const int wellIdx, const int phaseIdx) const {
                const WellControls* wc = wells().ctrls[wellIdx];
                if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
                    const double* distr = well_controls_get_current_distr(wc);
                    return wellVolumeFraction(wellIdx, phaseIdx) / distr[phaseIdx];
                }
                std::vector<double> g = {1,1,0.01};
                return (wellVolumeFraction(wellIdx, phaseIdx) / g[phaseIdx]);
            }



            template <class WellState>
            bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                     const WellState& well_state,
                                     const int well_number) const
            {
                const Opm::PhaseUsage& pu = phase_usage_;
                const int np = well_state.numPhases();

                if (econ_production_limits.onMinOilRate()) {
                    assert(active_[Oil]);
                    const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                    const double min_oil_rate = econ_production_limits.minOilRate();
                    if (std::abs(oil_rate) < min_oil_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinGasRate() ) {
                    assert(active_[Gas]);
                    const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
                    const double min_gas_rate = econ_production_limits.minGasRate();
                    if (std::abs(gas_rate) < min_gas_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinLiquidRate() ) {
                    assert(active_[Oil]);
                    assert(active_[Water]);
                    const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                    const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
                    const double liquid_rate = oil_rate + water_rate;
                    const double min_liquid_rate = econ_production_limits.minLiquidRate();
                    if (std::abs(liquid_rate) < min_liquid_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinReservoirFluidRate()) {
                    OpmLog::warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
                }

                return false;
            }


            using WellMapType = typename WellState::WellMapType;
            using WellMapEntryType = typename WellState::mapentry_t;

            // a tuple type for ratio limit check.
            // first value indicates whether ratio limit is violated, when the ratio limit is not violated, the following three
            // values should not be used.
            // second value indicates whehter there is only one connection left.
            // third value indicates the indx of the worst-offending connection.
            // the last value indicates the extent of the violation for the worst-offending connection, which is defined by
            // the ratio of the actual value to the value of the violated limit.
            using RatioCheckTuple = std::tuple<bool, bool, int, double>;

            enum ConnectionIndex {
                INVALIDCONNECTION = -10000
            };


            template <class WellState>
            RatioCheckTuple checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                                 const WellState& well_state,
                                                 const WellMapEntryType& map_entry) const
            {
                // TODO: not sure how to define the worst-offending connection when more than one
                //       ratio related limit is violated.
                //       The defintion used here is that we define the violation extent based on the
                //       ratio between the value and the corresponding limit.
                //       For each violated limit, we decide the worst-offending connection separately.
                //       Among the worst-offending connections, we use the one has the biggest violation
                //       extent.

                bool any_limit_violated = false;
                bool last_connection = false;
                int worst_offending_connection = INVALIDCONNECTION;
                double violation_extent = -1.0;

                if (econ_production_limits.onMaxWaterCut()) {
                    const RatioCheckTuple water_cut_return = checkMaxWaterCutLimit(econ_production_limits, well_state, map_entry);
                    bool water_cut_violated = std::get<0>(water_cut_return);
                    if (water_cut_violated) {
                        any_limit_violated = true;
                        const double violation_extent_water_cut = std::get<3>(water_cut_return);
                        if (violation_extent_water_cut > violation_extent) {
                            violation_extent = violation_extent_water_cut;
                            worst_offending_connection = std::get<2>(water_cut_return);
                            last_connection = std::get<1>(water_cut_return);
                        }
                    }
                }

                if (econ_production_limits.onMaxGasOilRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_GOR", "the support for max Gas-Oil ratio is not implemented yet!");
                }

                if (econ_production_limits.onMaxWaterGasRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_WGR", "the support for max Water-Gas ratio is not implemented yet!");
                }

                if (econ_production_limits.onMaxGasLiquidRatio()) {
                    OpmLog::warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
                }

                if (any_limit_violated) {
                    assert(worst_offending_connection >=0);
                    assert(violation_extent > 1.);
                }

                return std::make_tuple(any_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }

            template <class WellState>
            RatioCheckTuple checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                                  const WellState& well_state,
                                                  const WellMapEntryType& map_entry) const
            {
                bool water_cut_limit_violated = false;
                int worst_offending_connection = INVALIDCONNECTION;
                bool last_connection = false;
                double violation_extent = -1.0;

                const int np = well_state.numPhases();
                const Opm::PhaseUsage& pu = phase_usage_;
                const int well_number = map_entry[0];

                assert(active_[Oil]);
                assert(active_[Water]);

                const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
                const double liquid_rate = oil_rate + water_rate;
                double water_cut;
                if (std::abs(liquid_rate) != 0.) {
                    water_cut = water_rate / liquid_rate;
                } else {
                    water_cut = 0.0;
                }

                const double max_water_cut_limit = econ_production_limits.maxWaterCut();
                if (water_cut > max_water_cut_limit) {
                    water_cut_limit_violated = true;
                }

                if (water_cut_limit_violated) {
                    // need to handle the worst_offending_connection
                    const int perf_start = map_entry[1];
                    const int perf_number = map_entry[2];

                    std::vector<double> water_cut_perf(perf_number);
                    for (int perf = 0; perf < perf_number; ++perf) {
                        const int i_perf = perf_start + perf;
                        const double oil_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
                        const double water_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
                        const double liquid_perf_rate = oil_perf_rate + water_perf_rate;
                        if (std::abs(liquid_perf_rate) != 0.) {
                            water_cut_perf[perf] = water_perf_rate / liquid_perf_rate;
                        } else {
                            water_cut_perf[perf] = 0.;
                        }
                    }

                    last_connection = (perf_number == 1);
                    if (last_connection) {
                        worst_offending_connection = 0;
                        violation_extent = water_cut_perf[0] / max_water_cut_limit;
                        return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
                    }

                    double max_water_cut_perf = 0.;
                    for (int perf = 0; perf < perf_number; ++perf) {
                        if (water_cut_perf[perf] > max_water_cut_perf) {
                            worst_offending_connection = perf;
                            max_water_cut_perf = water_cut_perf[perf];
                        }
                    }

                    assert(max_water_cut_perf != 0.);
                    assert((worst_offending_connection >= 0) && (worst_offending_connection < perf_number));

                    violation_extent = max_water_cut_perf / max_water_cut_limit;
                }

                return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }





            template <class WellState>
            void updateWellStateWithTarget(const WellControls* wc,
                                           const int current,
                                           const int well_index,
                                           WellState& xw) const
            {
                // number of phases
                const int np = wells().number_of_phases;
                // Updating well state and primary variables.
                // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
                const double target = well_controls_iget_target(wc, current);
                const double* distr = well_controls_iget_distr(wc, current);
                switch (well_controls_iget_type(wc, current)) {
                case BHP:
                    xw.bhp()[well_index] = target;
                    break;

                case THP: {
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    const Opm::PhaseUsage& pu = phase_usage_;

                    if (active_[ Water ]) {
                        aqua = xw.wellRates()[well_index*np + pu.phase_pos[ Water ] ];
                    }
                    if (active_[ Oil ]) {
                        liquid = xw.wellRates()[well_index*np + pu.phase_pos[ Oil ] ];
                    }
                    if (active_[ Gas ]) {
                        vapour = xw.wellRates()[well_index*np + pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(wc, current);
                    const double& thp    = well_controls_iget_target(wc, current);
                    const double& alq    = well_controls_iget_alq(wc, current);

                    //Set *BHP* target by calculating bhp from THP
                    const WellType& well_type = wells().type[well_index];

                    // pick the density in the top layer
                    const int perf = wells().well_connpos[well_index];
                    const double rho = well_perforation_densities_[perf];

                    if (well_type == INJECTOR) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), well_index, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                    rho, gravity_);

                        xw.bhp()[well_index] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                    }
                    else if (well_type == PRODUCER) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), well_index, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                    rho, gravity_);

                        xw.bhp()[well_index] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                    }
                    else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }
                    break;
                }

                case RESERVOIR_RATE:
                    // No direct change to any observable quantity at
                    // surface condition.  In this case, use existing
                    // flow rates as initial conditions as reservoir
                    // rate acts only in aggregate.
                    break;

                case SURFACE_RATE:
                    // assign target value as initial guess for injectors and
                    // single phase producers (orat, grat, wrat)
                    const WellType& well_type = wells().type[well_index];
                    if (well_type == INJECTOR) {
                        for (int phase = 0; phase < np; ++phase) {
                            const double& compi = wells().comp_frac[np * well_index + phase];
                            // TODO: it was commented out from the master branch already.
                            //if (compi > 0.0) {
                            xw.wellRates()[np*well_index + phase] = target * compi;
                            //}
                        }
                    } else if (well_type == PRODUCER) {
                        // only set target as initial rates for single phase
                        // producers. (orat, grat and wrat, and not lrat)
                        // lrat will result in numPhasesWithTargetsUnderThisControl == 2
                        int numPhasesWithTargetsUnderThisControl = 0;
                        for (int phase = 0; phase < np; ++phase) {
                            if (distr[phase] > 0.0) {
                                numPhasesWithTargetsUnderThisControl += 1;
                            }
                        }
                        for (int phase = 0; phase < np; ++phase) {
                            if (distr[phase] > 0.0 && numPhasesWithTargetsUnderThisControl < 2 ) {
                                xw.wellRates()[np*well_index + phase] = target * distr[phase];
                            }
                        }
                    } else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }

                    break;
                } // end of switch


                std::vector<double> g = {1.0, 1.0, 0.01};
                if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                    for (int phase = 0; phase < np; ++phase) {
                        g[phase] = distr[phase];
                    }
                }

                // the number of wells
                const int nw = wells().number_of_wells;

                switch (well_controls_iget_type(wc, current)) {
                case THP:
                case BHP: {
                    const WellType& well_type = wells().type[well_index];
                    xw.wellSolutions()[nw*XvarWell + well_index] = 0.0;
                    if (well_type == INJECTOR) {
                        for (int p = 0; p < np; ++p) {
                            xw.wellSolutions()[nw*XvarWell + well_index] += xw.wellRates()[np*well_index + p] * wells().comp_frac[np*well_index + p];
                        }
                    } else {
                        for (int p = 0; p < np; ++p) {
                            xw.wellSolutions()[nw*XvarWell + well_index] += g[p] * xw.wellRates()[np*well_index + p];
                        }
                    }
                    break;
                }
                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                    xw.wellSolutions()[nw*XvarWell + well_index] = xw.bhp()[well_index];
                    break;
                } // end of switch

                double tot_well_rate = 0.0;
                for (int p = 0; p < np; ++p)  {
                    tot_well_rate += g[p] * xw.wellRates()[np*well_index + p];
                }
                if(std::abs(tot_well_rate) > 0) {
                    if (active_[ Water ]) {
                        xw.wellSolutions()[WFrac*nw + well_index] = g[Water] * xw.wellRates()[np*well_index + Water] / tot_well_rate;
                    }
                    if (active_[ Gas ]) {
                        xw.wellSolutions()[GFrac*nw + well_index] = g[Gas] * xw.wellRates()[np*well_index + Gas] / tot_well_rate ;
                    }
                 } else {
                    if (active_[ Water ]) {
                        xw.wellSolutions()[WFrac*nw + well_index] =  wells().comp_frac[np*well_index + Water];
                    }

                    if (active_[ Gas ]) {
                        xw.wellSolutions()[GFrac*nw + well_index] =  wells().comp_frac[np*well_index + Gas];
                    }
                }

            }

        };


} // namespace Opm

#include "StandardWellsDense_impl.hpp"
#endif
