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
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/LinearisedBlackoilResidual.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoilDense.hpp>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

namespace Opm {


        /// Class for handling the standard well model.
        template<typename FluidSystem, typename BlackoilIndices>
        class StandardWellsDense {
        public:

            // ---------      Types      ---------
            typedef DenseAd::Evaluation<double, /*size=*/6> EvalWell;
            typedef WellStateFullyImplicitBlackoilDense WellState;
            typedef BlackoilModelParameters ModelParameters;

            typedef double Scalar;
            typedef Dune::FieldVector<Scalar, 3    > VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, 3, 3 > MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType> Mat;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            // ---------  Public methods  ---------
            StandardWellsDense(const Wells* wells_arg,
                               const ModelParameters& param,
                               const bool terminal_output,
                               const std::vector<double>& pv)
                : wells_active_(wells_arg!=nullptr)
                , wells_(wells_arg)
                , param_(param)
                , terminal_output_(terminal_output)
                , fluid_(nullptr)
                , active_(nullptr)
                , vfp_properties_(nullptr)
                , well_perforation_densities_(wells_arg->well_connpos[wells_arg->number_of_wells])
                , well_perforation_pressure_diffs_(wells_arg->well_connpos[wells_arg->number_of_wells])
                , wellVariables_(wells_arg->number_of_wells * wells_arg->number_of_phases)
                , F0_(wells_arg->number_of_wells * wells_arg->number_of_phases)
              {
                invDuneD_.setBuildMode( Mat::row_wise );
                duneC_.setBuildMode( Mat::row_wise );
                duneB_.setBuildMode( Mat::row_wise );
              }

            void init(const BlackoilPropsAdInterface* fluid_arg,
                      const std::vector<bool>* active_arg,
                      const VFPProperties*  vfp_properties_arg,
                      const double gravity_arg,
                      const std::vector<double>& depth_arg,
                      const std::vector<double>& pv_arg)
            {

                if ( ! localWellsActive() ) {
                    return;
                }

                fluid_ = fluid_arg;
                active_ = active_arg;
                vfp_properties_ = vfp_properties_arg;
                gravity_ = gravity_arg;
                cell_depths_ = extractPerfData(depth_arg);
                pv_ = pv_arg;

                // setup sparsity pattern for the matrices
                //[A B^T    [x    =  [ res
                // C D] x_well]      res_well]

                const int nw = wells().number_of_wells;
                const int nperf = wells().well_connpos[nw];
                const int nc = numCells();

                // set invDuneD
                invDuneD_.setSize( nw, nw, nw );

                // set duneC
                duneC_.setSize( nw, nc, nperf );

                // set duneB
                duneB_.setSize( nw, nc, nperf );

                for(auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
                    // Add nonzeros for diagonal
                    row.insert(row.index());
                }

                for(auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
                    // Add nonzeros for diagonal
                    for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                        const int cell_idx = wells().well_cells[perf];
                        row.insert(cell_idx);
                    }
                }

                // make the B^T matrix
                for(auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
                    for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                        const int cell_idx = wells().well_cells[perf];
                        row.insert(cell_idx);
                    }
                }

                resWell_.resize( nw );
            }


            template <typename Simulator>
            IterationReport assemble(Simulator& ebosSimulator,
                                     const int iterationIdx,
                                     const double dt,
                                     WellState& well_state) {

                IterationReport iter_report = {false, false, 0, 0};
                if ( ! localWellsActive() ) {
                    return iter_report;
                }

                resetWellControlFromState(well_state);
                updateWellControls(well_state);
                // Set the primary variables for the wells
                setWellVariables(well_state);

                if (iterationIdx == 0) {
                    computeWellConnectionPressures(ebosSimulator, well_state);
                    computeAccumWells();
                }

                if (param_.solve_welleq_initially_ && iterationIdx == 0) {
                    // solve the well equations as a pre-processing step
                    iter_report = solveWellEq(ebosSimulator, dt, well_state);
                }
                assembleWellEq(ebosSimulator, dt, well_state, false);


                if (param_.compute_well_potentials_) {
                    //wellModel().computeWellPotentials(mob_perfcells, b_perfcells, state0, well_state);
                }
                return iter_report;
            }

            template <typename Simulator>
            void assembleWellEq(Simulator& ebosSimulator,
                                const double dt,
                                WellState& well_state,
                                bool only_wells) {
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;


                // clear all entries
                duneB_ = 0.0;
                duneC_ = 0.0;
                invDuneD_ = 0.0;
                resWell_ = 0.0;

                auto& ebosJac = ebosSimulator.model().linearizer().matrix();
                auto& ebosResid = ebosSimulator.model().linearizer().residual();

                const double volume = 0.002831684659200; // 0.1 cu ft;
                for (int w = 0; w < nw; ++w) {
                    bool allow_cf = allow_cross_flow(w, ebosSimulator);
                    for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {

                        const int cell_idx = wells().well_cells[perf];
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        std::vector<EvalWell> cq_s(np,0.0);
                        computeWellFlux(w, wells().WI[perf], intQuants, wellPerforationPressureDiffs()[perf], allow_cf, cq_s);

                        for (int p1 = 0; p1 < np; ++p1) {

                            if (!only_wells) {
                                // subtract sum of phase fluxes in the reservoir equation.
                                ebosResid[cell_idx][flowPhaseToEbosCompIdx(p1)] -= cq_s[p1].value();
                            }

                            // subtract sum of phase fluxes in the well equations.
                            resWell_[w][flowPhaseToEbosCompIdx(p1)] -= cq_s[p1].value();

                            // assemble the jacobians
                            for (int p2 = 0; p2 < np; ++p2) {
                                if (!only_wells) {
                                    ebosJac[cell_idx][cell_idx][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s[p1].derivative(p2);
                                    duneB_[w][cell_idx][flowToEbosPvIdx(p2)][flowPhaseToEbosCompIdx(p1)] -= cq_s[p1].derivative(p2+3); // intput in transformed matrix
                                    duneC_[w][cell_idx][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s[p1].derivative(p2);
                                }
                                invDuneD_[w][w][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s[p1].derivative(p2+3);
                            }

                            // Store the perforation phase flux for later usage.
                            well_state.perfPhaseRates()[perf*np + p1] = cq_s[p1].value();
                        }
                        // Store the perforation pressure for later usage.
                        well_state.perfPress()[perf] = well_state.bhp()[w] + wellPerforationPressureDiffs()[perf];
                    }

                    // add vol * dF/dt + Q to the well equations;
                    for (int p1 = 0; p1 < np; ++p1) {
                        EvalWell resWell_loc = (wellVolumeFraction(w, p1) - F0_[w + nw*p1]) * volume / dt;
                        resWell_loc += getQs(w, p1);
                        for (int p2 = 0; p2 < np; ++p2) {
                            invDuneD_[w][w][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] += resWell_loc.derivative(p2+3);
                        }
                        resWell_[w][flowPhaseToEbosCompIdx(p1)] += resWell_loc.value();
                    }
                }

                // do the local inversion of D.
                localInvert( invDuneD_ );

            }

            template <typename Simulator>
            bool allow_cross_flow(const int w, Simulator& ebosSimulator) const {

                if (wells().allow_cf[w]) {
                    return true;
                }
                // check for special case where all perforations have cross flow
                // then the wells must allow for cross flow
                for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    const int cell_idx = wells().well_cells[perf];
                    const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();
                    EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
                    EvalWell bhp = getBhp(w);
                    // Pressure drawdown (also used to determine direction of flow)
                    EvalWell well_pressure = bhp + wellPerforationPressureDiffs()[perf];
                    EvalWell drawdown = pressure - well_pressure;

                    if ( drawdown.value() < 0 && wells().type[w] == INJECTOR)  {
                        return false;
                    }
                    if ( drawdown.value() > 0 && wells().type[w] == PRODUCER)  {
                        return false;
                    }
                }
                return true;
            }

            void localInvert(Mat& istlA) const {
                for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
                    for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                        //std::cout << (*col) << std::endl;
                        (*col).invert();
                    }
                }
            }

            void print(Mat& istlA) const {
                for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
                    for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                        std::cout << row.index() << " " << col.index() << "/n \n"<<(*col) << std::endl;
                    }
                }
            }

            // substract Binv(D)rw from r;
            void apply( BVector& r) const {
                if ( ! localWellsActive() ) {
                    return;
                }

                BVector invDrw(invDuneD_.N());
                invDuneD_.mv(resWell_,invDrw);
                duneB_.mmtv(invDrw, r);
            }

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) {
                if ( ! localWellsActive() ) {
                    return;
                }
                BVector Cx(duneC_.N());
                duneC_.mv(x, Cx);
                BVector invDCx(invDuneD_.N());
                invDuneD_.mv(Cx, invDCx);
                duneB_.mmtv(invDCx,Ax);
            }

            // xw = inv(D)*(rw - C*x)
            void recoverVariable(const BVector& x, BVector& xw) const {
                if ( ! localWellsActive() ) {
                    return;
                }
                BVector resWell = resWell_;
                duneC_.mmv(x, resWell);
                invDuneD_.mv(resWell, xw);
            }


            int flowPhaseToEbosCompIdx( const int phaseIdx ) const
            {
                const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx };
                return phaseToComp[ phaseIdx ];
            }

            int flowToEbosPvIdx( const int flowPv ) const
            {
                const int flowToEbos[ 3 ] = {
                                              BlackoilIndices::pressureSwitchIdx,
                                              BlackoilIndices::waterSaturationIdx,
                                              BlackoilIndices::compositionSwitchIdx
                                            };
                return flowToEbos[ flowPv ];
            }

            int ebosCompToFlowPhaseIdx( const int compIdx ) const
            {
                const int compToPhase[ 3 ] = { Oil, Water, Gas };
                return compToPhase[ compIdx ];
            }



            std::vector<double>
            extractPerfData(const std::vector<double>& in) const {
                const int nw   = wells().number_of_wells;
                const int nperf = wells().well_connpos[nw];
                std::vector<double> out(nperf);
                for (int w = 0; w < nw; ++w) {
                    for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                        const int well_idx = wells().well_cells[perf];
                        out[perf] = in[well_idx];
                    }
                }
                return out;
            }

            int numPhases() const { return wells().number_of_phases; }

            int numCells() const { return pv_.size();}

            template<class WellState>
            void resetWellControlFromState(WellState xw) {
                const int        nw   = wells_->number_of_wells;
                for (int w = 0; w < nw; ++w) {
                    WellControls* wc = wells_->ctrls[w];
                    well_controls_set_current( wc, xw.currentControls()[w]);
                }
            }

            const Wells& wells() const
            {
                assert(wells_ != 0);
                return *(wells_);
            }

            const Wells* wellsPointer() const
            {
                return wells_;
            }

            /// return true if wells are available in the reservoir
            bool wellsActive() const
            {
                return wells_active_;
            }
            void setWellsActive(const bool wells_active)
            {
                wells_active_ = wells_active;
            }
            /// return true if wells are available on this process
            bool localWellsActive() const
            {
                return wells_ ? (wells_->number_of_wells > 0 ) : false;
            }

            int numWellVars() const
            {
                if ( !localWellsActive() )
                {
                    return 0;
                }

                // For each well, we have a bhp variable, and one flux per phase.
                const int nw = wells().number_of_wells;
                return (numPhases() + 1) * nw;
            }

            /// Density of each well perforation
            std::vector<double> wellPerforationDensities() const     {
                return well_perforation_densities_;
            }

            /// Diff to bhp for each well perforation.
            std::vector<double> wellPerforationPressureDiffs() const
            {
                return well_perforation_pressure_diffs_;
            }


            typedef DenseAd::Evaluation<double, /*size=*/3> Eval;
            EvalWell extendEval(Eval in) const {
                EvalWell out = 0.0;
                out.setValue(in.value());
                for(int i = 0;i<3;++i) {
                    out.setDerivative(i, in.derivative(flowToEbosPvIdx(i)));
                }
                return out;
            }

            void
            setWellVariables(const WellState& xw) {
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    for (int w = 0; w < nw; ++w) {
                        wellVariables_[w + nw*phaseIdx] = 0.0;
                        wellVariables_[w + nw*phaseIdx].setValue(xw.wellSolutions()[w + nw* phaseIdx]);
                        wellVariables_[w + nw*phaseIdx].setDerivative(np + phaseIdx, 1.0);
                    }
                }
            }

            void print(EvalWell in) const {
                std::cout << in.value() << std::endl;
                for (int i = 0; i < in.size; ++i) {
                    std::cout << in.derivative(i) << std::endl;
                }
            }

            void
            computeAccumWells() {
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    for (int w = 0; w < nw; ++w) {
                        F0_[w + nw * phaseIdx] = wellVolumeFraction(w,phaseIdx).value();
                    }
                }
            }



            template<typename intensiveQuants>
            void
            computeWellFlux(const int& w, const double& Tw, const intensiveQuants& intQuants, const double& cdp, const bool& allow_cf, std::vector<EvalWell>& cq_s)  const
            {
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                EvalWell bhp = getBhp(w);
                const int np = wells().number_of_phases;
                std::vector<EvalWell> cmix_s(np,0.0);
                for (int phase = 0; phase < np; ++phase) {
                    //int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                    cmix_s[phase] = wellVolumeFraction(w,phase);
                }
                const auto& fs = intQuants.fluidState();
                EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
                EvalWell rs = extendEval(fs.Rs());
                EvalWell rv = extendEval(fs.Rv());
                std::vector<EvalWell> b_perfcells_dense(np, 0.0);
                std::vector<EvalWell> mob_perfcells_dense(np, 0.0);
                for (int phase = 0; phase < np; ++phase) {
                    int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                    b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
                    mob_perfcells_dense[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
                }

                // Pressure drawdown (also used to determine direction of flow)
                EvalWell well_pressure = bhp + cdp;
                EvalWell drawdown = pressure - well_pressure;

                // injection perforations
                if ( drawdown.value() > 0 )  {

                    //Do nothing if crossflow is not allowed
                    if (!allow_cf && wells().type[w] == INJECTOR)
                        return;
                    // compute phase volumetric rates at standard conditions
                    std::vector<EvalWell> cq_ps(np, 0.0);
                    for (int phase = 0; phase < np; ++phase) {
                        const EvalWell cq_p = - Tw * (mob_perfcells_dense[phase] * drawdown);
                        cq_ps[phase] = b_perfcells_dense[phase] * cq_p;
                    }

                    if ((*active_)[Oil] && (*active_)[Gas]) {
                        const int oilpos = pu.phase_pos[Oil];
                        const int gaspos = pu.phase_pos[Gas];
                        const EvalWell cq_psOil = cq_ps[oilpos];
                        const EvalWell cq_psGas = cq_ps[gaspos];
                        cq_ps[gaspos] += rs * cq_psOil;
                        cq_ps[oilpos] += rv * cq_psGas;
                    }

                    // map to ADB
                    for (int phase = 0; phase < np; ++phase) {
                        cq_s[phase] = cq_ps[phase];
                    }

                } else {
                    //Do nothing if crossflow is not allowed
                    if (!allow_cf && wells().type[w] == PRODUCER)
                        return;

                    // Using total mobilities
                    EvalWell total_mob_dense = mob_perfcells_dense[0];
                    for (int phase = 1; phase < np; ++phase) {
                        total_mob_dense += mob_perfcells_dense[phase];
                    }
                    // injection perforations total volume rates
                    const EvalWell cqt_i = - Tw * (total_mob_dense * drawdown);

                    // compute volume ratio between connection at standard conditions
                    EvalWell volumeRatio = 0.0;
                    if ((*active_)[Water]) {
                        const int watpos = pu.phase_pos[Water];
                        volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
                    }

                    if ((*active_)[Oil] && (*active_)[Gas]) {
                        EvalWell well_temperature = extendEval(fs.temperature(FluidSystem::oilPhaseIdx));
                        EvalWell rsSatEval = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);
                        EvalWell rvSatEval = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);

                        const int oilpos = pu.phase_pos[Oil];
                        const int gaspos = pu.phase_pos[Gas];
                        EvalWell rvPerf = 0.0;
                        if (cmix_s[gaspos] > 0)
                            rvPerf = cmix_s[oilpos] / cmix_s[gaspos];

                        if (rvPerf.value() > rvSatEval.value()) {
                            rvPerf = rvSatEval;
                            //rvPerf.setValue(rvSatEval.value());
                        }

                        EvalWell rsPerf = 0.0;
                        if (cmix_s[oilpos] > 0)
                            rsPerf = cmix_s[gaspos] / cmix_s[oilpos];

                        if (rsPerf.value() > rsSatEval.value()) {
                            //rsPerf = 0.0;
                            rsPerf= rsSatEval;
                        }

                        // Incorporate RS/RV factors if both oil and gas active
                        const EvalWell d = 1.0 - rvPerf * rsPerf;

                        const EvalWell tmp_oil = (cmix_s[oilpos] - rvPerf * cmix_s[gaspos]) / d;
                        //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                        volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                        const EvalWell tmp_gas = (cmix_s[gaspos] - rsPerf * cmix_s[oilpos]) / d;
                        //std::cout << "tmp_gas " <<tmp_gas << std::endl;
                        volumeRatio += tmp_gas / b_perfcells_dense[gaspos];
                    }
                    else {
                        if ((*active_)[Oil]) {
                            const int oilpos = pu.phase_pos[Oil];
                            volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                        }
                        if ((*active_)[Gas]) {
                            const int gaspos = pu.phase_pos[Gas];
                            volumeRatio += cmix_s[gaspos] / b_perfcells_dense[gaspos];
                        }
                    }
                    // injecting connections total volumerates at standard conditions
                    EvalWell cqt_is = cqt_i/volumeRatio;
                    //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
                    for (int phase = 0; phase < np; ++phase) {
                        cq_s[phase] = cmix_s[phase] * cqt_is; // * b_perfcells_dense[phase];
                    }
                }
            }


            template <typename Simulator>
            IterationReport solveWellEq(Simulator& ebosSimulator,
                                        const double dt,
                                        WellState& well_state)
            {
                const int nw = wells().number_of_wells;
                WellState well_state0 = well_state;

                int it  = 0;
                bool converged;
                do {
                    assembleWellEq(ebosSimulator, dt, well_state, true);
                    converged = getWellConvergence(ebosSimulator, it);

                    if (converged) {
                        break;
                    }

                    ++it;
                    if( localWellsActive() )
                    {
                        BVector dx_well (nw);
                        invDuneD_.mv(resWell_, dx_well);

                        updateWellState(dx_well, well_state);
                        updateWellControls(well_state);
                        setWellVariables(well_state);
                    }
                } while (it < 15);

                if (!converged) {
                    well_state = well_state0;
                }

                const bool failed = false; // Not needed in this method.
                const int linear_iters = 0; // Not needed in this method
                return IterationReport{failed, converged, linear_iters, it};
            }

            void printIf(int c, double x, double y, double eps, std::string type) {
                if (std::abs(x-y) > eps) {
                    std::cout << type << " " << c << ": "<<x << " " << y << std::endl;
                }
            }


            std::vector<double> residual() {

                const int np = numPhases();
                const int nw = wells().number_of_wells;
                std::vector<double> res(np*nw);
                for( int p=0; p<np; ++p) {
                    const int ebosCompIdx = flowPhaseToEbosCompIdx(p);
                    for (int i = 0; i < nw; ++i) {
                        int idx = i + nw*p;
                        res[idx] = resWell_[ i ][ ebosCompIdx ];
                    }
                }
                return res;
            }


            template <typename Simulator>
            bool getWellConvergence(Simulator& ebosSimulator,
                                    const int iteration)
            {
                typedef std::vector< double > Vector;
                const int np = numPhases();
                const int nc = numCells();
                const double tol_wells = param_.tolerance_wells_;
                const double maxResidualAllowed = param_.max_residual_allowed_;

                Vector R_sum(np);
                Vector B_avg(np);
                Vector maxCoeff(np);
                Vector maxNormWell(np);

                std::vector< Vector > B( np, Vector( nc ) );
                std::vector< Vector > R2( np, Vector( nc ) );
                std::vector< Vector > tempV( np, Vector( nc ) );

                for ( int idx = 0; idx < np; ++idx )
                {
                    Vector& B_idx  = B[ idx ];
                    const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);

                    for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        const auto& fs = intQuants.fluidState();

                        B_idx [cell_idx] = 1 / fs.invB(ebosPhaseIdx).value();
                    }
                }

                detail::convergenceReduction(B, tempV, R2,
                                             R_sum, maxCoeff, B_avg, maxNormWell,
                                             nc, np, pv_, residual() );


                Vector well_flux_residual(np);

                bool converged_Well = true;
                // Finish computation
                for ( int idx = 0; idx < np; ++idx )
                {
                    well_flux_residual[idx] = B_avg[idx] * maxNormWell[idx];
                    converged_Well = converged_Well && (well_flux_residual[idx] < tol_wells);
                }

                // if one of the residuals is NaN, throw exception, so that the solver can be restarted
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                    if (std::isnan(well_flux_residual[phaseIdx])) {
                        OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                    }
                    if (well_flux_residual[phaseIdx] > maxResidualAllowed) {
                        OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                    }
                }

                if ( terminal_output_ )
                {
                    // Only rank 0 does print to std::cout
                    if (iteration == 0) {
                        std::string msg;
                        msg = "Iter";
                        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                            const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                            msg += "  W-FLUX(" + phaseName + ")";
                        }
                        OpmLog::note(msg);
                    }
                    std::ostringstream ss;
                    const std::streamsize oprec = ss.precision(3);
                    const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                    ss << std::setw(4) << iteration;
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        ss << std::setw(11) << well_flux_residual[phaseIdx];
                    }
                    ss.precision(oprec);
                    ss.flags(oflags);
                    OpmLog::note(ss.str());
                }
                return converged_Well;
            }



            template<typename Simulator>
            void
            computeWellConnectionPressures(const Simulator& ebosSimulator,
                                           const WellState& xw)
            {
                if( ! localWellsActive() ) return ;
                // 1. Compute properties required by computeConnectionPressureDelta().
                //    Note that some of the complexity of this part is due to the function
                //    taking std::vector<double> arguments, and not Eigen objects.
                std::vector<double> b_perf;
                std::vector<double> rsmax_perf;
                std::vector<double> rvmax_perf;
                std::vector<double> surf_dens_perf;
                computePropertiesForWellConnectionPressures(ebosSimulator, xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);
                computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, cell_depths_, gravity_);

            }

            template<typename Simulator, class WellState>
            void
            computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                        const WellState& xw,
                                                        std::vector<double>& b_perf,
                                                        std::vector<double>& rsmax_perf,
                                                        std::vector<double>& rvmax_perf,
                                                        std::vector<double>& surf_dens_perf)
            {
                const int nperf = wells().well_connpos[wells().number_of_wells];
                const int nw = wells().number_of_wells;
                const PhaseUsage& pu = fluid_->phaseUsage();
                const int np = fluid_->numPhases();
                b_perf.resize(nperf*np);
                rsmax_perf.resize(nperf);
                rvmax_perf.resize(nperf);
                surf_dens_perf.resize(nperf*np);

                // Compute the average pressure in each well block
                for (int w = 0; w < nw; ++w) {
                    for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {

                        const int cell_idx = wells().well_cells[perf];
                        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                        const auto& fs = intQuants.fluidState();

                        const double p_above = perf == wells().well_connpos[w] ? xw.bhp()[w] : xw.perfPress()[perf - 1];
                        const double p_avg = (xw.perfPress()[perf] + p_above)/2;
                        double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value();

                        if (pu.phase_used[BlackoilPhases::Aqua]) {
                            b_perf[ pu.phase_pos[BlackoilPhases::Aqua] + perf * pu.num_phases] =
                                    FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        }

                        if (pu.phase_used[BlackoilPhases::Vapour]) {
                            int gaspos = pu.phase_pos[BlackoilPhases::Vapour] + perf * pu.num_phases;
                            int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;

                            if (pu.phase_used[BlackoilPhases::Liquid]) {
                                int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
                                const double oilrate = std::abs(xw.wellRates()[oilpos_well]); //in order to handle negative rates in producers
                                rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), temperature, p_avg);
                                if (oilrate > 0) {
                                    const double gasrate = std::abs(xw.wellRates()[gaspos_well]);
                                    double rv = 0.0;
                                    if (gasrate > 0) {
                                        rv = oilrate / gasrate;
                                    }
                                    rv = std::min(rv, rvmax_perf[perf]);

                                    b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rv);
                                }
                                else {
                                    b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                                }

                            } else {
                                b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                            }
                        }

                        if (pu.phase_used[BlackoilPhases::Liquid]) {
                            int oilpos = pu.phase_pos[BlackoilPhases::Liquid] + perf * pu.num_phases;
                            int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
                            if (pu.phase_used[BlackoilPhases::Vapour]) {
                                rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                                int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;
                                const double gasrate = std::abs(xw.wellRates()[gaspos_well]);
                                if (gasrate > 0) {
                                    const double oilrate = std::abs(xw.wellRates()[oilpos_well]);
                                    double rs = 0.0;
                                    if (oilrate > 0) {
                                        rs = gasrate / oilrate;
                                    }
                                    rs = std::min(rs, rsmax_perf[perf]);
                                    b_perf[oilpos] = FluidSystem::oilPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rs);
                                } else {
                                    b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                                }
                            } else {
                                b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                            }
                        }

                        // Surface density.
                        for (int p = 0; p < pu.num_phases; ++p) {
                            surf_dens_perf[np*perf + p] = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( p ), fs.pvtRegionIndex());
                        }
                    }
                }
            }

            template <class WellState>
            void updateWellState(const BVector& dwells,
                                 WellState& well_state)
            {
                if( localWellsActive() )
                {
                    const int np = wells().number_of_phases;
                    const int nw = wells().number_of_wells;

                    double dFLimit = 0.2;
                    double dBHPLimit = 2.0;
                    std::vector<double> xvar_well_old = well_state.wellSolutions();

                    for (int w = 0; w < nw; ++w) {

                        // update the second and third well variable (The flux fractions)
                        std::vector<double> F(np,0.0);
                        const int sign2 = dwells[w][flowPhaseToEbosCompIdx(1)] > 0 ? 1: -1;
                        const double dx2_limited = sign2 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(1)]),dFLimit);
                        well_state.wellSolutions()[nw + w] = xvar_well_old[nw + w] - dx2_limited;
                        const int sign3 = dwells[w][flowPhaseToEbosCompIdx(2)] > 0 ? 1: -1;
                        const double dx3_limited = sign3 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(2)]),dFLimit);
                        well_state.wellSolutions()[2*nw + w] = xvar_well_old[2*nw + w] - dx3_limited;
                        F[Water] = well_state.wellSolutions()[nw + w];
                        F[Gas] = well_state.wellSolutions()[2*nw + w];
                        F[Oil] = 1.0 - F[Water] - F[Gas];

                        if (F[Water] < 0.0) {
                            F[Gas] /= (1.0 - F[Water]);
                            F[Oil] /= (1.0 - F[Water]);
                            F[Water] = 0.0;
                        }
                        if (F[Gas] < 0.0) {
                            F[Water] /= (1.0 - F[Gas]);
                            F[Oil] /= (1.0 - F[Gas]);
                            F[Gas] = 0.0;
                        }
                        if (F[Oil] < 0.0) {
                            F[Water] /= (1.0 - F[Oil]);
                            F[Gas] /= (1.0 - F[Oil]);
                            F[Oil] = 0.0;
                        }
                        well_state.wellSolutions()[nw + w] = F[Water];
                        well_state.wellSolutions()[2*nw + w] = F[Gas];
                        //std::cout << wells().name[w] << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;

                        // The interpretation of the first well variable depends on the well control
                        const WellControls* wc = wells().ctrls[w];

                        // The current control in the well state overrides
                        // the current control set in the Wells struct, which
                        // is instead treated as a default.
                        const int current = well_state.currentControls()[w];
                        const double target_rate = well_controls_iget_target(wc, current);

                        std::vector<double> g = {1,1,0.01};
                        if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                            const double* distr = well_controls_iget_distr(wc, current);
                            for (int p = 0; p < np; ++p) {
                                F[p] /= distr[p];
                            }
                        } else {
                            for (int p = 0; p < np; ++p) {
                                F[p] /= g[p];
                            }
                        }

                        switch (well_controls_iget_type(wc, current)) {
                        case THP: // The BHP and THP both uses the total rate as first well variable.
                        case BHP:
                        {
                            well_state.wellSolutions()[w] = xvar_well_old[w] - dwells[w][flowPhaseToEbosCompIdx(0)];

                            switch (wells().type[w]) {
                            case INJECTOR:
                                for (int p = 0; p < np; ++p) {
                                    const double comp_frac = wells().comp_frac[np*w + p];
                                    well_state.wellRates()[w*np + p] = comp_frac * well_state.wellSolutions()[w];
                                }
                                break;
                            case PRODUCER:
                                for (int p = 0; p < np; ++p) {
                                    well_state.wellRates()[w*np + p] = well_state.wellSolutions()[w] * F[p];
                                }
                                break;
                            }

                            if (well_controls_iget_type(wc, current) == THP) {

                                // Calculate bhp from thp control and well rates
                                double aqua = 0.0;
                                double liquid = 0.0;
                                double vapour = 0.0;

                                const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                                if ((*active_)[ Water ]) {
                                    aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                                }
                                if ((*active_)[ Oil ]) {
                                    liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                                }
                                if ((*active_)[ Gas ]) {
                                    vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                                }

                                const int vfp        = well_controls_iget_vfp(wc, current);
                                const double& thp    = well_controls_iget_target(wc, current);
                                const double& alq    = well_controls_iget_alq(wc, current);

                                //Set *BHP* target by calculating bhp from THP
                                const WellType& well_type = wells().type[w];
                                // pick the density in the top layer
                                const int perf = wells().well_connpos[w];
                                const double rho = well_perforation_densities_[perf];

                                if (well_type == INJECTOR) {
                                    double dp = wellhelpers::computeHydrostaticCorrection(
                                                wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                                rho, gravity_);

                                    well_state.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                                }
                                else if (well_type == PRODUCER) {
                                    double dp = wellhelpers::computeHydrostaticCorrection(
                                                wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                                rho, gravity_);

                                    well_state.bhp()[w] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                                }
                                else {
                                    OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                                }
                            }

                        }
                            break;
                        case SURFACE_RATE: // Both rate controls use bhp as first well variable
                        case RESERVOIR_RATE:
                        {
                            const int sign1 = dwells[w][flowPhaseToEbosCompIdx(0)] > 0 ? 1: -1;
                            const double dx1_limited = sign1 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(0)]),std::abs(xvar_well_old[w])*dBHPLimit);
                            well_state.wellSolutions()[w] = std::max(xvar_well_old[w] - dx1_limited,1e5);
                            well_state.bhp()[w] = well_state.wellSolutions()[w];

                            if (well_controls_iget_type(wc, current) == SURFACE_RATE) {
                                if (wells().type[w]==PRODUCER) {

                                    double F_target = 0.0;
                                    for (int p = 0; p < np; ++p) {
                                        F_target += wells().comp_frac[np*w + p] * F[p];
                                    }
                                    for (int p = 0; p < np; ++p) {
                                        well_state.wellRates()[np*w + p] = F[p] * target_rate /F_target;
                                    }
                                } else {

                                    for (int p = 0; p < np; ++p) {
                                        well_state.wellRates()[w*np + p] = wells().comp_frac[np*w + p] * target_rate;
                                    }
                                }
                            } else { // RESERVOIR_RATE
                                for (int p = 0; p < np; ++p) {
                                    well_state.wellRates()[np*w + p] = F[p] * target_rate;
                                }
                            }
                        }
                            break;
                        }
                    }
                }
            }



            template <class WellState>
            void updateWellControls(WellState& xw)
            {
                if( !localWellsActive() ) return ;

                std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
                // Find, for each well, if any constraints are broken. If so,
                // switch control to first broken constraint.
                const int np = wells().number_of_phases;
                const int nw = wells().number_of_wells;
        #pragma omp parallel for schedule(dynamic)
                for (int w = 0; w < nw; ++w) {
                    WellControls* wc = wells().ctrls[w];
                    // The current control in the well state overrides
                    // the current control set in the Wells struct, which
                    // is instead treated as a default.
                    int current = xw.currentControls()[w];
                    // Loop over all controls except the current one, and also
                    // skip any RESERVOIR_RATE controls, since we cannot
                    // handle those.
                    const int nwc = well_controls_get_num(wc);
                    int ctrl_index = 0;
                    for (; ctrl_index < nwc; ++ctrl_index) {
                        if (ctrl_index == current) {
                            // This is the currently used control, so it is
                            // used as an equation. So this is not used as an
                            // inequality constraint, and therefore skipped.
                            continue;
                        }
                        if (wellhelpers::constraintBroken(
                                xw.bhp(), xw.thp(), xw.wellRates(),
                                w, np, wells().type[w], wc, ctrl_index)) {
                            // ctrl_index will be the index of the broken constraint after the loop.
                            break;
                        }
                    }
                    if (ctrl_index != nwc) {
                        // Constraint number ctrl_index was broken, switch to it.
                        // We disregard terminal_ouput here as with it only messages
                        // for wells on one process will be printed.
                        std::ostringstream ss;
                        ss << "Switching control mode for well " << wells().name[w]
                              << " from " << modestring[well_controls_iget_type(wc, current)]
                              << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                        OpmLog::info(ss.str());
                        xw.currentControls()[w] = ctrl_index;
                        current = xw.currentControls()[w];
                        well_controls_set_current( wc, current);



                        // Updating well state and primary variables if constraint is broken

                        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
                        const double target = well_controls_iget_target(wc, current);
                        const double* distr = well_controls_iget_distr(wc, current);
                        switch (well_controls_iget_type(wc, current)) {
                        case BHP:
                            xw.bhp()[w] = target;
                            break;

                        case THP: {
                            double aqua = 0.0;
                            double liquid = 0.0;
                            double vapour = 0.0;

                            const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                            if ((*active_)[ Water ]) {
                                aqua = xw.wellRates()[w*np + pu.phase_pos[ Water ] ];
                            }
                            if ((*active_)[ Oil ]) {
                                liquid = xw.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                            }
                            if ((*active_)[ Gas ]) {
                                vapour = xw.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                            }

                            const int vfp        = well_controls_iget_vfp(wc, current);
                            const double& thp    = well_controls_iget_target(wc, current);
                            const double& alq    = well_controls_iget_alq(wc, current);

                            //Set *BHP* target by calculating bhp from THP
                            const WellType& well_type = wells().type[w];

                            // pick the density in the top layer
                            const int perf = wells().well_connpos[w];
                            const double rho = well_perforation_densities_[perf];

                            if (well_type == INJECTOR) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                            rho, gravity_);

                                xw.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                            }
                            else if (well_type == PRODUCER) {
                                double dp = wellhelpers::computeHydrostaticCorrection(
                                            wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                            rho, gravity_);

                                xw.bhp()[w] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
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
                            const WellType& well_type = wells().type[w];
                            if (well_type == INJECTOR) {
                                for (int phase = 0; phase < np; ++phase) {
                                    const double& compi = wells().comp_frac[np * w + phase];
                                    //if (compi > 0.0) {
                                        xw.wellRates()[np*w + phase] = target * compi;
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
                                        xw.wellRates()[np*w + phase] = target * distr[phase];
                                    }
                                }
                            } else {
                                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                            }


                            break;
                        }
                        std::vector<double> g = {1,1,0.01};
                        if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                            const double* distr = well_controls_iget_distr(wc, current);
                            for (int phase = 0; phase < np; ++phase) {
                                g[phase] = distr[phase];
                            }
                        }
                        switch (well_controls_iget_type(wc, current)) {
                        case THP:
                        case BHP:
                        {
                            const WellType& well_type = wells().type[w];
                            xw.wellSolutions()[w] = 0.0;
                            if (well_type == INJECTOR) {
                                for (int p = 0; p < np; ++p)  {
                                    xw.wellSolutions()[w] += xw.wellRates()[np*w + p] * wells().comp_frac[np*w + p];
                                }
                            } else {
                                for (int p = 0; p < np; ++p)  {
                                    xw.wellSolutions()[w] += g[p] * xw.wellRates()[np*w + p];
                                }
                            }
                        }
                            break;


                        case RESERVOIR_RATE: // Intentional fall-through
                        case SURFACE_RATE:
                        {
                            xw.wellSolutions()[w] = xw.bhp()[w];
                        }
                            break;
                        }

                        double tot_well_rate = 0.0;
                        for (int p = 0; p < np; ++p)  {
                            tot_well_rate += g[p] * xw.wellRates()[np*w + p];
                        }
                        if(std::abs(tot_well_rate) > 0) {
                            xw.wellSolutions()[nw + w] = g[Water] * xw.wellRates()[np*w + Water] / tot_well_rate;
                            xw.wellSolutions()[2*nw + w] = g[Gas] * xw.wellRates()[np*w + Gas] / tot_well_rate ;
                        } else {
                            xw.wellSolutions()[nw + w] =  wells().comp_frac[np*w + Water];
                            xw.wellSolutions()[2 * nw + w] =  wells().comp_frac[np*w + Gas];
                        }
                    }
                }
            }


            int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
            {
                const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
                return flowToEbos[ phaseIdx ];
            }

            /// upate the dynamic lists related to economic limits
            template<class WellState>
            void
            updateListEconLimited(const Schedule& schedule,
                                  const int current_step,
                                  const Wells* wells_struct,
                                  const WellState& well_state,
                                  DynamicListEconLimited& list_econ_limited) const
            {
                const int nw = wells_struct->number_of_wells;

                for (int w = 0; w < nw; ++w) {
                    // flag to check if the mim oil/gas rate limit is violated
                    bool rate_limit_violated = false;
                    const std::string& well_name = wells_struct->name[w];
                    const Well* well_ecl = schedule.getWell(well_name);
                    const WellEconProductionLimits& econ_production_limits = well_ecl->getEconProductionLimits(current_step);

                    // economic limits only apply for production wells.
                    if (wells_struct->type[w] != PRODUCER) {
                        continue;
                    }

                    // if no limit is effective here, then continue to the next well
                    if ( !econ_production_limits.onAnyEffectiveLimit() ) {
                        continue;
                    }
                    // for the moment, we only handle rate limits, not handling potential limits
                    // the potential limits should not be difficult to add
                    const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
                    if (quantity_limit == WellEcon::POTN) {
                        const std::string msg = std::string("POTN limit for well ") + well_name + std::string(" is not supported for the moment. \n")
                                              + std::string("All the limits will be evaluated based on RATE. ");
                        OpmLog::warning("NOT_SUPPORTING_POTN", msg);
                    }

                    const WellMapType& well_map = well_state.wellMap();
                    const typename WellMapType::const_iterator i_well = well_map.find(well_name);
                    assert(i_well != well_map.end()); // should always be found?
                    const WellMapEntryType& map_entry = i_well->second;
                    const int well_number = map_entry[0];

                    if (econ_production_limits.onAnyRateLimit()) {
                        rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state, well_number);
                    }

                    if (rate_limit_violated) {
                        if (econ_production_limits.endRun()) {
                            const std::string warning_message = std::string("ending run after well closed due to economic limits is not supported yet \n")
                                                              + std::string("the program will keep running after ") + well_name + std::string(" is closed");
                            OpmLog::warning("NOT_SUPPORTING_ENDRUN", warning_message);
                        }

                        if (econ_production_limits.validFollowonWell()) {
                            OpmLog::warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
                        }

                        if (well_ecl->getAutomaticShutIn()) {
                            list_econ_limited.addShutWell(well_name);
                            const std::string msg = std::string("well ") + well_name + std::string(" will be shut in due to economic limit");
                            OpmLog::info(msg);
                        } else {
                            list_econ_limited.addStoppedWell(well_name);
                            const std::string msg = std::string("well ") + well_name + std::string(" will be stopped due to economic limit");
                            OpmLog::info(msg);
                        }
                        // the well is closed, not need to check other limits
                        continue;
                    }

                    // checking for ratio related limits, mostly all kinds of ratio.
                    bool ratio_limits_violated = false;
                    RatioCheckTuple ratio_check_return;

                    if (econ_production_limits.onAnyRatioLimit()) {
                        ratio_check_return = checkRatioEconLimits(econ_production_limits, well_state, map_entry);
                        ratio_limits_violated = std::get<0>(ratio_check_return);
                    }

                    if (ratio_limits_violated) {
                        const bool last_connection = std::get<1>(ratio_check_return);
                        const int worst_offending_connection = std::get<2>(ratio_check_return);

                        const int perf_start = map_entry[1];

                        assert((worst_offending_connection >= 0) && (worst_offending_connection <  map_entry[2]));

                        const int cell_worst_offending_connection = wells_struct->well_cells[perf_start + worst_offending_connection];
                        list_econ_limited.addClosedConnectionsForWell(well_name, cell_worst_offending_connection);
                        const std::string msg = std::string("Connection ") + std::to_string(worst_offending_connection) + std::string(" for well ")
                                              + well_name + std::string(" will be closed due to economic limit");
                        OpmLog::info(msg);

                        if (last_connection) {
                            list_econ_limited.addShutWell(well_name);
                            const std::string msg2 = well_name + std::string(" will be shut due to the last connection closed");
                            OpmLog::info(msg2);
                        }
                    }

                }
            }

            template <class WellState>
            void computeWellConnectionDensitesPressures(const WellState& xw,
                                                        const std::vector<double>& b_perf,
                                                        const std::vector<double>& rsmax_perf,
                                                        const std::vector<double>& rvmax_perf,
                                                        const std::vector<double>& surf_dens_perf,
                                                        const std::vector<double>& depth_perf,
                                                        const double grav) {
                // Compute densities
                well_perforation_densities_ =
                        WellDensitySegmented::computeConnectionDensities(
                                wells(), xw, fluid_->phaseUsage(),
                                b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

                // Compute pressure deltas
                well_perforation_pressure_diffs_ =
                        WellDensitySegmented::computeConnectionPressureDelta(
                                wells(), depth_perf, well_perforation_densities_, grav);



            }

        protected:
            bool wells_active_;
            const Wells*   wells_;
            ModelParameters param_;
            bool terminal_output_;

            const BlackoilPropsAdInterface* fluid_;
            const std::vector<bool>*  active_;
            const VFPProperties* vfp_properties_;
            double gravity_;
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

                    const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                    if ((*active_)[ Water ]) {
                        aqua = getQs(wellIdx, pu.phase_pos[ Water]);
                    }
                    if ((*active_)[ Oil ]) {
                        liquid = getQs(wellIdx, pu.phase_pos[ Oil ]);
                    }
                    if ((*active_)[ Gas ]) {
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
                return wellVariables_[wellIdx];
            }

            EvalWell getQs(const int wellIdx, const int phaseIdx) const {
                EvalWell qs = 0.0;
                const WellControls* wc = wells().ctrls[wellIdx];
                const int np = wells().number_of_phases;
                const double target_rate = well_controls_get_current_target(wc);

                if (wells().type[wellIdx] == INJECTOR) {
                    const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
                    if (comp_frac == 0.0)
                        return qs;

                    if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                        return wellVariables_[wellIdx];
                    }
                    qs.setValue(target_rate);
                    return qs;
                }

                // Producers
                if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP ) {
                    return wellVariables_[wellIdx] * wellVolumeFractionScaled(wellIdx,phaseIdx);
                }
                if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                    const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];

                    if (comp_frac == 1.0) {
                        qs.setValue(target_rate);
                        return qs;
                    }
                    int currentControlIdx = 0;
                    for (int i = 0; i < np; ++i) {
                        currentControlIdx += wells().comp_frac[np*wellIdx + i] * i;
                    }

                    if (wellVolumeFractionScaled(wellIdx,currentControlIdx) == 0) {
                        return qs;
                    }
                    return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / wellVolumeFractionScaled(wellIdx,currentControlIdx));
                }
                // ReservoirRate
                return target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx);
            }

            EvalWell wellVolumeFraction(const int wellIdx, const int phaseIdx) const {
                assert(wells().number_of_phases == 3);
                const int nw = wells().number_of_wells;
                if (phaseIdx == Water) {
                   return wellVariables_[nw + wellIdx];
                }

                if (phaseIdx == Gas) {
                   return wellVariables_[2*nw + wellIdx];
                }

                // Oil
                return 1.0 - wellVariables_[nw + wellIdx] - wellVariables_[2 * nw + wellIdx];
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
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                const int np = well_state.numPhases();

                if (econ_production_limits.onMinOilRate()) {
                    assert((*active_)[Oil]);
                    const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
                    const double min_oil_rate = econ_production_limits.minOilRate();
                    if (std::abs(oil_rate) < min_oil_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinGasRate() ) {
                    assert((*active_)[Gas]);
                    const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
                    const double min_gas_rate = econ_production_limits.minGasRate();
                    if (std::abs(gas_rate) < min_gas_rate) {
                        return true;
                    }
                }

                if (econ_production_limits.onMinLiquidRate() ) {
                    assert((*active_)[Oil]);
                    assert((*active_)[Water]);
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
                const Opm::PhaseUsage& pu = fluid_->phaseUsage();
                const int well_number = map_entry[0];

                assert((*active_)[Oil]);
                assert((*active_)[Water]);

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

        };


} // namespace Opm
#endif
