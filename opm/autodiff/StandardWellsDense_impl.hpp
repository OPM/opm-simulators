


namespace Opm {

    template<typename FluidSystem, typename BlackoilIndices>
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    StandardWellsDense(const Wells* wells_arg,
                       WellCollection* well_collection,
                       const ModelParameters& param,
                       const bool terminal_output)
       : wells_active_(wells_arg!=nullptr)
       , wells_(wells_arg)
       , well_collection_(well_collection)
       , param_(param)
       , terminal_output_(terminal_output)
       , well_perforation_efficiency_factors_((wells_!=nullptr ? wells_->well_connpos[wells_->number_of_wells] : 0), 1.0)
       , well_perforation_densities_( wells_ ? wells_arg->well_connpos[wells_arg->number_of_wells] : 0)
       , well_perforation_pressure_diffs_( wells_ ? wells_arg->well_connpos[wells_arg->number_of_wells] : 0)
       , wellVariables_( wells_ ? (wells_arg->number_of_wells * wells_arg->number_of_phases) : 0)
       , F0_(wells_ ? (wells_arg->number_of_wells * wells_arg->number_of_phases) : 0 )
    {
        if( wells_ )
        {
            invDuneD_.setBuildMode( Mat::row_wise );
            duneC_.setBuildMode( Mat::row_wise );
            duneB_.setBuildMode( Mat::row_wise );
         }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    init(const PhaseUsage phase_usage_arg,
         const std::vector<bool>& active_arg,
         const VFPProperties*  vfp_properties_arg,
         const double gravity_arg,
         const std::vector<double>& depth_arg,
         const std::vector<double>& pv_arg,
         const RateConverterType* rate_converter)
    {
        if ( ! localWellsActive() ) {
            return;
        }

        phase_usage_ = phase_usage_arg;
        active_ = active_arg;
        vfp_properties_ = vfp_properties_arg;
        gravity_ = gravity_arg;
        cell_depths_ = extractPerfData(depth_arg);
        pv_ = pv_arg;
        rate_converter_ = rate_converter;

        calculateEfficiencyFactors();

        // setup sparsity pattern for the matrices
        //[A B^T    [x    =  [ res
        // C D] x_well]      res_well]

        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const int nc = numCells();

#ifndef NDEBUG
        const auto pu = phase_usage_;
        const int np = pu.num_phases;

        // assumes the gas fractions are stored after water fractions
        // WellVariablePositions needs to be changed for 2p runs
        assert (np == 3 || (np == 2 && !pu.phase_used[Gas]) );
#endif

        // set invDuneD
        invDuneD_.setSize( nw, nw, nw );

        // set duneC
        duneC_.setSize( nw, nc, nperf );

        // set duneB
        duneB_.setSize( nw, nc, nperf );

        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        // make the B^T matrix
        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize( nw );

        // resize temporary class variables
        Cx_.resize( duneC_.N() );
        invDrw_.resize( invDuneD_.N() );
    }






    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    SimulatorReport
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    assemble(Simulator& ebosSimulator,
             const int iterationIdx,
             const double dt,
             WellState& well_state)
    {
        SimulatorReport report;
        if ( ! localWellsActive() ) {
            return report;
        }

        if (param_.compute_well_potentials_) {
            computeWellPotentials(ebosSimulator, well_state);
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
            report = solveWellEq(ebosSimulator, dt, well_state);
        }
        assembleWellEq(ebosSimulator, dt, well_state, false);

        report.converged = true;
        return report;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    assembleWellEq(Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells)
    {
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
            const EvalWell bhp = getBhp(w);
            for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {

                const int cell_idx = wells().well_cells[perf];
                const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                std::vector<EvalWell> cq_s(np,0.0);
                computeWellFlux(w, wells().WI[perf], intQuants, bhp, wellPerforationPressureDiffs()[perf], allow_cf, cq_s);

                for (int p1 = 0; p1 < np; ++p1) {

                    // the cq_s entering mass balance equations need to consider the efficiency factors.
                    const EvalWell cq_s_effective = cq_s[p1] * well_perforation_efficiency_factors_[perf];

                    if (!only_wells) {
                        // subtract sum of phase fluxes in the reservoir equation.
                        // need to consider the efficiency factor
                        ebosResid[cell_idx][flowPhaseToEbosCompIdx(p1)] -= cq_s_effective.value();
                    }

                    // subtract sum of phase fluxes in the well equations.
                    resWell_[w][flowPhaseToEbosCompIdx(p1)] -= cq_s[p1].value();

                    // assemble the jacobians
                    for (int p2 = 0; p2 < np; ++p2) {
                        if (!only_wells) {
                            // also need to consider the efficiency factor when manipulating the jacobians.
                            ebosJac[cell_idx][cell_idx][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s_effective.derivative(p2);
                            duneB_[w][cell_idx][flowToEbosPvIdx(p2)][flowPhaseToEbosCompIdx(p1)] -= cq_s_effective.derivative(p2+blocksize); // intput in transformed matrix
                            duneC_[w][cell_idx][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s_effective.derivative(p2);
                        }
                        invDuneD_[w][w][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] -= cq_s[p1].derivative(p2+blocksize);
                    }

                    // add trivial equation for 2p cases (Only support water + oil)
                    if (np == 2) {
                        assert(!active_[ Gas ]);
                        invDuneD_[w][w][flowPhaseToEbosCompIdx(Gas)][flowToEbosPvIdx(Gas)] = 1.0;
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
                    invDuneD_[w][w][flowPhaseToEbosCompIdx(p1)][flowToEbosPvIdx(p2)] += resWell_loc.derivative(p2+blocksize);
                }
                resWell_[w][flowPhaseToEbosCompIdx(p1)] += resWell_loc.value();
            }
        }

        // do the local inversion of D.
        localInvert( invDuneD_ );
    }





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    allow_cross_flow(const int w, Simulator& ebosSimulator) const
    {
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

            if (drawdown.value() < 0 && wells().type[w] == INJECTOR)  {
                return false;
            }

            if (drawdown.value() > 0 && wells().type[w] == PRODUCER)  {
                return false;
            }
        }
        return true;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    localInvert(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                //std::cout << (*col) << std::endl;
                (*col).invert();
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    print(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                std::cout << row.index() << " " << col.index() << "/n \n"<<(*col) << std::endl;
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    apply( BVector& r) const
    {
        if ( ! localWellsActive() ) {
            return;
        }

        assert( invDrw_.size() == invDuneD_.N() );

        invDuneD_.mv(resWell_,invDrw_);
        duneB_.mmtv(invDrw_, r);
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    apply(const BVector& x, BVector& Ax)
    {
        if ( ! localWellsActive() ) {
            return;
        }

        assert( Cx_.size() == duneC_.N() );

        BVector& invDCx = invDrw_;
        assert( invDCx.size() == invDuneD_.N());

        duneC_.mv(x, Cx_);
        invDuneD_.mv(Cx_, invDCx);
        duneB_.mmtv(invDCx,Ax);
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax)
    {
        if ( ! localWellsActive() ) {
            return;
        }

        if( scaleAddRes_.size() != Ax.size() ) {
            scaleAddRes_.resize( Ax.size() );
        }

        scaleAddRes_ = 0.0;
        apply( x, scaleAddRes_ );
        Ax.axpy( alpha, scaleAddRes_ );
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    recoverVariable(const BVector& x, BVector& xw) const
    {
        if ( ! localWellsActive() ) {
             return;
        }
        BVector resWell = resWell_;
        duneC_.mmv(x, resWell);
        invDuneD_.mv(resWell, xw);
    }





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx };
        return phaseToComp[ phaseIdx ];
    }





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    flowToEbosPvIdx( const int flowPv ) const
    {
        const int flowToEbos[ 3 ] = {
                                     BlackoilIndices::pressureSwitchIdx,
                                     BlackoilIndices::waterSaturationIdx,
                                     BlackoilIndices::compositionSwitchIdx
                                    };
        return flowToEbos[ flowPv ];
    }





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    ebosCompToFlowPhaseIdx( const int compIdx ) const
    {
        const int compToPhase[ 3 ] = { Oil, Water, Gas };
        return compToPhase[ compIdx ];
    }





    template<typename FluidSystem, typename BlackoilIndices>
    std::vector<double>
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    extractPerfData(const std::vector<double>& in) const
    {
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





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    numPhases() const
    {
        return wells().number_of_phases;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    numCells() const
    {
        return pv_.size();
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    resetWellControlFromState(WellState xw)
    {
        const int        nw   = wells_->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells_->ctrls[w];
            well_controls_set_current( wc, xw.currentControls()[w]);
        }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    const Wells&
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }





    template<typename FluidSystem, typename BlackoilIndices>
    const Wells*
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wellsPointer() const
    {
        return wells_;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wellsActive() const
    {
        return wells_active_;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    numWellVars() const
    {
        if ( !localWellsActive() ) {
            return 0;
        }

        // For each well, we have a bhp variable, and one flux per phase.
        const int nw = wells().number_of_wells;
        return (numPhases() + 1) * nw;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    const std::vector<double>&
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wellPerforationDensities() const
    {
         return well_perforation_densities_;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    const std::vector<double>&
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }


    template<typename FluidSystem, typename BlackoilIndices>
    typename StandardWellsDense<FluidSystem, BlackoilIndices>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    extendEval(Eval in) const {
        EvalWell out = 0.0;
        out.setValue(in.value());
        for(int i = 0; i < blocksize;++i) {
            out.setDerivative(i, in.derivative(flowToEbosPvIdx(i)));
        }
        return out;
    }


    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    setWellVariables(const WellState& xw)
    {
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
            for (int w = 0; w < nw; ++w) {
                wellVariables_[w + nw*phaseIdx] = 0.0;
                wellVariables_[w + nw*phaseIdx].setValue(xw.wellSolutions()[w + nw* phaseIdx]);
                wellVariables_[w + nw*phaseIdx].setDerivative(blocksize + phaseIdx, 1.0);
            }
        }
    }

    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    print(EvalWell in) const
    {
        std::cout << in.value() << std::endl;
        for (int i = 0; i < in.size; ++i) {
            std::cout << in.derivative(i) << std::endl;
        }
    }

    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    computeAccumWells()
    {
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
            for (int w = 0; w < nw; ++w) {
                F0_[w + nw * phaseIdx] = wellVolumeFraction(w,phaseIdx).value();
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    template<typename intensiveQuants>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    computeWellFlux(const int& w, const double& Tw,
                    const intensiveQuants& intQuants,
                    const EvalWell& bhp, const double& cdp,
                    const bool& allow_cf, std::vector<EvalWell>& cq_s)  const
    {
        const Opm::PhaseUsage& pu = phase_usage_;
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
            if (!allow_cf && wells().type[w] == INJECTOR) {
                return;
            }

            // compute phase volumetric rates at standard conditions
            std::vector<EvalWell> cq_ps(np, 0.0);
            for (int phase = 0; phase < np; ++phase) {
                const EvalWell cq_p = - Tw * (mob_perfcells_dense[phase] * drawdown);
                cq_ps[phase] = b_perfcells_dense[phase] * cq_p;
            }

            if (active_[Oil] && active_[Gas]) {
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
            if (!allow_cf && wells().type[w] == PRODUCER) {
                return;
            }

            // Using total mobilities
            EvalWell total_mob_dense = mob_perfcells_dense[0];
            for (int phase = 1; phase < np; ++phase) {
                total_mob_dense += mob_perfcells_dense[phase];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - Tw * (total_mob_dense * drawdown);

            // compute volume ratio between connection at standard conditions
            EvalWell volumeRatio = 0.0;
            if (active_[Water]) {
                const int watpos = pu.phase_pos[Water];
                volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
            }

            if (active_[Oil] && active_[Gas]) {
                EvalWell well_temperature = extendEval(fs.temperature(FluidSystem::oilPhaseIdx));
                EvalWell rsSatEval = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);
                EvalWell rvSatEval = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);

                const int oilpos = pu.phase_pos[Oil];
                const int gaspos = pu.phase_pos[Gas];
                EvalWell rvPerf = 0.0;
                if (cmix_s[gaspos] > 0) {
                    rvPerf = cmix_s[oilpos] / cmix_s[gaspos];
                }

                if (rvPerf.value() > rvSatEval.value()) {
                    rvPerf = rvSatEval;
                    //rvPerf.setValue(rvSatEval.value());
                }

                EvalWell rsPerf = 0.0;
                if (cmix_s[oilpos] > 0) {
                    rsPerf = cmix_s[gaspos] / cmix_s[oilpos];
                }

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
                if (active_[Oil]) {
                    const int oilpos = pu.phase_pos[Oil];
                    volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                }
                if (active_[Gas]) {
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





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    SimulatorReport
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    solveWellEq(Simulator& ebosSimulator,
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

            // checking whether the group targets are converged
            if (wellCollection()->groupControlActive()) {
                converged = converged && wellCollection()->groupTargetConverged(well_state.wellRates());
            }

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

        SimulatorReport report;
        report.converged = converged;
        report.total_well_iterations = it;
        return report;
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    printIf(const int c, const double x, const double y, const double eps, const std::string type) const
    {
        if (std::abs(x-y) > eps) {
            std::cout << type << " " << c << ": "<<x << " " << y << std::endl;
        }
    }





    template<typename FluidSystem, typename BlackoilIndices>
    std::vector<double>
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    residual() const
    {
        if( ! wellsActive() )
        {
            return std::vector<double>();
        }

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





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    getWellConvergence(Simulator& ebosSimulator,
                       const int iteration) const
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





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
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





    template<typename FluidSystem, typename BlackoilIndices>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& xw,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf) const
    {
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;
        const PhaseUsage& pu = phase_usage_;
        const int np = phase_usage_.num_phases;
        b_perf.resize(nperf*np);
        surf_dens_perf.resize(nperf*np);

        //rs and rv are only used if both oil and gas is present
        if (pu.phase_used[BlackoilPhases::Vapour] && pu.phase_pos[BlackoilPhases::Liquid]) {
            rsmax_perf.resize(nperf);
            rvmax_perf.resize(nperf);
        }

        // Compute the average pressure in each well block
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {

                const int cell_idx = wells().well_cells[perf];
                const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                const auto& fs = intQuants.fluidState();

                const double p_above = perf == wells().well_connpos[w] ? xw.bhp()[w] : xw.perfPress()[perf - 1];
                const double p_avg = (xw.perfPress()[perf] + p_above)/2;
                const double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value();

                if (pu.phase_used[BlackoilPhases::Aqua]) {
                    b_perf[ pu.phase_pos[BlackoilPhases::Aqua] + perf * pu.num_phases] =
                    FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                }

                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    const int gaspos = pu.phase_pos[BlackoilPhases::Vapour] + perf * pu.num_phases;
                    const int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;

                    if (pu.phase_used[BlackoilPhases::Liquid]) {
                        const int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
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
                    const int oilpos = pu.phase_pos[BlackoilPhases::Liquid] + perf * pu.num_phases;
                    const int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
                    if (pu.phase_used[BlackoilPhases::Vapour]) {
                        rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        const int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;
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





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    updateWellState(const BVector& dwells,
                    WellState& well_state) const
    {
        if( !localWellsActive() ) return;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;

        double dFLimit = dWellFractionMax();
        double dBHPLimit = dbhpMaxRel();
        std::vector<double> xvar_well_old = well_state.wellSolutions();

        for (int w = 0; w < nw; ++w) {

            // update the second and third well variable (The flux fractions)
            std::vector<double> F(np,0.0);
            if (active_[ Water ]) {
                const int sign2 = dwells[w][flowPhaseToEbosCompIdx(WFrac)] > 0 ? 1: -1;
                const double dx2_limited = sign2 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(WFrac)]),dFLimit);
                well_state.wellSolutions()[WFrac*nw + w] = xvar_well_old[WFrac*nw + w] - dx2_limited;
            }

            if (active_[ Gas ]) {
                const int sign3 = dwells[w][flowPhaseToEbosCompIdx(GFrac)] > 0 ? 1: -1;
                const double dx3_limited = sign3 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(GFrac)]),dFLimit);
                well_state.wellSolutions()[GFrac*nw + w] = xvar_well_old[GFrac*nw + w] - dx3_limited;
            }

            assert(active_[ Oil ]);
            F[Oil] = 1.0;
            if (active_[ Water ]) {
                F[Water] = well_state.wellSolutions()[WFrac*nw + w];
                F[Oil] -= F[Water];
            }

            if (active_[ Gas ]) {
                F[Gas] = well_state.wellSolutions()[GFrac*nw + w];
                F[Oil] -= F[Gas];
            }

            if (active_[ Water ]) {
                if (F[Water] < 0.0) {
                    if (active_[ Gas ]) {
                        F[Gas] /= (1.0 - F[Water]);
                    }
                    F[Oil] /= (1.0 - F[Water]);
                    F[Water] = 0.0;
                }
            }
            if (active_[ Gas ]) {
                if (F[Gas] < 0.0) {
                    if (active_[ Water ]) {
                        F[Water] /= (1.0 - F[Gas]);
                    }
                    F[Oil] /= (1.0 - F[Gas]);
                    F[Gas] = 0.0;
                }
            }
            if (F[Oil] < 0.0) {
                if (active_[ Water ]) {
                     F[Water] /= (1.0 - F[Oil]);
                }
                if (active_[ Gas ]) {
                    F[Gas] /= (1.0 - F[Oil]);
                }
                F[Oil] = 0.0;
            }

            if (active_[ Water ]) {
                well_state.wellSolutions()[WFrac*nw + w] = F[Water];
            }
            if (active_[ Gas ]) {
                well_state.wellSolutions()[GFrac*nw + w] = F[Gas];
            }

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
                    well_state.wellSolutions()[nw*XvarWell + w] = xvar_well_old[nw*XvarWell + w] - dwells[w][flowPhaseToEbosCompIdx(XvarWell)];

                    switch (wells().type[w]) {
                    case INJECTOR:
                        for (int p = 0; p < np; ++p) {
                            const double comp_frac = wells().comp_frac[np*w + p];
                            well_state.wellRates()[w*np + p] = comp_frac * well_state.wellSolutions()[nw*XvarWell + w];
                        }
                        break;
                    case PRODUCER:
                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[w*np + p] = well_state.wellSolutions()[nw*XvarWell + w] * F[p];
                        }
                        break;
                    }

                    if (well_controls_iget_type(wc, current) == THP) {

                        // Calculate bhp from thp control and well rates
                        double aqua = 0.0;
                        double liquid = 0.0;
                        double vapour = 0.0;

                        const Opm::PhaseUsage& pu = phase_usage_;

                        if (active_[ Water ]) {
                            aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                        }
                        if (active_[ Oil ]) {
                            liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                        }
                        if (active_[ Gas ]) {
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
                             const double dp = wellhelpers::computeHydrostaticCorrection(
                                               wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                               rho, gravity_);

                             well_state.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                        }
                        else if (well_type == PRODUCER) {
                            const double dp = wellhelpers::computeHydrostaticCorrection(
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
                    const int sign1 = dwells[w][flowPhaseToEbosCompIdx(XvarWell)] > 0 ? 1: -1;
                    const double dx1_limited = sign1 * std::min(std::abs(dwells[w][flowPhaseToEbosCompIdx(XvarWell)]),std::abs(xvar_well_old[nw*XvarWell + w])*dBHPLimit);
                    well_state.wellSolutions()[nw*XvarWell + w] = std::max(xvar_well_old[nw*XvarWell + w] - dx1_limited,1e5);
                    well_state.bhp()[w] = well_state.wellSolutions()[nw*XvarWell + w];

                    if (well_controls_iget_type(wc, current) == SURFACE_RATE) {
                        if (wells().type[w]==PRODUCER) {

                            const double* distr = well_controls_iget_distr(wc, current);

                            double F_target = 0.0;
                            for (int p = 0; p < np; ++p) {
                                F_target += distr[p] * F[p];
                            }
                            for (int p = 0; p < np; ++p) {
                                well_state.wellRates()[np*w + p] = F[p] * target_rate / F_target;
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
            } // end of switch (well_controls_iget_type(wc, current))
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename FluidSystem, typename BlackoilIndices>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    updateWellControls(WellState& xw) const
    {
        if( !localWellsActive() ) return ;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;

        // keeping a copy of the current controls, to see whether control changes later.
        std::vector<int> old_control_index(nw, 0);
        for (int w = 0; w < nw; ++w) {
            old_control_index[w] = xw.currentControls()[w];
        }

        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
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
                xw.currentControls()[w] = ctrl_index;
                current = xw.currentControls()[w];
                well_controls_set_current( wc, current);
            }

            // update whether well is under group control
            if (wellCollection()->groupControlActive()) {
                // get well node in the well collection
                WellNode& well_node = well_collection_->findWellNode(std::string(wells().name[w]));

                // update whehter the well is under group control or individual control
                if (well_node.groupControlIndex() >= 0 && current == well_node.groupControlIndex()) {
                    // under group control
                    well_node.setIndividualControl(false);
                } else {
                    // individual control
                    well_node.setIndividualControl(true);
                }
            }
        }

        // upate the well targets following group controls
        if (wellCollection()->groupControlActive()) {
            applyVREPGroupControl(xw);
            wellCollection()->updateWellTargets(xw.wellRates());
        }

        // the new well control indices after all the related updates,
        std::vector<int> updated_control_index(nw, 0);
        for (int w = 0; w < nw; ++w) {
            updated_control_index[w] = xw.currentControls()[w];
        }

        // checking whether control changed
        wellhelpers::WellSwitchingLogger logger;
        for (int w = 0; w < nw; ++w) {
            if (updated_control_index[w] != old_control_index[w]) {
                WellControls* wc = wells().ctrls[w];
                logger.wellSwitched(wells().name[w],
                                    well_controls_iget_type(wc, old_control_index[w]),
                                    well_controls_iget_type(wc, updated_control_index[w]));
                updateWellStateWithTarget(wc, updated_control_index[w], w, xw);
            }
        }
    }


} // namespace Opm
