


namespace Opm {

    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    init(const PhaseUsage phase_usage_arg,
         const std::vector<bool>& active_arg,
         const VFPProperties*  vfp_properties_arg,
         const double gravity_arg,
         const std::vector<double>& depth_arg,
         const std::vector<double>& pv_arg,
         const RateConverterType* rate_converter,
         long int global_nc)
    {
        // has to be set always for the convergence check!
        global_nc_   = global_nc;

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






    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    SimulatorReport
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    assemble(Simulator& ebosSimulator,
             const int iterationIdx,
             const double dt,
             WellState& well_state)
    {

        if (iterationIdx == 0) {
            prepareTimeStep(ebosSimulator, well_state);
        }


        SimulatorReport report;
        if ( ! localWellsActive() ) {
            return report;
        }

        // resetWellControlFromState(well_state);
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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

                std::vector<EvalWell> mob(np, 0.0);
                getMobility(ebosSimulator, perf, cell_idx, mob);
                computeWellFlux(w, wells().WI[perf], intQuants.fluidState(), mob, bhp, wellPerforationPressureDiffs()[perf], allow_cf, cq_s);

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
                EvalWell resWell_loc = (wellSurfaceVolumeFraction(w, p1) - F0_[w + nw*p1]) * volume / dt;
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


    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    getMobility(const Simulator& ebosSimulator, const int perf, const int cell_idx, std::vector<EvalWell>& mob) const
    {

        const int np = wells().number_of_phases;
        assert (mob.size() == np);
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = wells().sat_table_id[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
            }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            Eval relativePerms[3];
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(relativePerms[ebosPhaseIdx] / intQuants.fluidState().viscosity(ebosPhaseIdx));
            }
        }
    }



    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    localInvert(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                //std::cout << (*col) << std::endl;
                (*col).invert();
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    print(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                std::cout << row.index() << " " << col.index() << "/n \n"<<(*col) << std::endl;
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    apply( BVector& r) const
    {
        if ( ! localWellsActive() ) {
            return;
        }

        assert( invDrw_.size() == invDuneD_.N() );

        invDuneD_.mv(resWell_,invDrw_);
        duneB_.mmtv(invDrw_, r);
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    recoverVariable(const BVector& x, BVector& xw) const
    {
        if ( ! localWellsActive() ) {
             return;
        }
        BVector resWell = resWell_;
        duneC_.mmv(x, resWell);
        invDuneD_.mv(resWell, xw);
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx };
        return phaseToComp[ phaseIdx ];
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    flowToEbosPvIdx( const int flowPv ) const
    {
        const int flowToEbos[ 3 ] = {
                                     BlackoilIndices::pressureSwitchIdx,
                                     BlackoilIndices::waterSaturationIdx,
                                     BlackoilIndices::compositionSwitchIdx
                                    };
        return flowToEbos[ flowPv ];
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
        return flowToEbos[ phaseIdx ];
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    ebosCompToFlowPhaseIdx( const int compIdx ) const
    {
        const int compToPhase[ 3 ] = { Oil, Water, Gas };
        return compToPhase[ compIdx ];
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    std::vector<double>
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    numPhases() const
    {
        return wells().number_of_phases;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    numCells() const
    {
        return pv_.size();
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    resetWellControlFromState(const WellState& xw) const
    {
        const int        nw   = wells_->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells_->ctrls[w];
            well_controls_set_current( wc, xw.currentControls()[w]);
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    const Wells&
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    const Wells*
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellsPointer() const
    {
        return wells_;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellsActive() const
    {
        return wells_active_;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    int
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    numWellVars() const
    {
        if ( !localWellsActive() ) {
            return 0;
        }

        // For each well, we have a bhp variable, and one flux per phase.
        const int nw = wells().number_of_wells;
        return (numPhases() + 1) * nw;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    const std::vector<double>&
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellPerforationDensities() const
    {
         return well_perforation_densities_;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    const std::vector<double>&
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    extendEval(Eval in) const {
        EvalWell out = 0.0;
        out.setValue(in.value());
        for(int i = 0; i < blocksize;++i) {
            out.setDerivative(i, in.derivative(flowToEbosPvIdx(i)));
        }
        return out;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    print(EvalWell in) const
    {
        std::cout << in.value() << std::endl;
        for (int i = 0; i < in.size; ++i) {
            std::cout << in.derivative(i) << std::endl;
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    computeAccumWells()
    {
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
            for (int w = 0; w < nw; ++w) {
                F0_[w + nw * phaseIdx] = wellSurfaceVolumeFraction(w, phaseIdx).value();
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template<typename FluidState>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    computeWellFlux(const int& w, const double& Tw,
                    const FluidState& fs,
                    const std::vector<EvalWell>& mob_perfcells_dense,
                    const EvalWell& bhp, const double& cdp,
                    const bool& allow_cf, std::vector<EvalWell>& cq_s)  const
    {
        const Opm::PhaseUsage& pu = phase_usage_;
        const int np = wells().number_of_phases;
        std::vector<EvalWell> cmix_s(np,0.0);
        for (int phase = 0; phase < np; ++phase) {
            //int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
            cmix_s[phase] = wellSurfaceVolumeFraction(w, phase);
        }

        EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        EvalWell rs = extendEval(fs.Rs());
        EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(np, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
            b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
        }

        // Pressure drawdown (also used to determine direction of flow)
        EvalWell well_pressure = bhp + cdp;
        EvalWell drawdown = pressure - well_pressure;

        // producing perforations
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

                // Incorporate RS/RV factors if both oil and gas active
                const EvalWell d = 1.0 - rv * rs;

                if (d.value() == 0.0) {
                    OPM_THROW(Opm::NumericalProblem, "Zero d value obtained for well " << wells().name[w] << " during flux calcuation"
                                                  << " with rs " << rs << " and rv " << rv);
                }

                const EvalWell tmp_oil = (cmix_s[oilpos] - rv * cmix_s[gaspos]) / d;
                //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                const EvalWell tmp_gas = (cmix_s[gaspos] - rs * cmix_s[oilpos]) / d;
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    SimulatorReport
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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
            // also recover the old well controls
            for (int w = 0; w < nw; ++w) {
                WellControls* wc = wells().ctrls[w];
                well_controls_set_current(wc, well_state.currentControls()[w]);
            }
        }

        SimulatorReport report;
        report.converged = converged;
        report.total_well_iterations = it;
        return report;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    printIf(const int c, const double x, const double y, const double eps, const std::string type) const
    {
        if (std::abs(x-y) > eps) {
            std::cout << type << " " << c << ": "<<x << " " << y << std::endl;
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    std::vector<double>
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    getWellConvergence(Simulator& ebosSimulator,
                       const int iteration) const
    {
        typedef double Scalar;
        typedef std::vector< Scalar > Vector;

        const int np = numPhases();
        const double tol_wells = param_.tolerance_wells_;
        const double maxResidualAllowed = param_.max_residual_allowed_;

        std::vector< Scalar > B_avg( np, Scalar() );
        std::vector< Scalar > maxNormWell(np, Scalar() );

        auto& grid = ebosSimulator.gridManager().grid();
        const auto& gridView = grid.leafGridView();
        ElementContext elemCtx(ebosSimulator);
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updatePrimaryStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            for ( int idx = 0; idx < np; ++idx )
            {
                auto& B  = B_avg[ idx ];
                const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);

                B += 1 / fs.invB(ebosPhaseIdx).value();
            }
        }

        // compute global average
        grid.comm().sum(B_avg.data(), B_avg.size());
        for(auto& bval: B_avg)
        {
            bval/=global_nc_;
        }
        
        auto res = residual();
        const int nw = res.size() / np;

        for ( int idx = 0; idx < np; ++idx )
        {
            for ( int w = 0; w < nw; ++w ) {
                maxNormWell[idx] = std::max(maxNormWell[idx], std::abs(res[nw*idx + w]));
            }
        }

        grid.comm().max(maxNormWell.data(), maxNormWell.size());

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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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
                    if (distr[p] > 0.) { // For injection wells, there only one non-zero distr value
                        F[p] /= distr[p];
                    } else {
                        F[p] = 0.;
                    }
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


        // for the wells having a THP constaint, we should update their thp value
        // If it is under THP control, it will be set to be the target value. Otherwise,
        // the thp value will be calculated based on the bhp value, assuming the bhp value is correctly calculated.
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            const int nwc = well_controls_get_num(wc);
            // Looping over all controls until we find a THP constraint
            int ctrl_index = 0;
            for ( ; ctrl_index < nwc; ++ctrl_index) {
                if (well_controls_iget_type(wc, ctrl_index) == THP) {
                    // the current control
                    const int current = well_state.currentControls()[w];
                    // If under THP control at the moment
                    if (current == ctrl_index) {
                        const double thp_target = well_controls_iget_target(wc, current);
                        well_state.thp()[w] = thp_target;
                    } else { // otherwise we calculate the thp from the bhp value
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

                        const double alq = well_controls_iget_alq(wc, ctrl_index);
                        const int table_id = well_controls_iget_vfp(wc, ctrl_index);

                        const WellType& well_type = wells().type[w];
                        const int perf = wells().well_connpos[w]; //first perforation.
                        if (well_type == INJECTOR) {
                            const double dp = wellhelpers::computeHydrostaticCorrection(
                                              wells(), w, vfp_properties_->getInj()->getTable(table_id)->getDatumDepth(),
                                              wellPerforationDensities()[perf], gravity_);

                            const double bhp = well_state.bhp()[w];
                            well_state.thp()[w] = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
                        } else if (well_type == PRODUCER) {
                            const double dp = wellhelpers::computeHydrostaticCorrection(
                                              wells(), w, vfp_properties_->getProd()->getTable(table_id)->getDatumDepth(),
                                              wellPerforationDensities()[perf], gravity_);

                            const double bhp = well_state.bhp()[w];
                            well_state.thp()[w] = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
                        } else {
                            OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                        }
                    }

                    // the THP control is found, we leave the loop now
                    break;
                }
            }  // end of for loop for seaching THP constraints

            // no THP constraint found
            if (ctrl_index == nwc) { // not finding a THP contstraints
                well_state.thp()[w] = 0.0;
            }
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
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

        // the new well control indices after all the related updates,
        std::vector<int> updated_control_index(nw, 0);
        for (int w = 0; w < nw; ++w) {
            updated_control_index[w] = xw.currentControls()[w];
        }

        // checking whether control changed
        wellhelpers::WellSwitchingLogger logger;
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            if (updated_control_index[w] != old_control_index[w]) {
                logger.wellSwitched(wells().name[w],
                                    well_controls_iget_type(wc, old_control_index[w]),
                                    well_controls_iget_type(wc, updated_control_index[w]));
            }

            if (updated_control_index[w] != old_control_index[w] || well_collection_->groupControlActive()) {
                updateWellStateWithTarget(wc, updated_control_index[w], w, xw);
            }
        }

        // upate the well targets following group controls
        // it will not change the control mode, only update the targets
        if (wellCollection()->groupControlActive()) {
            applyVREPGroupControl(xw);
            wellCollection()->updateWellTargets(xw.wellRates());
            for (int w = 0; w < nw; ++w) {
                const WellControls* wc = wells().ctrls[w];
                updateWellStateWithTarget(wc, updated_control_index[w], w, xw);
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    updateListEconLimited(const Schedule& schedule,
                          const int current_step,
                          const Wells* wells_struct,
                          const WellState& well_state,
                          DynamicListEconLimited& list_econ_limited) const
    {
        // With no wells (on process) wells_struct is a null pointer
        const int nw = (wells_struct)? wells_struct->number_of_wells : 0;

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

        } // for (int w = 0; w < nw; ++w)
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    computeWellConnectionDensitesPressures(const WellState& xw,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& surf_dens_perf,
                                           const std::vector<double>& depth_perf,
                                           const double grav)
    {
        // Compute densities
        well_perforation_densities_ =
                  WellDensitySegmented::computeConnectionDensities(
                          wells(), xw, phase_usage_,
                          b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        // Compute pressure deltas
        well_perforation_pressure_diffs_ =
                  WellDensitySegmented::computeConnectionPressureDelta(
                          wells(), depth_perf, well_perforation_densities_, grav);
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    template <typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials)  const
    {

        // number of wells and phases
        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;

        well_potentials.resize(nw * np, 0.0);

        for (int w = 0; w < nw; ++w) {

            bool is_thp_determined = wellHasTHPConstraints(w);

            if (!is_thp_determined) {

            // bhp needs to be determined for the well potential calculation
            // There can be more than one BHP/THP constraints.
            // TODO: there is an option to ignore the THP limit when calculating well potentials,
            // we are not handling it for the moment, while easy to incorporate

            // the bhp will be used to compute well potentials
            double bhp;

            // type of the well, INJECTOR or PRODUCER
            const WellType& well_type = wells().type[w];
            // initial bhp value, making the value not usable
            switch(well_type) {
                case INJECTOR:
                    bhp = std::numeric_limits<double>::max();
                    break;
                case PRODUCER:
                    bhp = -std::numeric_limits<double>::max();
                    break;
                default:
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << wells().name[w]);
            }

            // the well controls
            const WellControls* well_control = wells().ctrls[w];
            // The number of the well controls/constraints
            const int nwc = well_controls_get_num(well_control);

            for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
                // finding a BHP constraint
                if (well_controls_iget_type(well_control, ctrl_index) == BHP) {
                    // get the bhp constraint value, it should always be postive assummingly
                    const double bhp_target = well_controls_iget_target(well_control, ctrl_index);

                    switch(well_type) {
                        case INJECTOR: // using the lower bhp contraint from Injectors
                            if (bhp_target < bhp) {
                                bhp = bhp_target;
                            }
                            break;
                        case PRODUCER:
                            if (bhp_target > bhp) {
                                bhp = bhp_target;
                            }
                            break;
                        default:
                            OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << wells().name[w]);
                    } // end of switch
                }

                // finding a THP constraint
                if (well_controls_iget_type(well_control, ctrl_index) == THP) {
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

                    const int vfp        = well_controls_iget_vfp(well_control, ctrl_index);
                    const double& thp    = well_controls_iget_target(well_control, ctrl_index);
                    const double& alq    = well_controls_iget_alq(well_control, ctrl_index);

                    // Calculating the BHP value based on THP
                    const int first_perf = wells().well_connpos[w]; //first perforation

                    if (well_type == INJECTOR) {
                        const double dp = wellhelpers::computeHydrostaticCorrection(
                                          wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                          wellPerforationDensities()[first_perf], gravity_);
                        const double bhp_calculated = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                        // apply the strictest of the bhp controlls i.e. smallest bhp for injectors
                        if (bhp_calculated < bhp) {
                            bhp = bhp_calculated;
                        }
                    }
                    else if (well_type == PRODUCER) {
                        const double dp = wellhelpers::computeHydrostaticCorrection(
                                          wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                          wellPerforationDensities()[first_perf], gravity_);
                        const double bhp_calculated = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                        // apply the strictest of the bhp controlls i.e. largest bhp for producers
                        if (bhp_calculated > bhp) {
                            bhp = bhp_calculated;
                        }
                    } else {
                       OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }
                }
            }

            // there should be always some avaible bhp/thp constraints there
            assert(std::abs(bhp) != std::numeric_limits<double>::max());

            // Should we consider crossflow when calculating well potentionals?
            const bool allow_cf = allow_cross_flow(w, ebosSimulator);
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const int cell_index = wells().well_cells[perf];
                const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_index, /*timeIdx=*/ 0));
                std::vector<EvalWell> well_potentials_perf(np, 0.0);
                std::vector<EvalWell> mob(np, 0.0);
                getMobility(ebosSimulator, perf, cell_index, mob);
                computeWellFlux(w, wells().WI[perf], intQuants.fluidState(), mob, bhp, wellPerforationPressureDiffs()[perf], allow_cf, well_potentials_perf);
                for(int p = 0; p < np; ++p) {
                    well_potentials[w * np + p] += std::abs(well_potentials_perf[p].value());
                }
            }
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext>
    template<typename Simulator>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    prepareTimeStep(const Simulator& ebos_simulator,
                    WellState& well_state)
    {
        // after restarting, the well_controls can be modified while
        // the well_state still uses the old control index
        // we need to synchronize these two.
        const int nw = wells().number_of_wells;
        for (int w = 0; w < nw; ++w) {
            const int ctrl_index = well_state.currentControls()[w];
            WellControls* wc = wells().ctrls[w];
            const int ctrl_index_2 = well_controls_get_current(wc);
            if (ctrl_index_2 != ctrl_index) {
                well_controls_set_current(wc, ctrl_index);
            }

            // update whether well is under group control
            if (wellCollection()->groupControlActive()) {
                // get well node in the well collection
                WellNode& well_node = well_collection_->findWellNode(std::string(wells().name[w]));

                // update whehter the well is under group control or individual control
                if (well_node.groupControlIndex() >= 0 && ctrl_index == well_node.groupControlIndex()) {
                    // under group control
                    well_node.setIndividualControl(false);
                } else {
                    // individual control
                    well_node.setIndividualControl(true);
                }
            }
        }

        if (well_collection_->groupControlActive()) {
            // calculate the well potentials
            // two functions will probably be merged in the final version
            // and also the well potentials related parts in well state.
            if (param_.compute_well_potentials_) {

                // the following part should be made a function
                // const int nw = wells().number_of_wells;
                for (int w = 0; w < nw; ++w) {
                    WellControls* wc = wells().ctrls[w];
                    const int control = well_controls_get_current(wc);
                    well_state.currentControls()[w] = control;
                    // TODO: when we under defaulted BHP value here, it is not
                    // wise to update the WellState with this target.
                    // It should only be the case with `GRUP` while we have not
                    // applied group control.
                    // updateWellStateWithTarget(wc, control, w, well_state);
                }

                setWellVariables(well_state);
                computeWellConnectionPressures(ebos_simulator, well_state);

                // To store well potentials for each well
                std::vector<double> well_potentials;

                computeWellPotentials(ebos_simulator, well_state, well_potentials);

                // update/setup guide rates for each well based on the well_potentials
                well_collection_->setGuideRatesWithPotentials(wellsPointer(), phase_usage_, well_potentials);

                // handling the situation that wells does not have a valid control
                // it happens the well specified with GRUP and restarting due to non-convergencing
                // putting the well under group control for this situation
                if (wellCollection()->groupControlActive()) {
                    for (int w = 0; w < nw; ++w) {
                        WellControls* wc = wells().ctrls[w];
                        const int ctrl_index = well_controls_get_current(wc);
                        WellNode& well_node = well_collection_->findWellNode(std::string(wells().name[w]));

                        const int group_control_index = well_node.groupControlIndex();
                        if (group_control_index >= 0 && ctrl_index < 0) {
                            well_controls_set_current(wc, group_control_index);
                            well_state.currentControls()[w] = group_control_index;
                            well_node.setIndividualControl(false);
                        }
                    }
                }

            }
            applyVREPGroupControl(well_state);

            if (!wellCollection()->groupControlApplied()) {
                wellCollection()->applyGroupControls();
            } else {
                wellCollection()->updateWellTargets(well_state.wellRates());
            }
        }

        // since the controls are all updated, we should update well_state accordingly
        // const int nw = wells().number_of_wells;
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells().ctrls[w];
            const int control = well_controls_get_current(wc);
            well_state.currentControls()[w] = control;
            updateWellStateWithTarget(wc, control, w, well_state);
        }
    }






    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext>
    WellCollection*
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellCollection() const
    {
        return well_collection_;
    }




    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    const std::vector<double>&
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellPerfEfficiencyFactors() const
    {
        return well_perforation_efficiency_factors_;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    calculateEfficiencyFactors()
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    computeWellVoidageRates(const WellState& well_state,
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    applyVREPGroupControl(WellState& well_state) const
    {
        if ( wellCollection()->havingVREPGroups() ) {
            std::vector<double> well_voidage_rates;
            std::vector<double> voidage_conversion_coeffs;
            computeWellVoidageRates(well_state, well_voidage_rates, voidage_conversion_coeffs);
            wellCollection()->applyVREPGroupControls(well_voidage_rates, voidage_conversion_coeffs);

            // for the wells under group control, update the control index for the well_state and well_controls
            for (const WellNode* well_node : wellCollection()->getLeafNodes()) {
                if (well_node->isInjector() && !well_node->individualControl()) {
                    const int well_index = well_node->selfIndex();
                    well_state.currentControls()[well_index] = well_node->groupControlIndex();

                    WellControls* wc = wells().ctrls[well_index];
                    well_controls_set_current(wc, well_node->groupControlIndex());
                }
            }
        }
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    getBhp(const int wellIdx) const {
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    getQs(const int wellIdx, const int phaseIdx) const
    {
        EvalWell qs = 0.0;
        const WellControls* wc = wells().ctrls[wellIdx];
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const double target_rate = well_controls_get_current_target(wc);

        // TODO: the formulation for the injectors decides it only work with single phase
        // surface rate injection control. Improvement will be required.
        if (wells().type[wellIdx] == INJECTOR) {
            const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
            if (comp_frac == 0.0) {
                return qs;
            }

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

                // looking for the phase under control
                int phase_under_control = -1;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        phase_under_control = phase;
                        break;
                    }
                }

                assert(phase_under_control >= 0);

                if (phaseIdx == phase_under_control) {
                    qs.setValue(target_rate);
                    return qs;
                }

                // TODO: not sure why the single phase under control will have near zero fraction
                const double eps = 1e-6;
                if (wellVolumeFractionScaled(wellIdx, phase_under_control) < eps) {
                    return qs;
                }
                return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / wellVolumeFractionScaled(wellIdx, phase_under_control));
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

            // TODO: three phase surface rate control is not tested yet
            if (num_phases_under_rate_control == 3) {
                return target_rate * wellSurfaceVolumeFraction(wellIdx, phaseIdx);
            }
        } else if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            // ReservoirRate
            return target_rate * wellVolumeFractionScaled(wellIdx, phaseIdx);
        } else {
            OPM_THROW(std::logic_error, "Unknown control type for well " << wells().name[wellIdx]);
        }

        // avoid warning of condition reaches end of non-void function
        return qs;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellVolumeFraction(const int wellIdx, const int phaseIdx) const
    {
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellVolumeFractionScaled(const int wellIdx, const int phaseIdx) const
    {
        const WellControls* wc = wells().ctrls[wellIdx];
        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            const double* distr = well_controls_get_current_distr(wc);
            if (distr[phaseIdx] > 0.) {
                return wellVolumeFraction(wellIdx, phaseIdx) / distr[phaseIdx];
            } else {
                // TODO: not sure why return EvalWell(0.) causing problem here
                // Probably due to the wrong Jacobians.
                return wellVolumeFraction(wellIdx, phaseIdx);
            }
        }
        std::vector<double> g = {1,1,0.01};
        return (wellVolumeFraction(wellIdx, phaseIdx) / g[phaseIdx]);
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::EvalWell
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    wellSurfaceVolumeFraction(const int well_index, const int phase) const
    {
        EvalWell sum_volume_fraction_scaled = 0.;
        const int np = wells().number_of_phases;
        for (int p = 0; p < np; ++p) {
            sum_volume_fraction_scaled += wellVolumeFractionScaled(well_index, p);
        }

        assert(sum_volume_fraction_scaled.value() != 0.);

        return wellVolumeFractionScaled(well_index, phase) / sum_volume_fraction_scaled;
    }





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::RatioCheckTuple
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    typename StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::RatioCheckTuple
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
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





    template<typename FluidSystem, typename BlackoilIndices, typename ElementContext, typename MaterialLaw>
    void
    StandardWellsDense<FluidSystem, BlackoilIndices, ElementContext, MaterialLaw>::
    updateWellStateWithTarget(const WellControls* wc,
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
            // TODO: similar to the way below to handle THP
            // we should not something related to thp here when there is thp constraint
            break;

        case THP: {
            xw.thp()[well_index] = target;

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
                const double dp = wellhelpers::computeHydrostaticCorrection(
                                  wells(), well_index, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                  rho, gravity_);

                xw.bhp()[well_index] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
            }
            else if (well_type == PRODUCER) {
                const double dp = wellhelpers::computeHydrostaticCorrection(
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
            // break;

        case SURFACE_RATE:
            // checking the number of the phases under control
            int numPhasesWithTargetsUnderThisControl = 0;
            for (int phase = 0; phase < np; ++phase) {
                if (distr[phase] > 0.0) {
                    numPhasesWithTargetsUnderThisControl += 1;
                }
            }

            assert(numPhasesWithTargetsUnderThisControl > 0);

            const WellType& well_type = wells().type[well_index];
            if (well_type == INJECTOR) {
                // assign target value as initial guess for injectors
                // only handles single phase control at the moment
                assert(numPhasesWithTargetsUnderThisControl == 1);

                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.) {
                        xw.wellRates()[np*well_index + phase] = target / distr[phase];
                    } else {
                        xw.wellRates()[np * well_index + phase] = 0.;
                    }
                }
            } else if (well_type == PRODUCER) {

                // update the rates of phases under control based on the target,
                // and also update rates of phases not under control to keep the rate ratio,
                // assuming the mobility ratio does not change for the production wells
                double orignal_rates_under_phase_control = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        orignal_rates_under_phase_control += xw.wellRates()[np * well_index + phase] * distr[phase];
                    }
                }

                if (orignal_rates_under_phase_control != 0.0 ) {
                    double scaling_factor = target / orignal_rates_under_phase_control;

                    for (int phase = 0; phase < np; ++phase) {
                        xw.wellRates()[np * well_index + phase] *= scaling_factor;
                    }
                } else { // scaling factor is not well defied when orignal_rates_under_phase_control is zero
                    // separating targets equally between phases under control
                    const double target_rate_devided = target / numPhasesWithTargetsUnderThisControl;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            xw.wellRates()[np * well_index + phase] = target_rate_devided / distr[phase];
                        } else {
                            // this only happens for SURFACE_RATE control
                            xw.wellRates()[np * well_index + phase] = target_rate_devided;
                        }
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
            const WellType& well_type = wells().type[well_index];
            if (well_type == INJECTOR) {
                // only single phase injection handled
                if (active_[Water]) {
                    if (distr[Water] > 0.0) {
                        xw.wellSolutions()[WFrac * nw + well_index] = 1.0;
                    } else {
                        xw.wellSolutions()[WFrac * nw + well_index] = 0.0;
                    }
                }

                if (active_[Gas]) {
                    if (distr[Gas] > 0.0) {
                        xw.wellSolutions()[GFrac * nw + well_index] = 1.0;
                    } else {
                        xw.wellSolutions()[GFrac * nw + well_index] = 0.0;
                    }
                }

                // TODO: it is possible to leave injector as a oil well,
                // when F_w and F_g both equals to zero, not sure under what kind of circumstance
                // this will happen.
            } else if (well_type == PRODUCER) { // producers
                if (active_[Water]) {
                    xw.wellSolutions()[WFrac * nw + well_index] = 1.0 / np;
                }
                if (active_[Gas]) {
                    xw.wellSolutions()[GFrac * nw + well_index] = 1.0 / np;
                }
            } else {
                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
            }
        }

    }





    template<typename FluidSystem, typename BlackoilIndices>
    bool
    StandardWellsDense<FluidSystem, BlackoilIndices>::
    wellHasTHPConstraints(const int well_index) const
    {
        const WellType& well_type = wells().type[well_index];
        const WellControls* well_control = wells().ctrls[well_index];
        const int nwc = well_controls_get_num(well_control);
        for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
            if (well_controls_iget_type(well_control, ctrl_index) == THP) {
                return true;
            }
        }
        return false;
    }

} // namespace Opm
