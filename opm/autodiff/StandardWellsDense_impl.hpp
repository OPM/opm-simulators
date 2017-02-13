


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
            for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {

                const int cell_idx = wells().well_cells[perf];
                const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                std::vector<EvalWell> cq_s(np,0.0);
                const EvalWell bhp = getBhp(w);
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


} // namespace Opm
