


namespace Opm {

    template<typename TypeTag>
    StandardWellsDense<TypeTag>::
    StandardWellsDense(const Wells* wells_arg,
                       WellCollection* well_collection,
                       const std::vector< const Well* >& wells_ecl,
                       const ModelParameters& param,
                       const RateConverterType& rate_converter,
                       const bool terminal_output,
                       const int current_timeIdx)
       : wells_active_(wells_arg!=nullptr)
       , wells_(wells_arg)
       , wells_ecl_(wells_ecl)
       , number_of_wells_(wells_arg ? (wells_arg->number_of_wells) : 0)
       , well_collection_(well_collection)
       , param_(param)
       , terminal_output_(terminal_output)
       , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
       , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
       , current_timeIdx_(current_timeIdx)
       , rate_converter_(rate_converter)
       , well_perforation_efficiency_factors_((wells_!=nullptr ? wells_->well_connpos[wells_->number_of_wells] : 0), 1.0)
       , well_perforation_densities_( wells_ ? wells_arg->well_connpos[wells_arg->number_of_wells] : 0)
       , well_perforation_pressure_diffs_( wells_ ? wells_arg->well_connpos[wells_arg->number_of_wells] : 0)
       , wellVariables_( wells_ ? (wells_arg->number_of_wells * numWellEq) : 0)
       , F0_(wells_ ? (wells_arg->number_of_wells * numWellEq) : 0 )
    {
       createWellContainer(wells_arg);
        if( wells_ )
        {
            invDuneD_.setBuildMode( Mat::row_wise );
            duneB_.setBuildMode( Mat::row_wise );
            duneC_.setBuildMode( Mat::row_wise );
         }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    init(const PhaseUsage phase_usage_arg,
         const std::vector<bool>& active_arg,
         const double gravity_arg,
         const std::vector<double>& depth_arg,
         const std::vector<double>& pv_arg,
         long int global_nc,
         const Grid& grid)
    {
        // has to be set always for the convergence check!
        global_nc_   = global_nc;

        if ( ! localWellsActive() ) {
            return;
        }

        phase_usage_ = phase_usage_arg;
        active_ = active_arg;
        gravity_ = gravity_arg;
        cell_depths_ = extractPerfData(depth_arg);
        pv_ = pv_arg;

        calculateEfficiencyFactors();

        // setup sparsity pattern for the matrices
        //[A C^T    [x    =  [ res
        // B D] x_well]      res_well]

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

        // set duneB
        duneB_.setSize( nw, nc, nperf );

        // set duneC
        duneC_.setSize( nw, nc, nperf );

        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        // make the C^T matrix
        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize( nw );

        // resize temporary class variables
        Bx_.resize( duneB_.N() );
        invDrw_.resize( invDuneD_.N() );

        if (has_polymer_)
        {
            if (PolymerModule::hasPlyshlog()) {
                computeRepRadiusPerfLength(grid);
            }
        }

        // do the initialization work

        // do the initialization for all the wells
        // TODO: to see whether we can postpone of the intialization of the well containers to
        // optimize the usage of the following several member variables
        for (auto& well : well_container_) {
            well->init(&phase_usage_, &active_, vfp_properties_, gravity_, nc);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setVFPProperties(const VFPProperties*  vfp_properties_arg)
    {
        vfp_properties_ = vfp_properties_arg;
    }






    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    createWellContainer(const Wells* wells_arg)
    {
        well_container_.clear();
        // There might be no wells in the process
        if (localWellsActive()) {
            const int nw = number_of_wells_;

            well_container_.reserve(nw);

            // With the following way, it will have the same order with wells struct
            // Hopefully, it can generate the same residual history with master branch
            for (int w = 0; w < nw; ++w) {
                const std::string well_name = std::string(wells_arg->name[w]);

                // finding the location of the well in wells_ecl
                const int nw_wells_ecl = wells_ecl_.size();
                int index_well = 0;
                for (; index_well < nw_wells_ecl; ++index_well) {
                    if (well_name == wells_ecl_[index_well]->name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }

                const Well* well_ecl = wells_ecl_[index_well];
                if (well_ecl->getStatus(current_timeIdx_) == WellCommon::SHUT) {
                    continue;
                }

                if (well_ecl->isMultiSegment(current_timeIdx_)) {
                    OPM_THROW(Opm::NumericalProblem, "Not handling Multisegment Wells for now");
                }

                // Basically, we are handling all the wells as StandardWell for the moment
                // TODO: to be changed when we begin introducing MultisegmentWell
                well_container_.push_back(std::make_shared<StandardWell<TypeTag> >(well_ecl, current_timeIdx_, wells_arg) );
            }
        }
    }





    template<typename TypeTag>
    SimulatorReport
    StandardWellsDense<TypeTag>::
    assemble(Simulator& ebosSimulator,
             const int iterationIdx,
             const double dt,
             WellState& well_state)
    {

        if (iterationIdx == 0) {
            prepareTimeStep(ebosSimulator, well_state);
        }

        SimulatorReport report;
        if ( ! wellsActive() ) {
            return report;
        }

        updateWellControls(well_state);
        updateGroupControls(well_state);
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    assembleWellEq(Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells)
    {
        for (int w = 0; w < number_of_wells_; ++w) {
            well_container_[w]->assembleWellEq(ebosSimulator, dt, well_state, only_wells);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag >::
    getMobility(const Simulator& ebosSimulator, const int w, const int perf, const int cell_idx, std::vector<EvalWell>& mob) const
    {

        const int np = wells().number_of_phases;
        assert (int(mob.size()) == numComponents());
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
            if (has_solvent_) {
                mob[solventSaturationIdx] = extendEval(intQuants.solventMobility());
            }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            Eval relativePerms[3] = { 0.0, 0.0, 0.0 };
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(relativePerms[ebosPhaseIdx] / intQuants.fluidState().viscosity(ebosPhaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if (has_solvent_) {
                OPM_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent");
            }
        }

        // modify the water mobility if polymer is present
        if (has_polymer_) {
            // assume fully mixture for wells.
            EvalWell polymerConcentration = extendEval(intQuants.polymerConcentration());

            if (wells().type[w] == INJECTOR) {
                const auto& viscosityMultiplier = PolymerModule::plyviscViscosityMultiplierTable(intQuants.pvtRegionIndex());
                mob[ Water ] /= (extendEval(intQuants.waterViscosityCorrection()) * viscosityMultiplier.eval(polymerConcentration, /*extrapolate=*/true) );
            }

            if (PolymerModule::hasPlyshlog()) {
                // compute the well water velocity with out shear effects.
                const int numComp = numComponents();
                const bool allow_cf = well_container_[w]->crossFlowAllowed(ebosSimulator);
                const EvalWell& bhp = getBhp(w);
                std::vector<EvalWell> cq_s(numComp,0.0);
                computeWellFlux(w, wells().WI[perf], intQuants, mob, bhp, wellPerforationPressureDiffs()[perf], allow_cf, cq_s);
                double area = 2 * M_PI * wells_rep_radius_[perf] * wells_perf_length_[perf];
                const auto& materialLawManager = ebosSimulator.problem().materialLawManager();
                const auto& scaledDrainageInfo =
                        materialLawManager->oilWaterScaledEpsInfoDrainage(cell_idx);
                const Scalar& Swcr = scaledDrainageInfo.Swcr;
                const EvalWell poro = extendEval(intQuants.porosity());
                const EvalWell Sw = extendEval(intQuants.fluidState().saturation(flowPhaseToEbosPhaseIdx(Water)));
                // guard against zero porosity and no water
                const EvalWell denom = Opm::max( (area * poro * (Sw - Swcr)), 1e-12);
                EvalWell waterVelocity = cq_s[ Water ] / denom * extendEval(intQuants.fluidState().invB(flowPhaseToEbosPhaseIdx(Water)));

                if (PolymerModule::hasShrate()) {
                    // TODO Use the same conversion as for the reservoar equations.
                    // Need the "permeability" of the well?
                    // For now use the same formula as in legacy.
                    waterVelocity *= PolymerModule::shrate( intQuants.pvtRegionIndex() ) / wells_bore_diameter_[perf];
                }
                EvalWell polymerConcentration = extendEval(intQuants.polymerConcentration());
                EvalWell shearFactor = PolymerModule::computeShearFactor(polymerConcentration,
                                                                         intQuants.pvtRegionIndex(),
                                                                         waterVelocity);

                // modify the mobility with the shear factor and recompute the well fluxes.
                mob[ Water ] /= shearFactor;
            }
        }

    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    localInvert(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                //std::cout << (*col) << std::endl;
                (*col).invert();
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    print(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                std::cout << row.index() << " " << col.index() << "/n \n"<<(*col) << std::endl;
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    apply( BVector& r) const
    {
        if ( ! localWellsActive() ) {
            return;
        }

        assert( invDrw_.size() == invDuneD_.N() );

        invDuneD_.mv(resWell_,invDrw_);
        duneC_.mmtv(invDrw_, r);
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        if ( ! localWellsActive() ) {
            return;
        }

        assert( Bx_.size() == duneB_.N() );

        BVector& invDBx = invDrw_;
        assert( invDBx.size() == invDuneD_.N());

        duneB_.mv(x, Bx_);
        invDuneD_.mv(Bx_, invDBx);
        duneC_.mmtv(invDBx,Ax);
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    recoverVariable(const BVector& x, BVector& xw) const
    {
        if ( ! localWellsActive() ) {
             return;
        }
        BVector resWell = resWell_;
        duneB_.mmv(x, resWell);
        invDuneD_.mv(resWell, xw);
    }




    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    flowToEbosPvIdx( const int flowPv ) const
    {
        const int flowToEbos[ 3 ] = {
            BlackoilIndices::pressureSwitchIdx,
            BlackoilIndices::waterSaturationIdx,
            BlackoilIndices::compositionSwitchIdx
        };

        if (flowPv > 2 )
            return flowPv;

        return flowToEbos[ flowPv ];
    }

    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx};
        if (phaseIdx > 2 )
            return phaseIdx;
        return phaseToComp[ phaseIdx ];
    }





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        assert(phaseIdx < 3);
        const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
        return flowToEbos[ phaseIdx ];
    }

    template<typename TypeTag>
    std::vector<double>
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    numPhases() const
    {
        return wells().number_of_phases;
    }





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    numCells() const
    {
        return pv_.size();
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    resetWellControlFromState(const WellState& xw) const
    {
        const int        nw   = wells_->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells_->ctrls[w];
            well_controls_set_current( wc, xw.currentControls()[w]);
        }
    }





    template<typename TypeTag>
    const Wells&
    StandardWellsDense<TypeTag>::
    wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }





    template<typename TypeTag>
    const Wells*
    StandardWellsDense<TypeTag>::
    wellsPointer() const
    {
        return wells_;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    wellsActive() const
    {
        return wells_active_;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    numWellVars() const
    {
        if ( !localWellsActive() ) {
            return 0;
        }

        const int nw = wells().number_of_wells;
        return numWellEq * nw;
    }





    template<typename TypeTag>
    const std::vector<double>&
    StandardWellsDense<TypeTag>::
    wellPerforationDensities() const
    {
         return well_perforation_densities_;
    }





    template<typename TypeTag>
    const std::vector<double>&
    StandardWellsDense<TypeTag>::
    wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    extendEval(const Eval& in) const {
        EvalWell out = 0.0;
        out.setValue(in.value());
        for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
            out.setDerivative(eqIdx, in.derivative(flowToEbosPvIdx(eqIdx)));
        }
        return out;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setWellVariables(const WellState& xw)
    {
        for (auto& well_interface : well_container_) {
            well_interface->setWellVariables(xw);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    print(const EvalWell& in) const
    {
        std::cout << in.value() << std::endl;
        for (int i = 0; i < in.size; ++i) {
            std::cout << in.derivative(i) << std::endl;
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeAccumWells()
    {
        for (auto& well : well_container_) {
            well->computeAccumWell();
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellFlux(const int& w, const double& Tw,
                    const IntensiveQuantities& intQuants,
                    const std::vector<EvalWell>& mob_perfcells_dense,
                    const EvalWell& bhp, const double& cdp,
                    const bool& allow_cf, std::vector<EvalWell>& cq_s)  const
    {
        const Opm::PhaseUsage& pu = phase_usage_;
        const int np = wells().number_of_phases;
        const int numComp = numComponents();
        std::vector<EvalWell> cmix_s(numComp,0.0);
        for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
            cmix_s[componentIdx] = wellSurfaceVolumeFraction(w, componentIdx);
        }
        auto& fs = intQuants.fluidState();

        EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        EvalWell rs = extendEval(fs.Rs());
        EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(numComp, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
            b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
        }
        if (has_solvent_) {
            b_perfcells_dense[solventSaturationIdx] = extendEval(intQuants.solventInverseFormationVolumeFactor());
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

            // compute component volumetric rates at standard conditions
            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                const EvalWell cq_p = - Tw * (mob_perfcells_dense[componentIdx] * drawdown);
                cq_s[componentIdx] = b_perfcells_dense[componentIdx] * cq_p;
            }

            if (active_[Oil] && active_[Gas]) {
                const int oilpos = pu.phase_pos[Oil];
                const int gaspos = pu.phase_pos[Gas];
                const EvalWell cq_sOil = cq_s[oilpos];
                const EvalWell cq_sGas = cq_s[gaspos];
                cq_s[gaspos] += rs * cq_sOil;
                cq_s[oilpos] += rv * cq_sGas;
            }

        } else {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && wells().type[w] == PRODUCER) {
                return;
            }

            // Using total mobilities
            EvalWell total_mob_dense = mob_perfcells_dense[0];
            for (int componentIdx = 1; componentIdx < numComp; ++componentIdx) {
                total_mob_dense += mob_perfcells_dense[componentIdx];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - Tw * (total_mob_dense * drawdown);

            // compute volume ratio between connection at standard conditions
            EvalWell volumeRatio = 0.0;
            if (active_[Water]) {
                const int watpos = pu.phase_pos[Water];
                volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
            }

            if (has_solvent_) {
                volumeRatio += cmix_s[solventSaturationIdx] / b_perfcells_dense[solventSaturationIdx];
            }

            if (active_[Oil] && active_[Gas]) {
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
            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is; // * b_perfcells_dense[phase];
            }
        }
    }





    template<typename TypeTag>
    SimulatorReport
    StandardWellsDense<TypeTag>::
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
                for (auto& well : well_container_) {
                    well->wellEqIteration(ebosSimulator, param_, well_state);
                }
            }
            // updateWellControls uses communication
            // Therefore the following is executed if there
            // are active wells anywhere in the global domain.
            if( wellsActive() )
            {
                updateWellControls(well_state);
                updateGroupControls(well_state);
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    printIf(const int c, const double x, const double y, const double eps, const std::string type) const
    {
        if (std::abs(x-y) > eps) {
            std::cout << type << " " << c << ": "<<x << " " << y << std::endl;
        }
    }





    template<typename TypeTag>
    std::vector<double>
    StandardWellsDense<TypeTag>::
    residual() const
    {
        if( ! wellsActive() )
        {
            return std::vector<double>();
        }

        const int nw = wells().number_of_wells;
        const int numComp = numComponents();
        std::vector<double> res(numEq*nw, 0.0);
        for( int compIdx = 0; compIdx < numComp; ++compIdx) {

            for (int wellIdx = 0; wellIdx < nw; ++wellIdx) {
                int idx = wellIdx + nw*compIdx;
                res[idx] = resWell_[ wellIdx ][ compIdx ];
            }
        }
        return res;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    getWellConvergence(Simulator& ebosSimulator,
                       const int iteration) const
    {
        typedef double Scalar;
        typedef std::vector< Scalar > Vector;

        const int numComp = numComponents();

        std::vector< Scalar > B_avg( numComp, Scalar() );
        computeAverageFormationFactor(ebosSimulator, B_avg);

        bool converged_well = true;

        // TODO: to check the strategy here
        // currently, if there is any well not converged, we consider the well eqautions do not get converged
        for (const auto& well : well_container_) {
            if ( !well->getWellConvergence(ebosSimulator, B_avg, param_) ) {
                converged_well = false;
                // break; // TODO: no need to check other wells?
            }
        }

        // TODO: to think about the output here.
        /* if ( terminal_output_ )
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
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                ss << std::setw(11) << well_flux_residual[compIdx];
            }
            ss.precision(oprec);
            ss.flags(oflags);
            OpmLog::note(ss.str());
        } */
        return converged_well;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellConnectionPressures(const Simulator& ebosSimulator,
                                   const WellState& xw)
    {
         if( ! localWellsActive() ) return ;

         for (auto& well : well_container_) {
             well->computeWellConnectionPressures(ebosSimulator, xw);
         }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateWellState(const BVector& dwells,
                    WellState& well_state) const
    {
        // TODO: the interface of the function should change.
        // the current plan is to make different wells have different matrix
        // and residual system.
        if( !localWellsActive() ) return;
        /* for (auto& well : well_container_) {
            well->updateWellState(dwells, param_, well_state);
        } */
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateWellControls(WellState& xw) const
    {
        // Even if there no wells active locally, we cannot
        // return as the Destructor of the WellSwitchingLogger
        // uses global communication. For no well active globally
        // we simply return.
        if( !wellsActive() ) return ;

        for (const auto& well : well_container_) {
            well->updateWellControl(xw);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials) const
    {

        // number of wells and phases
        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;

        well_potentials.resize(nw * np, 0.0);

        for (int w = 0; w < nw; ++w) {

            // get the bhp value based on the bhp constraints
            const double bhp = mostStrictBhpFromBhpLimits(w);

            // does the well have a THP related constraint?
            const bool has_thp_control = wellHasTHPConstraints(w);

            std::vector<double> potentials(np);

            if ( !has_thp_control ) {

                assert(std::abs(bhp) != std::numeric_limits<double>::max());

                computeWellRatesWithBhp(ebosSimulator, bhp, w, potentials);

            } else { // the well has a THP related constraint
                // checking whether a well is newly added, it only happens at the beginning of the report step
                if ( !well_state.isNewWell(w) ) {
                    for (int p = 0; p < np; ++p) {
                        // This is dangerous for new added well
                        // since we are not handling the initialization correctly for now
                        potentials[p] = well_state.wellRates()[w * np + p];
                    }
                } else {
                    // We need to generate a reasonable rates to start the iteration process
                    computeWellRatesWithBhp(ebosSimulator, bhp, w, potentials);
                    for (double& value : potentials) {
                        // make the value a little safer in case the BHP limits are default ones
                        // TODO: a better way should be a better rescaling based on the investigation of the VFP table.
                        const double rate_safety_scaling_factor = 0.00001;
                        value *= rate_safety_scaling_factor;
                    }
                }

                potentials = computeWellPotentialWithTHP(ebosSimulator, w, bhp, potentials);
            }

            // putting the sucessfully calculated potentials to the well_potentials
            for (int p = 0; p < np; ++p) {
                well_potentials[w * np + p] = std::abs(potentials[p]);
            }
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    prepareTimeStep(const Simulator& ebos_simulator,
                    WellState& well_state)
    {
        const int nw = wells().number_of_wells;
        for (int w = 0; w < nw; ++w) {
            // after restarting, the well_controls can be modified while
            // the well_state still uses the old control index
            // we need to synchronize these two.
            resetWellControlFromState(well_state);

            if (wellCollection()->groupControlActive()) {
                WellControls* wc = wells().ctrls[w];
                WellNode& well_node = well_collection_->findWellNode(std::string(wells().name[w]));

                // handling the situation that wells do not have a valid control
                // it happens the well specified with GRUP and restarting due to non-convergencing
                // putting the well under group control for this situation
                int ctrl_index = well_controls_get_current(wc);

                const int group_control_index = well_node.groupControlIndex();
                if (group_control_index >= 0 && ctrl_index < 0) {
                    // put well under group control
                    well_controls_set_current(wc, group_control_index);
                    well_state.currentControls()[w] = group_control_index;
                }

                // Final step, update whehter the well is under group control or individual control
                // updated ctrl_index from the well control
                ctrl_index = well_controls_get_current(wc);
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
            if (well_collection_->requireWellPotentials()) {

                // calculate the well potentials
                setWellVariables(well_state);
                computeWellConnectionPressures(ebos_simulator, well_state);

                // To store well potentials for each well
                std::vector<double> well_potentials;
                computeWellPotentials(ebos_simulator, well_state, well_potentials);

                // update/setup guide rates for each well based on the well_potentials
                well_collection_->setGuideRatesWithPotentials(wellsPointer(), phase_usage_, well_potentials);
            }

            applyVREPGroupControl(well_state);

            if (!wellCollection()->groupControlApplied()) {
                wellCollection()->applyGroupControls();
            } else {
                wellCollection()->updateWellTargets(well_state.wellRates());
            }
        }

        // since the controls are all updated, we should update well_state accordingly
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells().ctrls[w];
            const int control = well_controls_get_current(wc);
            well_state.currentControls()[w] = control;
            updateWellStateWithTarget(wc, control, w, well_state);

            // The wells are not considered to be newly added
            // for next time step
            if (well_state.isNewWell(w) ) {
                well_state.setNewWell(w, false);
            }
        }  // end of for (int w = 0; w < nw; ++w)
    }





    template<typename TypeTag>
    WellCollection*
    StandardWellsDense<TypeTag>::
    wellCollection() const
    {
        return well_collection_;
    }




    template<typename TypeTag>
    const std::vector<double>&
    StandardWellsDense<TypeTag>::
    wellPerfEfficiencyFactors() const
    {
        return well_perforation_efficiency_factors_;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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

                rate_converter_.calcCoeff(well_rates, fipreg, convert_coeff);
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
                rate_converter_.calcCoeff(well_rates, fipreg, convert_coeff);
                std::copy(convert_coeff.begin(), convert_coeff.end(),
                          voidage_conversion_coeffs.begin() + np * w);
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateGroupControls(WellState& well_state) const
    {
        if (wellCollection()->groupControlActive()) {
            applyVREPGroupControl(well_state);
            wellCollection()->updateWellTargets(well_state.wellRates());

            // TODO: group control has to be applied in the level of the all wells
            // upate the well targets following group controls
            // it will not change the control mode, only update the targets
            for (int w = 0; w < number_of_wells_; ++w) {
                // TODO: check whether we need current argument in updateWellStateWithTarget
                // maybe there is some circumstances that the current is different from the one
                // in the WellState.
                // while probalby, the current argument can be removed
                const int current = well_state.currentControls()[w];
                well_container_[w]->updateWellStateWithTarget(current, well_state);
            }
        }
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    getBhp(const int wellIdx) const {
        // return well_container_(wellIdx)->getBhp();
        return 0.0; // TODO: for debugging
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    getQs(const int wellIdx, const int compIdx) const
    {
        // TODO: incoporate the change from the new PR to the getQs
	// in StandardWell
        return well_container_(wellIdx)->getQs(compIdx);
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    wellVolumeFraction(const int wellIdx, const int compIdx) const
    {
        const int nw = wells().number_of_wells;
        if (compIdx == Water) {
            return wellVariables_[WFrac * nw + wellIdx];
        }

        if (compIdx == Gas) {
            return wellVariables_[GFrac * nw + wellIdx];
        }

        if (has_solvent_ && compIdx == solventSaturationIdx) {
            return wellVariables_[SFrac * nw + wellIdx];
        }

        // Oil fraction
        EvalWell well_fraction = 1.0;
        if (active_[Water]) {
            well_fraction -= wellVariables_[WFrac * nw + wellIdx];
        }

        if (active_[Gas]) {
            well_fraction -= wellVariables_[GFrac * nw + wellIdx];
        }
        if (has_solvent_) {
            well_fraction -= wellVariables_[SFrac * nw + wellIdx];
        }
        return well_fraction;
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    wellVolumeFractionScaled(const int wellIdx, const int compIdx) const
    {
        const WellControls* wc = wells().ctrls[wellIdx];
        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {

            if (has_solvent_ && compIdx == solventSaturationIdx) {
                return wellVolumeFraction(wellIdx, compIdx);
            }
            const double* distr = well_controls_get_current_distr(wc);
            assert(compIdx < 3);
            if (distr[compIdx] > 0.) {
                return wellVolumeFraction(wellIdx, compIdx) / distr[compIdx];
            } else {
                // TODO: not sure why return EvalWell(0.) causing problem here
                // Probably due to the wrong Jacobians.
                return wellVolumeFraction(wellIdx, compIdx);
            }
        }
        std::vector<double> g = {1,1,0.01,0.01};
        return (wellVolumeFraction(wellIdx, compIdx) / g[compIdx]);
    }





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::EvalWell
    StandardWellsDense<TypeTag>::
    wellSurfaceVolumeFraction(const int well_index, const int compIdx) const
    {
        EvalWell sum_volume_fraction_scaled = 0.;
        const int numComp = numComponents();
        for (int idx = 0; idx < numComp; ++idx) {
            sum_volume_fraction_scaled += wellVolumeFractionScaled(well_index, idx);
        }

        assert(sum_volume_fraction_scaled.value() != 0.);

        return wellVolumeFractionScaled(well_index, compIdx) / sum_volume_fraction_scaled;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::RatioCheckTuple
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    typename StandardWellsDense<TypeTag>::RatioCheckTuple
    StandardWellsDense<TypeTag>::
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





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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

        case RESERVOIR_RATE: // intentional fall-through
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
                double original_rates_under_phase_control = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        original_rates_under_phase_control += xw.wellRates()[np * well_index + phase] * distr[phase];
                    }
                }

                if (original_rates_under_phase_control != 0.0 ) {
                    double scaling_factor = target / original_rates_under_phase_control;

                    for (int phase = 0; phase < np; ++phase) {
                        xw.wellRates()[np * well_index + phase] *= scaling_factor;
                    }
                } else { // scaling factor is not well defied when original_rates_under_phase_control is zero
                    // separating targets equally between phases under control
                    const double target_rate_divided = target / numPhasesWithTargetsUnderThisControl;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            xw.wellRates()[np * well_index + phase] = target_rate_divided / distr[phase];
                        } else {
                            // this only happens for SURFACE_RATE control
                            xw.wellRates()[np * well_index + phase] = target_rate_divided;
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
                xw.wellSolutions()[GFrac*nw + well_index] = g[Gas] * (xw.wellRates()[np*well_index + Gas] - xw.solventWellRate(well_index)) / tot_well_rate ;
            }
            if (has_solvent_) {
                xw.wellSolutions()[SFrac*nw + well_index] = g[Gas] * xw.solventWellRate(well_index) / tot_well_rate ;
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
                        xw.wellSolutions()[GFrac * nw + well_index] = 1.0 - wsolvent(well_index);
                        if (has_solvent_) {
                            xw.wellSolutions()[SFrac * nw + well_index] = wsolvent(well_index);
                        }
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





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    wellHasTHPConstraints(const int well_index) const
    {
        const WellControls* well_control = wells().ctrls[well_index];
        const int nwc = well_controls_get_num(well_control);
        for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
            if (well_controls_iget_type(well_control, ctrl_index) == THP) {
                return true;
            }
        }
        return false;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const EvalWell& bhp,
                            const int well_index,
                            std::vector<double>& well_flux) const
    {
        const int np = wells().number_of_phases;
        const int numComp = numComponents();
        well_flux.resize(np, 0.0);

        const bool allow_cf = well_container_[well_index]->crossFlowAllowed(ebosSimulator);
        for (int perf = wells().well_connpos[well_index]; perf < wells().well_connpos[well_index + 1]; ++perf) {
            const int cell_index = wells().well_cells[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_index, /*timeIdx=*/ 0));
            // flux for each perforation
            std::vector<EvalWell> cq_s(numComp, 0.0);
            std::vector<EvalWell> mob(numComp, 0.0);
            getMobility(ebosSimulator, well_index, perf, cell_index, mob);
            computeWellFlux(well_index, wells().WI[perf], intQuants, mob, bhp,
                            wellPerforationPressureDiffs()[perf], allow_cf, cq_s);

            for(int p = 0; p < np; ++p) {
                well_flux[p] += cq_s[p].value();
            }
        }
    }






    template<typename TypeTag>
    double
    StandardWellsDense<TypeTag>::
    mostStrictBhpFromBhpLimits(const int well_index) const
    {
        double bhp;

        // type of the well, INJECTOR or PRODUCER
        const WellType& well_type = wells().type[well_index];
        // initial bhp value, making the value not usable
        switch(well_type) {
        case INJECTOR:
            bhp = std::numeric_limits<double>::max();
            break;
        case PRODUCER:
            bhp = -std::numeric_limits<double>::max();
            break;
        default:
            OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << wells().name[well_index]);
        }

        // the well controls
        const WellControls* well_control = wells().ctrls[well_index];
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
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << wells().name[well_index]);
                } // end of switch
            }
        }

        return bhp;
    }





    template<typename TypeTag>
    std::vector<double>
    StandardWellsDense<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                const int well_index,
                                const double initial_bhp, // bhp from BHP constraints
                                const std::vector<double>& initial_potential) const
    {
        // TODO: pay attention to the situation that finally the potential is calculated based on the bhp control
        // TODO: should we consider the bhp constraints during the iterative process?
        const int np = wells().number_of_phases;

        assert( np == int(initial_potential.size()) );

        std::vector<double> potentials = initial_potential;
        std::vector<double> old_potentials = potentials; // keeping track of the old potentials

        double bhp = initial_bhp;
        double old_bhp = bhp;

        bool converged = false;
        const int max_iteration = 1000;
        const double bhp_tolerance = 1000.; // 1000 pascal

        int iteration = 0;

        while ( !converged && iteration < max_iteration ) {
            // for each iteration, we calculate the bhp based on the rates/potentials with thp constraints
            // with considering the bhp value from the bhp limits. At the beginning of each iteration,
            // we initialize the bhp to be the bhp value from the bhp limits. Then based on the bhp values calculated
            // from the thp constraints, we decide the effective bhp value for well potential calculation.
            bhp = initial_bhp;

            // the well controls
            const WellControls* well_control = wells().ctrls[well_index];
            // The number of the well controls/constraints
            const int nwc = well_controls_get_num(well_control);

            for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
                if (well_controls_iget_type(well_control, ctrl_index) == THP) {
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    const Opm::PhaseUsage& pu = phase_usage_;

                    if (active_[ Water ]) {
                        aqua = potentials[pu.phase_pos[ Water ] ];
                    }
                    if (active_[ Oil ]) {
                        liquid = potentials[pu.phase_pos[ Oil ] ];
                    }
                    if (active_[ Gas ]) {
                        vapour = potentials[pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(well_control, ctrl_index);
                    const double thp    = well_controls_iget_target(well_control, ctrl_index);
                    const double alq    = well_controls_iget_alq(well_control, ctrl_index);

                    // Calculating the BHP value based on THP
                    // TODO: check whether it is always correct to do calculation based on the depth of the first perforation.
                    const int first_perf = wells().well_connpos[well_index]; //first perforation

                    const WellType& well_type = wells().type[well_index];
                    if (well_type == INJECTOR) {
                        const double dp = wellhelpers::computeHydrostaticCorrection(
                                          wells(), well_index, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                          wellPerforationDensities()[first_perf], gravity_);
                        const double bhp_calculated = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                        // apply the strictest of the bhp controlls i.e. smallest bhp for injectors
                        if (bhp_calculated < bhp) {
                            bhp = bhp_calculated;
                        }
                    }
                    else if (well_type == PRODUCER) {
                        const double dp = wellhelpers::computeHydrostaticCorrection(
                                          wells(), well_index, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
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

            // there should be always some available bhp/thp constraints there
            if (std::isinf(bhp) || std::isnan(bhp)) {
                OPM_THROW(std::runtime_error, "Unvalid bhp value obtained during the potential calculation for well " << wells().name[well_index]);
            }

            converged = std::abs(old_bhp - bhp) < bhp_tolerance;

            computeWellRatesWithBhp(ebosSimulator, bhp, well_index, potentials);

            // checking whether the potentials have valid values
            for (const double value : potentials) {
                if (std::isinf(value) || std::isnan(value)) {
                    OPM_THROW(std::runtime_error, "Unvalid potential value obtained during the potential calculation for well " << wells().name[well_index]);
                }
            }

            if (!converged) {
                old_bhp = bhp;
                for (int p = 0; p < np; ++p) {
                    // TODO: improve the interpolation, will it always be valid with the way below?
                    // TODO: finding better paramters, better iteration strategy for better convergence rate.
                    const double potential_update_damping_factor = 0.001;
                    potentials[p] = potential_update_damping_factor * potentials[p] + (1.0 - potential_update_damping_factor) * old_potentials[p];
                    old_potentials[p] = potentials[p];
                }
            }

            ++iteration;
        }

        if (!converged) {
            OPM_THROW(std::runtime_error, "Failed in getting converged for the potential calculation for well " << wells().name[well_index]);
        }

        return potentials;
    }

    template<typename TypeTag>
    double
    StandardWellsDense<TypeTag>::
    wsolvent(const int well_index) const {

        if (!has_solvent_) {
            return 0.0;
        }

        // loop over all wells until we find the well with the matching name
        for (const auto&  well : wells_ecl_) {
            if (well->getStatus( current_timeIdx_ ) == WellCommon::SHUT) {
                continue;
            }

            WellInjectionProperties injection = well->getInjectionProperties(current_timeIdx_);
            if (injection.injectorType == WellInjector::GAS) {

                double solventFraction = well->getSolventFraction(current_timeIdx_);

                // Look until we find the correct well
                if (well->name() == wells().name[well_index]) {
                    return solventFraction;
                }
            }
        }
        // we didn't find it return 0;
        assert(false);
        return 0.0;
    }


    template<typename TypeTag>
    double
    StandardWellsDense<TypeTag>::
    wpolymer(const int well_index) const {

        if (!has_polymer_) {
            return 0.0;
        }

        // loop over all wells until we find the well with the matching name
        for (const auto&  well : wells_ecl_) {
            if (well->getStatus( current_timeIdx_ ) == WellCommon::SHUT) {
                continue;
            }

            WellInjectionProperties injection = well->getInjectionProperties(current_timeIdx_);
            WellPolymerProperties polymer = well->getPolymerProperties(current_timeIdx_);
            if (injection.injectorType == WellInjector::WATER) {

                double polymerFraction = polymer.m_polymerConcentration;

                // Look until we find the correct well
                if (well->name() == wells().name[well_index]) {
                    return polymerFraction;
                }
            }
        }
        // we didn't find it return 0;
        assert(false);
        return 0.0;
    }


    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setupCompressedToCartesian(const int* global_cell, int number_of_cells, std::map<int,int>& cartesian_to_compressed ) const
    {
        if (global_cell) {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
            }
        }
        else {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(i, i));
            }
        }

    }

    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeRepRadiusPerfLength(const Grid& grid)
    {

        // TODO, the function does not work for parallel running
        // to be fixed later.
        int number_of_cells = Opm::UgGridHelpers::numCells(grid);
        const int* global_cell = Opm::UgGridHelpers::globalCell(grid);
        const int* cart_dims = Opm::UgGridHelpers::cartDims(grid);
        auto cell_to_faces = Opm::UgGridHelpers::cell2Faces(grid);
        auto begin_face_centroids = Opm::UgGridHelpers::beginFaceCentroids(grid);

        if (wells_ecl_.size() == 0) {
            OPM_MESSAGE("No wells specified in Schedule section, "
                        "initializing no wells");
            return;
        }

        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];

        const size_t timeStep = current_timeIdx_;

        wells_rep_radius_.clear();
        wells_perf_length_.clear();
        wells_bore_diameter_.clear();

        wells_rep_radius_.reserve(nperf);
        wells_perf_length_.reserve(nperf);
        wells_bore_diameter_.reserve(nperf);

        std::map<int,int> cartesian_to_compressed;

        setupCompressedToCartesian(global_cell, number_of_cells,
                                    cartesian_to_compressed);

        int well_index = 0;

        for (auto wellIter= wells_ecl_.begin(); wellIter != wells_ecl_.end(); ++wellIter) {
             const auto* well = (*wellIter);

             if (well->getStatus(timeStep) == WellCommon::SHUT) {
                 continue;
             }
             {   // COMPDAT handling
                 const auto& completionSet = well->getCompletions(timeStep);
                 for (size_t c=0; c<completionSet.size(); c++) {
                     const auto& completion = completionSet.get(c);
                     if (completion.getState() == WellCompletion::OPEN) {
                         int i = completion.getI();
                         int j = completion.getJ();
                         int k = completion.getK();

                         const int* cpgdim = cart_dims;
                         int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                         std::map<int, int>::const_iterator cgit = cartesian_to_compressed.find(cart_grid_indx);
                         if (cgit == cartesian_to_compressed.end()) {
                             OPM_THROW(std::runtime_error, "Cell with i,j,k indices " << i << ' ' << j << ' '
                                       << k << " not found in grid (well = " << well->name() << ')');
                         }
                         int cell = cgit->second;

                         {
                             double radius = 0.5*completion.getDiameter();
                             if (radius <= 0.0) {
                                 radius = 0.5*unit::feet;
                                 OPM_MESSAGE("**** Warning: Well bore internal radius set to " << radius);
                             }

                             const std::array<double, 3> cubical =
                             WellsManagerDetail::getCubeDim<3>(cell_to_faces, begin_face_centroids, cell);

                             WellCompletion::DirectionEnum direction = completion.getDirection();

                             double re; // area equivalent radius of the grid block
                             double perf_length; // the length of the well perforation

                             switch (direction) {
                                 case Opm::WellCompletion::DirectionEnum::X:
                                     re = std::sqrt(cubical[1] * cubical[2] / M_PI);
                                     perf_length = cubical[0];
                                     break;
                                 case Opm::WellCompletion::DirectionEnum::Y:
                                     re = std::sqrt(cubical[0] * cubical[2] / M_PI);
                                     perf_length = cubical[1];
                                     break;
                                 case Opm::WellCompletion::DirectionEnum::Z:
                                     re = std::sqrt(cubical[0] * cubical[1] / M_PI);
                                     perf_length = cubical[2];
                                     break;
                                 default:
                                     OPM_THROW(std::runtime_error, " Dirtecion of well is not supported ");
                             }

                             double repR = std::sqrt(re * radius);
                             wells_rep_radius_.push_back(repR);
                             wells_perf_length_.push_back(perf_length);
                             wells_bore_diameter_.push_back(2. * radius);
                         }
                     } else {
                         if (completion.getState() != WellCompletion::SHUT) {
                             OPM_THROW(std::runtime_error, "Completion state: " << WellCompletion::StateEnum2String( completion.getState() ) << " not handled");
                         }
                     }

                 }
            }
            well_index++;
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeAverageFormationFactor(Simulator& ebosSimulator,
                                  std::vector<double>& B_avg) const
    {
        const int np = numPhases();
        const int numComp = numComponents();

        const auto& grid = ebosSimulator.gridManager().grid();
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

            for ( int phaseIdx = 0; phaseIdx < np; ++phaseIdx )
            {
                auto& B  = B_avg[ phaseIdx ];
                const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phaseIdx);

                B += 1 / fs.invB(ebosPhaseIdx).value();
            }
            if (has_solvent_) {
                auto& B  = B_avg[solventSaturationIdx];
                B += 1 / intQuants.solventInverseFormationVolumeFactor().value();
            }
        }

        // compute global average
        grid.comm().sum(B_avg.data(), B_avg.size());
        for(auto& bval: B_avg)
        {
            bval/=global_nc_;
        }
    }




    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    outputWellState(const WellState& well_state) const
    {
        std::cout << " output the bhp " << std::endl;
        for (const double bhp : well_state.bhp()) {
            std::cout << bhp << " ";
        }
        std::cout << std::endl;

        std::cout << " output the well rates " << std::endl;
        for (const double rate : well_state.wellRates()) {
            std::cout << rate << " ";
        }
        std::cout << std::endl;

        std::cout << " output the wellSolutions " << std::endl;
        for (const double solution : well_state.wellSolutions()) {
            std::cout << solution << " ";
        }
        std::cout << std::endl;
    }




} // namespace Opm
