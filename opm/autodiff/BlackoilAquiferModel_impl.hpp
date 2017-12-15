namespace Opm {


    template<typename TypeTag>
    BlackoilAquiferModel<TypeTag>::
    BlackoilAquiferModel(Simulator& ebosSimulator,
                      const ModelParameters& param,
                      const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , param_(param)
        , terminal_output_(terminal_output)
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
    {
        const auto& eclState = ebosSimulator_.gridManager().eclState();
        phase_usage_ = phaseUsageFromDeck(eclState);

        active_.resize(phase_usage_.MaxNumPhases, false);
        for (int p = 0; p < phase_usage_.MaxNumPhases; ++p) {
            active_[ p ] = phase_usage_.phase_used[ p ] != 0;
        }

        const auto& gridView = ebosSimulator_.gridView();

        // calculate the number of elements of the compressed sequential grid. this needs
        // to be done in two steps because the dune communicator expects a reference as
        // argument for sum()
        number_of_cells_ = gridView.size(/*codim=*/0);
        global_nc_ = gridView.comm().sum(number_of_cells_);
        gravity_ = ebosSimulator_.problem().gravity()[2];
        init(ebosSimulator_, aquifers_);
    }



    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: beginTimeStep() 
    {
        // Right now it doesn't do shit.
    }

    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: timeStepSucceeded()
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            aquifer->after_time_step();
        }
    }

    // called at the beginning of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: beginReportStep(const int time_step) 
    {
        // Right now it doesn't do shit.
    }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: endReportStep() 
    {
        // Right now it just spits out the constants for each aquifers
        // We are using the simple integer indexing for the aquifers
        for (int i = 0; i < numAquifers(); ++i)
        {
            std::cout << "Aquifer[" << i << "]"
                      << " : Tc = " << aquifers()[i].time_constant()
                      << ", beta = " << aquifers()[i].aquifer_influx_constant() << std::endl;
        }
    }

    // Get the last report step
    template<typename TypeTag>
    const SimulatorReport& 
    BlackoilAquiferModel<TypeTag>:: lastReport() const 
    {
        for (auto i = aquifers_.begin(); i != aquifers_.end(); ++i){
            (*i).print_private_members();
        }
        return last_report_;
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    assemble( const SimulatorTimerInterface& timer,
              const int iterationIdx                )
    {
        last_report_ = SimulatorReport();
        
        // We need to update the reservoir pressures connected to the aquifer
        updateConnectionIntensiveQuantities();

        if (iterationIdx == 0) {
            // We can do the Table check and coefficients update in this function
            // For now, it does nothing!
            prepareTimeStep(timer);
        }

        if (iterationIdx == 0) {
            calculateExplicitQuantities();
        }

        if (param_.solve_aquifereq_initially_ && iterationIdx == 0) {
            // solve the aquifer equations as a pre-processing step
            last_report_ = solveAquiferEq(timer);
        }

        assembleAquiferEq(timer);

        last_report_.converged = true;
    }

    // Protected function: Update the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: updatePrimaryVariables() 
    {
        // Right now it doesn't do shit.
    }

    // Protected function: Init the primary variables
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: initPrimaryVariablesEvaluation() const 
    {
        // Right now it doesn't do shit.
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: updateConnectionIntensiveQuantities() const
    {
        ElementContext elemCtx(ebosSimulator_);
        const auto& gridView = ebosSimulator_.gridView();
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            elemCtx.updatePrimaryStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: calculateExplicitQuantities()
    {
        // for (auto aqui = aquifers_.begin(); aqui!= aquifers_.end(); ++aqui)
        // {
        //     std::cout << "calculateExplicitQuantities: Aquifer id = " << aqui->aquiferID() << std::endl;
        //     aqui->calculateExplicitQuantities(ebosSimulator_);
        // }
    }


    template<typename TypeTag>
    SimulatorReport 
    BlackoilAquiferModel<TypeTag>:: solveAquiferEq(const SimulatorTimerInterface& timer)
    {
        // We need to solve the equilibrium equation first to
        // obtain the initial pressure of water in the aquifer
        SimulatorReport report = SimulatorReport();
        return report;
    }

    // Protected function: Return number of components in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numComponents() const
    {
        if (numPhases() == 2) {
            return 2;
        }
        int numComp = FluidSystem::numComponents;
        if (has_solvent_) {
            numComp ++;
        }

        return numComp;
    }

    // Protected function: Return number of aquifers in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numAquifers() const
    {
        return aquifers_.size();
    }

    // Protected function: Return number of phases in the model.
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: numPhases() const
    {
        // Not implemented yet!!!!!!!!!!!!
        const auto& pu = phase_usage_;
        return pu.num_phases;
    }


    // Protected function: returns the phase index in ebos
    template<typename TypeTag>
    int
    BlackoilAquiferModel<TypeTag>:: flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const auto& pu = phase_usage_;
        if (active_[Water] && pu.phase_pos[Water] == phaseIdx)
            return FluidSystem::waterPhaseIdx;
        if (active_[Oil] && pu.phase_pos[Oil] == phaseIdx)
            return FluidSystem::oilPhaseIdx;
        if (active_[Gas] && pu.phase_pos[Gas] == phaseIdx)
            return FluidSystem::gasPhaseIdx;

        assert(phaseIdx < 3);
        // for other phases return the index
        return phaseIdx;
    }

    // Protected function which calls the individual aquifer models
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: assembleAquiferEq(const SimulatorTimerInterface& timer)
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            std::cout << "assembleAquiferEq: Aquifer id = " << aquifer->aquiferID() << std::endl;
            aquifer->assembleAquiferEq(ebosSimulator_, timer);
        }
    }

    // Protected function
    // some preparation work, mostly related to group control and RESV,
    // at the beginning of each time step (Not report step)
    template<typename TypeTag>
    void BlackoilAquiferModel<TypeTag>:: prepareTimeStep(const SimulatorTimerInterface& timer)
    {
        // Here we can ask each carter tracy aquifers to get the current previous time step's pressure
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            aquifer->before_time_step(ebosSimulator_, timer);
        }
    }

    // Protected function: Returns a reference to the aquifers members in the model
    template<typename TypeTag>
    const std::vector< AquiferCarterTracy<TypeTag> >&
    BlackoilAquiferModel<TypeTag>:: aquifers()
    {
        return aquifers_;
    }


    // Initialize the aquifers in the deck
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: init(const Simulator& ebosSimulator, std::vector< AquiferCarterTracy<TypeTag> >& aquifers)//, std::vector< AquiferCarterTracy<TypeTag> >& aquifers)
    {
        const auto& deck = ebosSimulator.gridManager().deck();
        const auto& eclState = ebosSimulator.gridManager().eclState();
        
        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        AquiferCT aquiferct = AquiferCT(eclState,deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
        std::vector<AquiferCT::AQUANCON_data> aquanconData = aquiferct.getAquancon();

        // for (auto aquiferData = aquifersData.begin(); aquiferData != aquifersData.end(); ++aquiferData)
        // {
            
        // }

        auto ita = aquifersData.cbegin();
        auto f_lambda = [&] (AquiferCT::AQUANCON_data i) {
            aquifers.push_back( AquiferCarterTracy<TypeTag> (*ita++, i, numComponents(), gravity_ ) );
        };
        std::for_each( aquanconData.cbegin(), aquanconData.cend(), f_lambda );
    }

    // Begin the hack to initialize the aquifers in the deck
    template<typename TypeTag>
    std::vector< AquiferCarterTracy<TypeTag> >
    BlackoilAquiferModel<TypeTag>:: hack_init(const Simulator& ebosSimulator)//, std::vector< AquiferCarterTracy<TypeTag> >& aquifers)
    {
        std::vector< AquiferCarterTracy<TypeTag> > aquifers;
        /** Begin hack!!!!! */
        const auto& deck = ebosSimulator.gridManager().deck();
        const auto& eclState = ebosSimulator.gridManager().eclState();
        
        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        AquiferCT aquiferct = AquiferCT(eclState,deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();

        for (auto aquiferData = aquifersData.begin(); aquiferData != aquifersData.end(); ++aquiferData)
        {
            aquifers.push_back( AquiferCarterTracy<TypeTag> (*aquiferData, numComponents(), gravity_ ) );
        }
    }

} // namespace Opm