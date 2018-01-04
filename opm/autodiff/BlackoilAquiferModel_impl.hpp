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
        const auto& eclState = ebosSimulator_.vanguard().eclState();
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

    }

    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: timeStepSucceeded(const SimulatorTimerInterface& timer)
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            aquifer->after_time_step(timer);
        }
    }

    // called at the beginning of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: beginReportStep(const int time_step) 
    {

    }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: endReportStep() 
    {

    }

    // Get the last report step
    template<typename TypeTag>
    const SimulatorReport& 
    BlackoilAquiferModel<TypeTag>:: lastReport() const 
    {
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

        if (param_.solve_aquifereq_initially_ && iterationIdx == 0) {
            // solve the aquifer equations as a pre-processing step
            last_report_ = solveAquiferEq(timer);
        }
        assembleAquiferEq(timer);
        last_report_.converged = true;
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
        const auto& pu = phase_usage_;
        return pu.num_phases;
    }

    // Protected function which calls the individual aquifer models
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: assembleAquiferEq(const SimulatorTimerInterface& timer)
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
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
    BlackoilAquiferModel<TypeTag>:: init(const Simulator& ebosSimulator, std::vector< AquiferCarterTracy<TypeTag> >& aquifers)
    {
        updateConnectionIntensiveQuantities();
        const auto& deck = ebosSimulator.vanguard().deck();
        const auto& eclState = ebosSimulator.vanguard().eclState();
        
        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        AquiferCT aquiferct = AquiferCT(eclState,deck);
        Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
        std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

        assert( aquifersData.size() == aquifer_connect.size() );


        for (int i = 0; i < aquifersData.size(); ++i)
        {
            aquifers.push_back( 
                                 AquiferCarterTracy<TypeTag> (aquifersData.at(i), aquifer_connection.at(i), numComponents(), gravity_, ebosSimulator_) 
                              );
        }
    }

} // namespace Opm