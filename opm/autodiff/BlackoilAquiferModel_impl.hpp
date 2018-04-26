namespace Opm {


    template<typename TypeTag>
    BlackoilAquiferModel<TypeTag>::
    BlackoilAquiferModel(Simulator& ebosSimulator)
        : ebosSimulator_(ebosSimulator)
    {
        init();
    }


    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: timeStepSucceeded(const SimulatorTimerInterface& timer)
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            aquifer->afterTimeStep(timer);
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::
    assemble( const SimulatorTimerInterface& timer,
              const int iterationIdx                )
    {
        // We need to update the reservoir pressures connected to the aquifer
        updateConnectionIntensiveQuantities();

        if (iterationIdx == 0) {
            // We can do the Table check and coefficients update in this function
            // For now, it does nothing!
            prepareTimeStep(timer);
        }

        assembleAquiferEq(timer);
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

    // Protected function which calls the individual aquifer models
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: assembleAquiferEq(const SimulatorTimerInterface& timer)
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer)
        {
            aquifer->assembleAquiferEq(timer);
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
            aquifer->beforeTimeStep(timer);
        }
    }

    // Initialize the aquifers in the deck
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: init()
    {
        updateConnectionIntensiveQuantities();
        const auto& deck = ebosSimulator_.vanguard().deck();
        const auto& eclState = ebosSimulator_.vanguard().eclState();
        
        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        const AquiferCT aquiferct = AquiferCT(eclState,deck);
        const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
        std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

        assert( aquifersData.size() == aquifer_connection.size() );


        for (size_t i = 0; i < aquifersData.size(); ++i)
        {
            aquifers_.push_back( 
                                  AquiferCarterTracy<TypeTag> (aquifersData.at(i), aquifer_connection.at(i), ebosSimulator_)
                               );
        }
    }

} // namespace Opm