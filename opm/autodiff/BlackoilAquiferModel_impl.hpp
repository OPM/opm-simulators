namespace Opm {


    template<typename TypeTag>
    BlackoilAquiferModel<TypeTag>::
    BlackoilAquiferModel(Simulator& simulator)
        : simulator_(simulator)
    {
        init();
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::initialSolutionApplied()
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer) {
            aquifer->initialSolutionApplied();
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::beginEpisode()
    { }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::beginTimeStep()
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer) {
            aquifer->beginTimeStep();
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::beginIteration()
    { }

    template<typename TypeTag>
    template <class Context>
    void
    BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                               const Context& context,
                                               unsigned spaceIdx,
                                               unsigned timeIdx) const
    {
        for (auto& aquifer: aquifers_) {
            aquifer.addToSource(rates, context, spaceIdx, timeIdx);
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::endIteration()
    { }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::endTimeStep()
    {
        for (auto aquifer = aquifers_.begin(); aquifer != aquifers_.end(); ++aquifer) {
            aquifer->endTimeStep();
        }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::endEpisode()
    { }

    // Initialize the aquifers in the deck
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: init()
    {
        const auto& deck = this->simulator_.vanguard().deck();

        if ( !deck.hasKeyword("AQUCT") ) {
            return ;
        }

        //updateConnectionIntensiveQuantities();
        const auto& eclState = this->simulator_.vanguard().eclState();

        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        const AquiferCT aquiferct = AquiferCT(eclState,deck);
        const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
        std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

        assert( aquifersData.size() == aquifer_connection.size() );


        for (size_t i = 0; i < aquifersData.size(); ++i)
        {
            aquifers_.push_back(
                                  AquiferCarterTracy<TypeTag> (aquifersData.at(i), aquifer_connection.at(i), this->simulator_)
                               );
        }
    }

    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: aquiferActive() const
    {
        return !aquifers_.empty();
    }

} // namespace Opm
