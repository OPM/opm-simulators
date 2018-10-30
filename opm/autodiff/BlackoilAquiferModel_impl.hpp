namespace Opm {


    template<typename TypeTag>
    BlackoilAquiferModel<TypeTag>::
    BlackoilAquiferModel(Simulator& simulator)
        : simulator_(simulator)
    {
        init();
    }

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
