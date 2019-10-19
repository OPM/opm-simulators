// #include <sstream>
#include <opm/grid/utility/cartesianToCompressed.hpp>
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
    if(aquiferCarterTracyActive())
    {
      for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
      {
        aquifer->initialSolutionApplied();
      }
    }
    if(aquiferFetkovichActive())
    {
      for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
      {
        aquifer->initialSolutionApplied();
      }
    }
  }

  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>::beginEpisode()
  { }

  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>::beginIteration()
  { }

  template<typename TypeTag>
  void BlackoilAquiferModel<TypeTag>:: beginTimeStep()
  {
    if(aquiferCarterTracyActive())
    {
      for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
      {
        aquifer->beginTimeStep();
      }
    }
    if(aquiferFetkovichActive())
    {
      for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
      {
        aquifer->beginTimeStep();
      }
    }
  }

  template<typename TypeTag>
  template<class Context>
  void BlackoilAquiferModel<TypeTag>:: addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
  {
    if(aquiferCarterTracyActive())
    {
      for (auto& aquifer : aquifers_CarterTracy)
      {
        aquifer.addToSource(rates, context, spaceIdx, timeIdx);
      }
    }
    if(aquiferFetkovichActive())
    {
      for (auto& aquifer : aquifers_Fetkovich)
      {
        aquifer.addToSource(rates, context, spaceIdx, timeIdx);
      }
    }
  }

  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>::endIteration()
  { }

  template<typename TypeTag>
  void BlackoilAquiferModel<TypeTag>:: endTimeStep()
  {
    if(aquiferCarterTracyActive())
    {
      for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
      {
        aquifer->endTimeStep();
      }
    }
    if(aquiferFetkovichActive())
    {
      for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
      {
        aquifer->endTimeStep();
      }
    }
  }
  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>::endEpisode()
  { }

    template <typename TypeTag>
    template <class Restarter>
    void
    BlackoilAquiferModel<TypeTag>::serialize(Restarter& /* res */)
    {
        // TODO (?)
        throw std::logic_error("BlackoilAquiferModel::serialize() is not yet implemented");
    }

    template<typename TypeTag>
    template <class Restarter>
    void
    BlackoilAquiferModel<TypeTag>::deserialize(Restarter& /* res */)
    {
        // TODO (?)
        throw std::logic_error("BlackoilAquiferModel::deserialize() is not yet implemented");
    }

  // Initialize the aquifers in the deck
  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>:: init()
  {
    const auto& deck = this->simulator_.vanguard().deck();
    if (deck.hasKeyword("AQUCT")) {
      //updateConnectionIntensiveQuantities();
      const auto& eclState = this->simulator_.vanguard().eclState();

      // Get all the carter tracy aquifer properties data and put it in aquifers vector
      const AquiferCT aquiferct = AquiferCT(eclState,deck);
      const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

      std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
      std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

      assert( aquifersData.size() == aquifer_connection.size() );
      const auto& ugrid = simulator_.vanguard().grid();
      const auto& gridView = simulator_.gridView();
      const int number_of_cells = gridView.size(0);

      cartesian_to_compressed_ = cartesianToCompressed(number_of_cells,
                                                       Opm::UgGridHelpers::globalCell(ugrid));

      for (size_t i = 0; i < aquifersData.size(); ++i)
      {
        aquifers_CarterTracy.push_back(
          AquiferCarterTracy<TypeTag>  (aquifer_connection.at(i), cartesian_to_compressed_, this->simulator_ , aquifersData.at(i))
        );
      }
    }
    if(deck.hasKeyword("AQUFETP"))
    {
      //updateConnectionIntensiveQuantities();
      const auto& eclState = this->simulator_.vanguard().eclState();

      // Get all the carter tracy aquifer properties data and put it in aquifers vector
      const Aquifetp aquifetp = Aquifetp(deck);
      const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

      const std::vector<Aquifetp::AQUFETP_data>& aquifersData = aquifetp.getAquifers();
      const std::vector<Aquancon::AquanconOutput>& aquifer_connection = aquifer_connect.getAquOutput();

      // it can happen that some aquifer defined while there are no valid connections associated with it
      assert( aquifersData.size() >= aquifer_connection.size() );
      std::cout << " aquifersData.size() " << aquifersData.size() << std::endl;
      std::cout << " aquifer_connection.size() " << aquifer_connection.size() << std::endl;
      std::cout << " aquifer IDs in aquifersData ";
      for (const auto& a : aquifersData) {
          std::cout << " " << a.aquiferID;
      }
      std::cout << std::endl;

      std::cout << " aquifer IDs in aquifer_connection ";
      for (const auto& c : aquifer_connection) {
          std::cout << " " << c.aquiferID;
      }
      std::cout << std::endl;

      const auto& ugrid = simulator_.vanguard().grid();
      const auto& gridView = simulator_.gridView();
      const int number_of_cells = gridView.size(0);

      cartesian_to_compressed_ = cartesianToCompressed(number_of_cells,
                                                       Opm::UgGridHelpers::globalCell(ugrid));

      size_t idx_aquifer_data = 0;
      for (size_t i = 0; i < aquifer_connection.size(); ++i, ++idx_aquifer_data) {
          const int aquifer_id = aquifer_connection[i].aquiferID;

          for (; idx_aquifer_data < aquifersData.size(); ++idx_aquifer_data) {
              if (aquifersData[idx_aquifer_data].aquiferID == aquifer_id) {
                  break;
              }
          }

          if (idx_aquifer_data == aquifersData.size()) {
              std::ostringstream msg;
              msg << " Connections for aquifer " << aquifer_id << " are defined, while there are no AQUFETP keyword"
                  << " specified for this aquifer ";
              throw std::logic_error(msg.str());
          }

          aquifers_Fetkovich.push_back(
            AquiferFetkovich<TypeTag> (aquifer_connection.at(i), cartesian_to_compressed_, this->simulator_ , aquifersData.at(idx_aquifer_data)) );

          std::cout << " aquifer " << aquifer_id << " connection idx " << i << " idx_aquifer_data " << idx_aquifer_data << std::endl;
      }
    }
  }
  template<typename TypeTag>
  bool
  BlackoilAquiferModel<TypeTag>:: aquiferActive() const
  {
    return (aquiferCarterTracyActive() || aquiferFetkovichActive());
  }
  template<typename TypeTag>
  bool
  BlackoilAquiferModel<TypeTag>:: aquiferCarterTracyActive() const
  {
    return !aquifers_CarterTracy.empty();
  }
  template<typename TypeTag>
  bool
  BlackoilAquiferModel<TypeTag>:: aquiferFetkovichActive() const
  {
    return !aquifers_Fetkovich.empty();
  }
} // namespace Opm
