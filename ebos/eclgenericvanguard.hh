// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::EclBaseVanguard
 */
#ifndef EWOMS_ECL_GENERIC_VANGUARD_HH
#define EWOMS_ECL_GENERIC_VANGUARD_HH

#include <opm/grid/common/GridEnums.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
#include <dune/common/parallel/communication.hh>
#else
#include <dune/common/parallel/collectivecommunication.hh>
#endif

#include <array>
#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


namespace Opm {

namespace Action { class State; }
class Deck;
class EclipseState;
class ErrorGuard;
struct NumericalAquiferCell;
class ParseContext;
class Schedule;
class Python;
class SummaryConfig;
class SummaryState;
class UDQState;
class WellTestState;

class EclGenericVanguard {
public:
    using ParallelWellStruct = std::vector<std::pair<std::string,bool>>;

    /*!
     * \brief Constructor.
     * \details Needs to be in compile unit.
     */
    EclGenericVanguard();

    /*!
     * \brief Destructor.
     * \details Empty, but needs to be in compile unit.
     */
    ~EclGenericVanguard();

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static std::string canonicalDeckPath(const std::string& caseName);

    /*!
     * \brief Creates an Opm::parseContext object assuming that the parameters are ready.
     */
    static std::unique_ptr<ParseContext> createParseContext(const std::string& ignoredKeywords,
                                                            bool eclStrictParsing);

    /*!
     * \brief Set the wall time which was spend externally to set up the external data structures
     *
     * i.e., the objects specified via the other setExternal*() methods.
     */
    static void setExternalSetupTime(double t)
    { externalSetupTime_ = t; }

    /*!
     * \brief Returns the wall time required to set up the simulator before it was born.
     */
    static double externalSetupTime()
    { return externalSetupTime_; }

    /*!
     * \brief Set the Opm::ParseContext object which ought to be used for parsing the deck and creating the Opm::EclipseState object.
     */
    static void setExternalParseContext(std::unique_ptr<ParseContext> parseContext);

    /*!
     * \brief Set the Opm::ErrorGuard object which ought to be used for parsing the deck and creating the Opm::EclipseState object.
     */
    static void setExternalErrorGuard(std::unique_ptr<ErrorGuard> errorGuard);

    /*!
     * \brief Set the Opm::Deck object which ought to be used when the simulator vanguard
     *        is instantiated.
     *
     * This is basically an optimization: In cases where the ECL input deck must be
     * examined to decide which simulator ought to be used, this avoids having to parse
     * the input twice. When this method is used, the caller is responsible for lifetime
     * management of these two objects, i.e., they are not allowed to be deleted as long
     * as the simulator vanguard object is alive.
     */
    static void setExternalDeck(std::shared_ptr<Deck> deck);
    static void setExternalDeck(std::unique_ptr<Deck> deck);

    /*!
     * \brief Set the Opm::EclipseState object which ought to be used when the simulator
     *        vanguard is instantiated.
     */
    static void setExternalEclState(std::shared_ptr<EclipseState> eclState);
    static void setExternalEclState(std::unique_ptr<EclipseState> eclState);

    /*!
     * \brief Set the schedule object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSchedule(std::shared_ptr<Schedule> schedule);
    static void setExternalSchedule(std::unique_ptr<Schedule> schedule);

    /*!
     * \brief Set the summary configuration object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSummaryConfig(std::shared_ptr<SummaryConfig> summaryConfig);
    static void setExternalSummaryConfig(std::unique_ptr<SummaryConfig> summaryConfig);

    static void setExternalUDQState(std::unique_ptr<UDQState> udqState);
    static void setExternalActionState(std::unique_ptr<Action::State> actionState);
    static void setExternalWTestState(std::unique_ptr<WellTestState> wtestState);


    /*!
     * \brief Return a reference to the parsed ECL deck.
     */
    const Deck& deck() const
    { return *deck_; }

    Deck& deck()
    { return *deck_; }

    /*!
     * \brief Return a reference to the internalized ECL deck.
     */
    const EclipseState& eclState() const
    { return *eclState_; }

    EclipseState& eclState()
    { return *eclState_; }

    /*!
     * \brief Return a reference to the object that managages the ECL schedule.
     */
    const Schedule& schedule() const
    { return *eclSchedule_; }

    Schedule& schedule()
    { return *eclSchedule_; }

    /*!
     * \brief Return a reference to the object that determines which quantities ought to
     *        be put into the ECL summary output.
     */
    const SummaryConfig& summaryConfig() const
    { return *eclSummaryConfig_; }

    /*!
    * \brief Returns the summary state
    *
    * The summary state is a small container object for
    * computed, ready to use summary values. The values will typically be used by
    * the UDQ, WTEST and ACTIONX calculations.
    */
    SummaryState& summaryState()
    { return *summaryState_; }

    const SummaryState& summaryState() const
    { return *summaryState_; }

    /*!
     * \brief Returns the action state
     *
     * The ActionState keeps track of how many times the various actions have run.
     */
    Action::State& actionState()
    { return *actionState_; }

    const Action::State& actionState() const
    { return *actionState_; }

    /*!
     * \brief Returns the udq state
     *
     * The UDQState keeps track of the result of UDQ evaluations.
     */
    UDQState& udqState()
    { return *udqState_; }

    const UDQState& udqState() const
    { return *udqState_; }

    WellTestState transferWTestState() {
        return *this->wtestState_.release();
    }


    /*!
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string& caseName() const
    { return caseName_; }

    /*!
     * \brief Parameter deciding the edge-weight strategy of the load balancer.
     */
    Dune::EdgeWeightMethod edgeWeightsMethod() const
    { return edgeWeightsMethod_; }

    /*!
     * \brief Parameter that decide if cells owned by rank are ordered before ghost cells.
     */
    bool ownersFirst() const
    { return ownersFirst_; }

    /*!
     * \brief Parameter that decides if partitioning for parallel runs
     *        should be performed on a single process only.
     */
    bool serialPartitioning() const
    { return serialPartitioning_; }

    /*!
     * \brief Parameter that sets the zoltan imbalance tolarance.
     */
    double zoltanImbalanceTol() const
    { return zoltanImbalanceTol_; }

    /*!
     * \brief Whether perforations of a well might be distributed.
     */
    bool enableDistributedWells() const
    { return enableDistributedWells_; }

    /*!
     * \brief Returns vector with name and whether the has local perforated cells
     *        for all wells.
     *
     * Will only have usable values for CpGrid.
     */
    const ParallelWellStruct& parallelWells() const
    { return parallelWells_; }

    //! \brief Set global communication.
    static void setCommunication(std::unique_ptr<Opm::Parallel::Communication> comm)
    { comm_ = std::move(comm); }

    //! \brief Obtain global communicator.
    static Parallel::Communication& comm()
    {
        assert(comm_);
        return *comm_;
    }

protected:
    void updateOutputDir_(std::string outputDir,
                          bool enableEclCompatFile);

    bool drsdtconEnabled() const;

    std::unordered_map<std::size_t, const NumericalAquiferCell*> allAquiferCells() const;

    void init();

    static double externalSetupTime_;
    static std::unique_ptr<ParseContext> externalParseContext_;
    static std::unique_ptr<ErrorGuard> externalErrorGuard_;

    // These variables may be owned by both Python and the simulator
    static std::shared_ptr<Deck> externalDeck_;
    static std::shared_ptr<EclipseState> externalEclState_;
    static std::shared_ptr<Schedule> externalEclSchedule_;
    static std::shared_ptr<SummaryConfig> externalEclSummaryConfig_;

    static bool externalDeckSet_;
    static std::unique_ptr<UDQState> externalUDQState_;
    static std::unique_ptr<Action::State> externalActionState_;
    static std::unique_ptr<WellTestState> externalWTestState_;
    static std::unique_ptr<Parallel::Communication> comm_;

    std::string caseName_;
    std::string fileName_;
    Dune::EdgeWeightMethod edgeWeightsMethod_;
    bool ownersFirst_;
    bool serialPartitioning_;
    double zoltanImbalanceTol_;
    bool enableDistributedWells_;
    std::string ignoredKeywords_;
    bool eclStrictParsing_;
    std::optional<int> outputInterval_;
    bool useMultisegmentWell_;
    bool enableExperiments_;

    std::unique_ptr<SummaryState> summaryState_;
    std::unique_ptr<UDQState> udqState_;
    std::unique_ptr<Action::State> actionState_;

    // Observe that this instance is handled differently from the other state
    // variables, it will only be initialized for a restart run. While
    // initializing a restarted run this instance is transferred to the WGState
    // member in the well model.
    std::unique_ptr<WellTestState> wtestState_;

    // these attributes point  either to the internal  or to the external version of the
    // parser objects.
    std::unique_ptr<ParseContext> parseContext_;
    std::unique_ptr<ErrorGuard> errorGuard_;
    std::shared_ptr<Python> python;
    // These variables may be owned by both Python and the simulator
    std::shared_ptr<Deck> deck_;
    std::shared_ptr<EclipseState> eclState_;
    std::shared_ptr<Schedule> eclSchedule_;
    std::shared_ptr<SummaryConfig> eclSummaryConfig_;

    /*! \brief Information about wells in parallel
     *
     * For each well in the model there is an entry with its name
     * and a boolean indicating whether it perforates local cells.
     */
    ParallelWellStruct parallelWells_;
};

} // namespace Opm

#endif
