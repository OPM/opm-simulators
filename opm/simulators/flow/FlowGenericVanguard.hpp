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
 * \copydoc Opm::FlowBaseVanguard
 */
#ifndef OPM_FLOW_GENERIC_VANGUARD_HPP
#define OPM_FLOW_GENERIC_VANGUARD_HPP

#include <dune/common/parallel/communication.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/common/GridEnums.hpp>

#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cassert>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm::Parameters {

struct AllowDistributedWells { static constexpr bool value = false; };
struct AllowSplittingInactiveWells { static constexpr bool value = true; };

struct EclOutputInterval { static constexpr int value = -1; };
struct EdgeWeightsMethod  { static constexpr auto value = "transmissibility"; };
struct EnableDryRun { static constexpr auto value = "auto"; };
struct EnableEclOutput { static constexpr auto value = true; };
struct EnableOpmRstFile { static constexpr bool value = false; };
struct ExternalPartition { static constexpr auto* value = ""; };

template<class Scalar>
struct ImbalanceTol { static constexpr Scalar value = 1.1; };

struct IgnoreKeywords { static constexpr auto value = ""; };
struct InputSkipMode { static constexpr auto value = "100"; };
struct MetisParams { static constexpr auto value = "default"; };

#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
struct NumJacobiBlocks { static constexpr int value = 0; };
#endif // HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA

struct OwnerCellsFirst { static constexpr bool value = true; };
struct ParsingStrictness { static constexpr auto value = "normal"; };
struct ActionParsingStrictness { static constexpr auto value = "normal"; };

/// 0: simple, 1: Zoltan, 2: METIS, 3: Zoltan with a all cells of a well
/// represented by one vertex in the graph, see GridEnums.hpp
struct PartitionMethod { static constexpr auto value = "zoltanwell"; };
struct AddCorners { static constexpr bool value = false; };
struct NumOverlap { static constexpr int value = 1; };

struct SchedRestart{ static constexpr bool value = false; };
struct SerialPartitioning{ static constexpr bool value = false; };

template<class Scalar>
struct ZoltanImbalanceTol { static constexpr Scalar value = 1.1; };

struct ZoltanPhgEdgeSizeThreshold { static constexpr auto value = 0.35; };

struct ZoltanParams { static constexpr auto value = "graph"; };

} // namespace Opm::Parameters

namespace Opm {

namespace Action { class State; }
class Deck;
class EclipseState;
struct NumericalAquiferCell;
class ParseContext;
class Schedule;
class Python;
class SummaryConfig;
class SummaryState;
class UDQState;
class WellTestState;

class FlowGenericVanguard {
public:
    using ParallelWellStruct = std::vector<std::pair<std::string,bool>>;

    struct SimulationModelParams {
        double setupTime_;
        std::unique_ptr<UDQState> udqState_;
        std::unique_ptr<Action::State> actionState_;
        std::unique_ptr<SummaryState> summaryState_;
        std::unique_ptr<WellTestState> wtestState_;
        std::shared_ptr<EclipseState> eclState_;
        std::shared_ptr<Schedule> eclSchedule_;
        std::shared_ptr<SummaryConfig> eclSummaryConfig_;
    };

    static SimulationModelParams modelParams_;

    /*!
     * \brief Constructor.
     * \details Needs to be in compile unit.
     */
    FlowGenericVanguard();
    explicit FlowGenericVanguard(SimulationModelParams&& params);

    /*!
     * \brief Destructor.
     * \details Empty, but needs to be in compile unit.
     */
    ~FlowGenericVanguard();

    static SimulationModelParams serializationTestParams();

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static std::string canonicalDeckPath(const std::string& caseName);

    /*!
     * \brief Returns the wall time required to set up the simulator before it was born.
     */
    double setupTime()
    { return setupTime_; }

    /*!
     * \brief Read a deck.
     * \param filename file to read
     */
    static void readDeck(const std::string& filename);

    /*!
     * \brief Set the simulation configuration objects.
     */
    void defineSimulationModel(SimulationModelParams&& params);

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

    std::unique_ptr<WellTestState> transferWTestState() {
        return std::move(this->wtestState_);
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
     * \brief Number of blocks in the Block-Jacobi preconditioner.
     */
    int numJacobiBlocks() const
    {
#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
        return numJacobiBlocks_;
#else
        return 0;
#endif // HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
    }

    /*!
     * \brief Parameter that decide if cells owned by rank are ordered before ghost cells.
     */
    bool ownersFirst() const
    { return ownersFirst_; }

#if HAVE_MPI
    bool addCorners() const
    { return addCorners_; }

    int numOverlap() const
    { return numOverlap_; }

    /*!
     * \brief Parameter deciding which partition method to use
     */
    Dune::PartitionMethod partitionMethod() const
    { return partitionMethod_; }

    /*!
     * \brief Parameter that decides if partitioning for parallel runs
     *        should be performed on a single process only.
     */
    bool serialPartitioning() const
    { return serialPartitioning_; }

    /*!
     * \brief Parameter that sets the imbalance tolarance, depending on the
     * chosen partition method
     */
    double imbalanceTol() const
    {
        if (zoltanImbalanceTolSet_) {
            OpmLog::info("The parameter --zoltan-imbalance-tol is deprecated "
                         "and has been renamed to --imbalance-tol, please "
                         "adjust your calls and scripts!");
            return zoltanImbalanceTol_;
        } else {
            return imbalanceTol_;
        }
    }

    const std::string& externalPartitionFile() const
    {
        return this->externalPartitionFile_;
    }
#endif // HAVE_MPI

    /*!
     * \brief Whether perforations of a well might be distributed.
     */
    bool enableDistributedWells() const
    { return enableDistributedWells_; }

    /*!
     * \brief Whether or not to emit result files that are compatible with
     * a commercial reservoir simulator.
     */
    bool enableEclOutput() const
    { return enableEclOutput_; }

    /*!
     * \brief Retrieve collection (a vector of pairs) of well names and
     * whether or not the corresponding well objects are perforated on the
     * current rank.
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
        assert(comm_ != nullptr);
        return *comm_;
    }

    // Private to avoid pulling schedule in header.
    // Static state is not serialized, only use for restart.
    template<class Serializer>
    void serializeOp(Serializer& serializer);

    // Only compares dynamic state.
    bool operator==(const FlowGenericVanguard& rhs) const;

protected:
    void updateOutputDir_(std::string outputDir,
                          bool enableEclCompatFile);

    void updateNOSIM_(std::string_view enableDryRun);

    bool drsdtconEnabled() const;

    std::unordered_map<std::size_t, const NumericalAquiferCell*> allAquiferCells() const;

    void init();

    template<class Scalar>
    static void registerParameters_();

    double setupTime_;

    // These variables may be owned by both Python and the simulator
    static std::unique_ptr<Parallel::Communication> comm_;

    std::string caseName_;
    std::string fileName_;
    Dune::EdgeWeightMethod edgeWeightsMethod_;

#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
    int numJacobiBlocks_{0};
#endif // HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA

    bool ownersFirst_;
#if HAVE_MPI
    bool addCorners_;
    int numOverlap_;

    Dune::PartitionMethod partitionMethod_;
    bool serialPartitioning_;
    double imbalanceTol_;

    bool zoltanImbalanceTolSet_;
    double zoltanImbalanceTol_;
    double zoltanPhgEdgeSizeThreshold_;
    std::string zoltanParams_;

    std::string metisParams_;

    std::string externalPartitionFile_{};
#endif // HAVE_MPI

    bool enableDistributedWells_;
    bool enableEclOutput_;
    bool allow_splitting_inactive_wells_;

    std::string ignoredKeywords_;
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
    std::shared_ptr<Python> python;
    // These variables may be owned by both Python and the simulator
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

#endif // OPM_FLOW_GENERIC_VANGUARD_HPP
