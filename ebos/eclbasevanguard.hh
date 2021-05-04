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
#ifndef EWOMS_ECL_BASE_VANGUARD_HH
#define EWOMS_ECL_BASE_VANGUARD_HH

#include <opm/models/io/basevanguard.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/parser/eclipse/Python/Python.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQConfig.hpp>
#include <opm/common/utility/TimeService.hpp>


#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <array>
#include <chrono>
#include <unordered_set>
#include <vector>

namespace Opm {
template <class TypeTag>
class EclBaseVanguard;
}

namespace Opm::Properties {

namespace TTag {
struct EclBaseVanguard {};
}

// declare the properties required by the for the ecl simulator vanguard
template<class TypeTag, class MyTypeTag>
struct EquilGrid {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EclDeckFileName {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableOpmRstFile {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EclStrictParsing {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SchedRestart {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EclOutputInterval {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct IgnoreKeywords {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EdgeWeightsMethod {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct OwnerCellsFirst {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SerialPartitioning {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct ZoltanImbalanceTol {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct AllowDistributedWells {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct IgnoreKeywords<TypeTag, TTag::EclBaseVanguard> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EclDeckFileName<TypeTag, TTag::EclBaseVanguard> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EclOutputInterval<TypeTag, TTag::EclBaseVanguard> {
    static constexpr int value = -1;
};
template<class TypeTag>
struct EnableOpmRstFile<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EclStrictParsing<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct SchedRestart<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EdgeWeightsMethod<TypeTag, TTag::EclBaseVanguard> {
    static constexpr int value = 1;
};
template<class TypeTag>
struct OwnerCellsFirst<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct SerialPartitioning<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct ZoltanImbalanceTol<TypeTag, TTag::EclBaseVanguard> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.1;
};

template<class TypeTag>
struct AllowDistributedWells<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};

template<class T1, class T2>
struct UseMultisegmentWell;

// Same as in BlackoilModelParametersEbos.hpp but for here.
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = true;
};
} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 */
template <class TypeTag>
class EclBaseVanguard : public BaseVanguard<TypeTag>
{
    using ParentType = BaseVanguard<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Vanguard>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

protected:
    static const int dimension = Grid::dimension;
    using Element = typename GridView::template Codim<0>::Entity;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;


public:
    /*!
     * \brief Register the common run-time parameters for all ECL simulator vanguards.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclDeckFileName,
                             "The name of the file which contains the ECL deck to be simulated");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclOutputInterval,
                             "The number of report steps that ought to be skipped between two writes of ECL results");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableOpmRstFile,
                             "Include OPM-specific keywords in the ECL restart file to enable restart of OPM simulators from these files");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, IgnoreKeywords,
                             "List of Eclipse keywords which should be ignored. As a ':' separated string.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclStrictParsing,
                             "Use strict mode for parsing - all errors are collected before the applicaton exists.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SchedRestart,
                             "When restarting: should we try to initialize wells and groups from historical SCHEDULE section.");
        EWOMS_REGISTER_PARAM(TypeTag, int, EdgeWeightsMethod,
                             "Choose edge-weighing strategy: 0=uniform, 1=trans, 2=log(trans).");
        EWOMS_REGISTER_PARAM(TypeTag, bool, OwnerCellsFirst,
                             "Order cells owned by rank before ghost/overlap cells.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SerialPartitioning,
                             "Perform partitioning for parallel runs on a single process.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, ZoltanImbalanceTol,
                             "Tolerable imbalance of the loadbalancing provided by Zoltan (default: 1.1).");
        EWOMS_REGISTER_PARAM(TypeTag, bool, AllowDistributedWells,
                             "Allow the perforations of a well to be distributed to interior of multiple processes");
        // register here for the use in the tests without BlackoildModelParametersEbos
        EWOMS_REGISTER_PARAM(TypeTag, bool, UseMultisegmentWell, "Use the well model for multi-segment wells instead of the one for single-segment wells");

    }

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static std::string canonicalDeckPath(const std::string& caseName)
    {
        const auto fileExists = [](const filesystem::path& f) -> bool
            {
                if (!filesystem::exists(f))
                    return false;

                if (filesystem::is_regular_file(f))
                    return true;

                return filesystem::is_symlink(f) && filesystem::is_regular_file(filesystem::read_symlink(f));
            };

        auto simcase = filesystem::path(caseName);
        if (fileExists(simcase))
            return simcase.string();

        for (const auto& ext : { std::string("data"), std::string("DATA") }) {
            if (fileExists(simcase.replace_extension(ext)))
                return simcase.string();
        }

        throw std::invalid_argument("Cannot find input case '"+caseName+"'");
    }

    /*!
     * \brief Creates an Opm::parseContext object assuming that the parameters are ready.
     */
    static std::unique_ptr<ParseContext> createParseContext(const std::string& ignoredKeywords,
                                                            bool eclStrictParsing)
    {
        typedef std::pair<std::string, InputError::Action> ParseModePair;
        typedef std::vector<ParseModePair> ParseModePairs;
        ParseModePairs tmp;
        tmp.emplace_back(ParseContext::PARSE_RANDOM_SLASH, InputError::IGNORE);
        tmp.emplace_back(ParseContext::PARSE_MISSING_DIMS_KEYWORD, InputError::WARN);
        tmp.emplace_back(ParseContext::SUMMARY_UNKNOWN_WELL, InputError::WARN);
        tmp.emplace_back(ParseContext::SUMMARY_UNKNOWN_GROUP, InputError::WARN);
        tmp.emplace_back(ParseContext::PARSE_EXTRA_RECORDS, InputError::WARN);

        auto parseContext = std::make_unique<ParseContext>(tmp);

        if (ignoredKeywords.size() > 0) {
            size_t pos;
            size_t offset = 0;
            while (true) {
                pos = ignoredKeywords.find(':', offset);
                if (pos == std::string::npos) {
                    parseContext->ignoreKeyword(ignoredKeywords.substr(offset));
                    break;
                }
                parseContext->ignoreKeyword(ignoredKeywords.substr(offset, pos - offset));
                offset = pos + 1;
            }
        }

        if (eclStrictParsing)
            parseContext->update(InputError::DELAYED_EXIT1);

        return parseContext;
    }

    /*!
     * \brief Set the wall time which was spend externally to set up the external data structures
     *
     * i.e., the objects specified via the other setExternal*() methods.
     */
    static void setExternalSetupTime(Scalar t)
    { externalSetupTime_ = t; }

    /*!
     * \brief Returns the wall time required to set up the simulator before it was born.
     */
    static Scalar externalSetupTime()
    { return externalSetupTime_; }

    /*!
     * \brief Set the Opm::ParseContext object which ought to be used for parsing the deck and creating the Opm::EclipseState object.
     */
    static void setExternalParseContext(std::unique_ptr<ParseContext> parseContext)
    { externalParseContext_ = std::move(parseContext); }

    /*!
     * \brief Set the Opm::ErrorGuard object which ought to be used for parsing the deck and creating the Opm::EclipseState object.
     */
    static void setExternalErrorGuard(std::unique_ptr<ErrorGuard> errorGuard)
    { externalErrorGuard_ = std::move(errorGuard); }

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
    static void setExternalDeck(std::unique_ptr<Deck> deck)
    { externalDeck_ = std::move(deck); externalDeckSet_ = true; }
    /*!
     * \brief Set the Opm::EclipseState object which ought to be used when the simulator
     *        vanguard is instantiated.
     */
    static void setExternalEclState(std::unique_ptr<EclipseState> eclState)
    { externalEclState_ = std::move(eclState); }

    /*!
     * \brief Create the grid for problem data files which use the ECL file format.
     *
     * This is the file format used by the commercial ECLiPSE simulator. Usually it uses
     * a cornerpoint description of the grid.
     */
    EclBaseVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        int myRank = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        edgeWeightsMethod_   = Dune::EdgeWeightMethod(EWOMS_GET_PARAM(TypeTag, int, EdgeWeightsMethod));
        ownersFirst_ = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
        serialPartitioning_ = EWOMS_GET_PARAM(TypeTag, bool, SerialPartitioning);
        zoltanImbalanceTol_ = EWOMS_GET_PARAM(TypeTag, Scalar, ZoltanImbalanceTol);
        enableDistributedWells_ = EWOMS_GET_PARAM(TypeTag, bool, AllowDistributedWells);

        // Make proper case name.
        {
            if (fileName == "")
                throw std::runtime_error("No input deck file has been specified as a command line argument,"
                                        " or via '--ecl-deck-file-name=CASE.DATA'");

            fileName = canonicalDeckPath(fileName);

            // compute the base name of the input file name
            const char directorySeparator = '/';
            long int i;
            for (i = fileName.size(); i >= 0; -- i)
                if (fileName[i] == directorySeparator)
                    break;
            std::string baseName = fileName.substr(i + 1, fileName.size());

            // remove the extension from the input file
            for (i = baseName.size(); i >= 0; -- i)
                if (baseName[i] == '.')
                    break;
            std::string rawCaseName;
            if (i < 0)
                rawCaseName = baseName;
            else
                rawCaseName = baseName.substr(0, i);

            // transform the result to ALL_UPPERCASE
            caseName_ = rawCaseName;
            std::transform(caseName_.begin(), caseName_.end(), caseName_.begin(), ::toupper);
        }

        std::unique_ptr<ErrorGuard> errorGuard = nullptr;

        // Check that we are in one of the known configurations for external variables
        // and move them to internal
        if (externalDeck_)
        {
            deck_ = std::move(externalDeck_);

            if (externalParseContext_ && externalErrorGuard_ )
            {
                parseContext_ = std::move(externalParseContext_);
                errorGuard = std::move(externalErrorGuard_);
            }
            else if(externalEclState_ && externalEclSchedule_ && externalEclSummaryConfig_)
            {
                eclState_ = std::move(externalEclState_);
                eclSchedule_ = std::move(externalEclSchedule_);
                eclSummaryConfig_ = std::move(externalEclSummaryConfig_);
            }
            else
            {
                OPM_THROW(std::logic_error, "Either parse context and error guard or ECL state, schedule, and summary config need to be"
                          << " set externally.");
            }
        }
        else if (externalParseContext_)
        {
            parseContext_ = std::move(externalParseContext_);
        }
        else
        {
            const std::string ignoredKeywords = EWOMS_GET_PARAM(TypeTag, std::string, IgnoreKeywords);
            bool eclStrictParsing = EWOMS_GET_PARAM(TypeTag, bool, EclStrictParsing);
            parseContext_ = createParseContext(ignoredKeywords, eclStrictParsing);
        }

        std::optional<int> outputInterval;
        int output_param = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (output_param >= 0)
            outputInterval = output_param;

        readDeck(myRank, fileName, deck_, eclState_, eclSchedule_,
                 eclSummaryConfig_, std::move(errorGuard), python,
                 std::move(parseContext_), /* initFromRestart = */ false,
                 /* checkDeck = */ enableExperiments, outputInterval);

        this->summaryState_ = std::make_unique<SummaryState>( TimeService::from_time_t(this->eclSchedule_->getStartTime() ));
        this->udqState_ = std::make_unique<UDQState>( this->eclSchedule_->getUDQConfig(0).params().undefinedValue() );
        this->actionState_ = std::make_unique<Action::State>() ;

        // Initialize parallelWells with all local wells
        const auto& schedule_wells = schedule().getWellsatEnd();
        parallelWells_.reserve(schedule_wells.size());

        for (const auto& well: schedule_wells)
        {
            parallelWells_.emplace_back(well.name(), true);
        }
        std::sort(parallelWells_.begin(), parallelWells_.end());

        // Check whether allowing distribute wells makes sense
        if (enableDistributedWells() )
        {
            int hasMsWell = false;

            if (EWOMS_GET_PARAM(TypeTag, bool, UseMultisegmentWell))
            {
                if (myRank == 0)
                {
                    const auto& wells = this->schedule().getWellsatEnd();
                    for ( const auto& well: wells)
                    {
                        hasMsWell = hasMsWell || well.isMultiSegment();
                    }
                }
            }
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
            const auto& comm = Dune::MPIHelper::getCommunication();
#else
            const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
#endif
            hasMsWell = comm.max(hasMsWell);

            if (hasMsWell)
            {
                if (myRank == 0)
                {
                    std::string message =
                        std::string("Option --allow-distributed-wells=true is only allowed if model\n")
                        + "only has only standard wells. You need to provide option \n"
                        + " with --enable-multisegement-wells=false to treat existing \n"
                        + "multisegment wells as standard wells.";
                    OpmLog::error(message);
                }
                comm.barrier();
                OPM_THROW(std::invalid_argument, "All wells need to be standard wells!");
            }
        }
    }

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
     * \brief Set the schedule object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSchedule(std::unique_ptr<Schedule> schedule)
    { externalEclSchedule_ = std::move(schedule); }

    /*!
     * \brief Return a reference to the object that determines which quantities ought to
     *        be put into the ECL summary output.
     */
    const SummaryConfig& summaryConfig() const
    { return *eclSummaryConfig_; }

    /*!
     * \brief Set the summary configuration object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSummaryConfig(std::unique_ptr<SummaryConfig> summaryConfig)
    { externalEclSummaryConfig_ = std::move(summaryConfig); }


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
    Scalar zoltanImbalanceTol() const
    { return zoltanImbalanceTol_; }

    /*!
     * \brief Whether perforations of a well might be distributed.
     */
    bool enableDistributedWells() const
    { return enableDistributedWells_; }

    /*!
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string& caseName() const
    { return caseName_; }

    const CartesianIndexMapper& cartesianMapper() const
    {  return asImp_().cartesianIndexMapper(); }

    /*!
     * \brief Returns the number of logically Cartesian cells in each direction
     */
    const std::array<int, dimension>& cartesianDimensions() const
    { return asImp_().cartesianIndexMapper().cartesianDimensions(); }

    /*!
     * \brief Returns the overall number of cells of the logically Cartesian grid
     */
    int cartesianSize() const
    { return asImp_().cartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the overall number of cells of the logically EquilCartesian grid
     */
    int equilCartesianSize() const
    { return asImp_().equilCartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the Cartesian cell id for identifaction with ECL data
     */
    unsigned cartesianIndex(unsigned compressedCellIdx) const
    { return asImp_().cartesianIndexMapper().cartesianIndex(compressedCellIdx); }

    /*!
     * \brief Return the index of the cells in the logical Cartesian grid
     */
    unsigned cartesianIndex(const std::array<int,dimension>& coords) const
    {
        unsigned cartIndex = coords[0];
        int factor = cartesianDimensions()[0];
        for (unsigned i = 1; i < dimension; ++i) {
            cartIndex += coords[i]*factor;
            factor *= cartesianDimensions()[i];
        }

        return cartIndex;
    }

    /*!
     * \brief Return compressed index from cartesian index
     *
     */
    int compressedIndex(int cartesianCellIdx) const
    {
        int index = cartesianToCompressed_[cartesianCellIdx];
        return index;
    }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void cartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().cartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }

    /*!
     * \brief Returns the Cartesian cell id given an element index for the grid used for equilibration
     */
    unsigned equilCartesianIndex(unsigned compressedEquilCellIdx) const
    { return asImp_().equilCartesianIndexMapper().cartesianIndex(compressedEquilCellIdx); }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell of the grid used for EQUIL.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void equilCartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().equilCartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }


    /*!
     * \brief Returns vector with name and whether the has local perforated cells
     *        for all wells.
     *
     * Will only have usable values for CpGrid.
     */
    const std::vector<std::pair<std::string,bool>>& parallelWells() const
    { return parallelWells_; }

    /*!
     * \brief Get the cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     */
    const std::vector<double>& cellCentroids() const
    {
        return centroids_;
    }

    /*!
     * \brief Returns the depth of a degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    Scalar cellCenterDepth(unsigned globalSpaceIdx) const
    {
        return cellCenterDepth_[globalSpaceIdx];
    }

    /*!
     * \brief Returns the thickness of a degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depths of the top surface
     * corners minus the average of the depths of the bottom surface corners
     * The cell thickness is computed only when needed.
     */
    Scalar cellThickness(unsigned globalSpaceIdx) const
    {
        assert(!cellThickness_.empty());
        return cellThickness_[globalSpaceIdx];
    }

    /*!
     * \brief Get the number of cells in the global leaf grid view.
     * \warn This is a collective operation that needs to be called
     * on all ranks.
     */
    std::size_t globalNumCells() const
    {
        const auto& grid = asImp_().grid();
        if (grid.comm().size() == 1)
        {
            return grid.leafGridView().size(0);
        }
        const auto& gridView = grid.leafGridView();
        constexpr int codim = 0;
        constexpr auto Part = Dune::Interior_Partition;
        auto local_cells = std::distance(gridView.template begin<codim, Part>(),
                                         gridView.template end<codim, Part>());
        return grid.comm().sum(local_cells);
    }

protected:
    void callImplementationInit()
    {
        asImp_().createGrids_();
        asImp_().filterConnections_();
        asImp_().updateOutputDir_();
        asImp_().finalizeInit_();
    }
    void updateCartesianToCompressedMapping_()
    {
        size_t num_cells = asImp_().grid().leafGridView().size(0);
        cartesianToCompressed_.resize(cartesianSize(), -1);
        for (unsigned i = 0; i < num_cells; ++i) {
            unsigned cartesianCellIdx = cartesianIndex(i);
            cartesianToCompressed_[cartesianCellIdx] = i;
        }
    }

    void updateCellDepths_()
    {
        int numCells = this->gridView().size(/*codim=*/0);
        cellCenterDepth_.resize(numCells);

        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
        auto elemIt = this->gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = this->gridView().template end</*codim=*/0>();

        const auto num_aqu_cells = this->eclState_->aquifer().numericalAquifers().allAquiferCells();

        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);
            cellCenterDepth_[elemIdx] = cellCenterDepth(element);

            if (!num_aqu_cells.empty()) {
                const unsigned int global_index = cartesianIndex(elemIdx);
                const auto search = num_aqu_cells.find(global_index);
                if (search != num_aqu_cells.end()) {
                    // updating the cell depth using aquifer cell depth
                    cellCenterDepth_[elemIdx] = search->second->depth;
                }
            }
        }
    }
    void updateCellThickness_()
    {
        bool drsdtcon = false;
        auto schIt = this->schedule().begin();
        const auto& schEndIt = this->schedule().end();
        for(; schIt != schEndIt; ++schIt) {
            const auto& oilVaporizationControl = schIt->oilvap();
            if(oilVaporizationControl.getType() == OilVaporizationProperties::OilVaporization::DRSDTCON) {
                drsdtcon = true;
                break;
            }
        }
        if (!drsdtcon)
            return;

        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());

        int numElements = this->gridView().size(/*codim=*/0);
        cellThickness_.resize(numElements);

        auto elemIt = this->gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = this->gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);
            cellThickness_[elemIdx] = asImp_().computeCellThickness(element);
        }
    }
    Scalar computeCellThickness(const Element&) const {
        OPM_THROW(std::runtime_error, "cellThickness not implemented for this grid!");
    }


private:
    void updateOutputDir_()
    {
        // update the location for output
        std::string outputDir = EWOMS_GET_PARAM(TypeTag, std::string, OutputDir);
        auto& ioConfig = eclState_->getIOConfig();
        if (outputDir == "")
            // If no output directory parameter is specified, use the output directory
            // which Opm::IOConfig thinks that should be used. Normally this is the
            // directory in which the input files are located.
            outputDir = ioConfig.getOutputDir();

        // ensure that the output directory exists and that it is a directory
        if (!filesystem::is_directory(outputDir)) {
            try {
                filesystem::create_directories(outputDir);
            }
            catch (...) {
                 throw std::runtime_error("Creation of output directory '"+outputDir+"' failed\n");
            }
        }

        // specify the directory output. This is not a very nice mechanism because
        // the eclState is supposed to be immutable here, IMO.
        ioConfig.setOutputDir(outputDir);

        ioConfig.setEclCompatibleRST(!EWOMS_GET_PARAM(TypeTag, bool, EnableOpmRstFile));
    }


    // computed from averaging cell corner depths
    Scalar cellCenterDepth(const Element& element) const
    {
        typedef typename Element::Geometry Geometry;
        static constexpr int zCoord = Element::dimension - 1;
        Scalar zz = 0.0;

        const Geometry& geometry = element.geometry();
        const int corners = geometry.corners();
        for (int i=0; i < corners; ++i)
            zz += geometry.corner(i)[zCoord];

        return zz/Scalar(corners);
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::string caseName_;

    static double externalSetupTime_;

    static std::unique_ptr<ParseContext> externalParseContext_;
    static std::unique_ptr<ErrorGuard> externalErrorGuard_;
    static std::unique_ptr<Deck> externalDeck_;
    static bool externalDeckSet_;
    static std::unique_ptr<EclipseState> externalEclState_;
    static std::unique_ptr<Schedule> externalEclSchedule_;
    static std::unique_ptr<SummaryConfig> externalEclSummaryConfig_;

    std::unique_ptr<SummaryState> summaryState_;
    std::unique_ptr<Action::State> actionState_;
    std::unique_ptr<UDQState> udqState_;

    // these attributes point  either to the internal  or to the external version of the
    // parser objects.
    std::unique_ptr<ParseContext> parseContext_;
    std::unique_ptr<ErrorGuard> errorGuard_;
    std::unique_ptr<Deck> deck_;
    std::unique_ptr<EclipseState> eclState_;
    std::unique_ptr<Schedule> eclSchedule_;
    std::unique_ptr<SummaryConfig> eclSummaryConfig_;
    std::shared_ptr<Python> python = std::make_shared<Python>();

    Dune::EdgeWeightMethod edgeWeightsMethod_;
    bool ownersFirst_;
    bool serialPartitioning_;
    Scalar zoltanImbalanceTol_;
    bool enableDistributedWells_;

protected:
    /*! \brief The cell centroids after loadbalance was called.
     * Empty otherwise. Used by EclTransmissibilty.
     */
    std::vector<double> centroids_;

    /*! \brief Mapping between cartesian and compressed cells.
     *  It is initialized the first time it is called
     */
    std::vector<int> cartesianToCompressed_;

    /*! \brief Cell center depths
     */
    std::vector<Scalar> cellCenterDepth_;

    /*! \brief Cell thichness
     */
    std::vector<Scalar> cellThickness_;

    /*! \brief information about wells in parallel
     *
     * For each well in the model there is an entry with its name
     * and a boolean indicating whether it perforates local cells.
     */
    std::vector<std::pair<std::string,bool>> parallelWells_;


};

template <class TypeTag>
double EclBaseVanguard<TypeTag>::externalSetupTime_ = 0.0;

template <class TypeTag>
std::unique_ptr<ParseContext> EclBaseVanguard<TypeTag>::externalParseContext_ = nullptr;

template <class TypeTag>
std::unique_ptr<ErrorGuard> EclBaseVanguard<TypeTag>::externalErrorGuard_ = nullptr;

template <class TypeTag>
std::unique_ptr<Deck> EclBaseVanguard<TypeTag>::externalDeck_ = nullptr;

template <class TypeTag>
bool EclBaseVanguard<TypeTag>::externalDeckSet_ = false;

template <class TypeTag>
std::unique_ptr<EclipseState> EclBaseVanguard<TypeTag>::externalEclState_;

template <class TypeTag>
std::unique_ptr<Schedule> EclBaseVanguard<TypeTag>::externalEclSchedule_ = nullptr;

template <class TypeTag>
std::unique_ptr<SummaryConfig> EclBaseVanguard<TypeTag>::externalEclSummaryConfig_ = nullptr;

} // namespace Opm

#endif
