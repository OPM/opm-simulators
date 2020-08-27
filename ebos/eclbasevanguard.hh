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
#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>

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

NEW_TYPE_TAG(EclBaseVanguard);

// declare the properties required by the for the ecl simulator vanguard
NEW_PROP_TAG(EquilGrid);
NEW_PROP_TAG(EclDeckFileName);
NEW_PROP_TAG(EnableOpmRstFile);
NEW_PROP_TAG(EclStrictParsing);
NEW_PROP_TAG(SchedRestart);
NEW_PROP_TAG(EclOutputInterval);
NEW_PROP_TAG(IgnoreKeywords);
NEW_PROP_TAG(EdgeWeightsMethod);
NEW_PROP_TAG(OwnerCellsFirst);

SET_STRING_PROP(EclBaseVanguard, IgnoreKeywords, "");
SET_STRING_PROP(EclBaseVanguard, EclDeckFileName, "");
SET_INT_PROP(EclBaseVanguard, EclOutputInterval, -1); // use the deck-provided value
SET_BOOL_PROP(EclBaseVanguard, EnableOpmRstFile, false);
SET_BOOL_PROP(EclBaseVanguard, EclStrictParsing, false);
SET_BOOL_PROP(EclBaseVanguard, SchedRestart, false);
SET_INT_PROP(EclBaseVanguard, EdgeWeightsMethod, 1);
SET_BOOL_PROP(EclBaseVanguard, OwnerCellsFirst, true);

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

    enum { enableExperiments = GET_PROP_VALUE(TypeTag, EnableExperiments) };

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

protected:
    static const int dimension = Grid::dimension;

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
    }

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static Opm::filesystem::path canonicalDeckPath(const std::string& caseName)
    {
        const auto fileExists = [](const Opm::filesystem::path& f) -> bool
            {
                if (!Opm::filesystem::exists(f))
                    return false;

                if (Opm::filesystem::is_regular_file(f))
                    return true;

                return Opm::filesystem::is_symlink(f) && Opm::filesystem::is_regular_file(Opm::filesystem::read_symlink(f));
            };

        auto simcase = Opm::filesystem::path(caseName);
        if (fileExists(simcase))
            return simcase;

        for (const auto& ext : { std::string("data"), std::string("DATA") }) {
            if (fileExists(simcase.replace_extension(ext)))
                return simcase;
        }

        throw std::invalid_argument("Cannot find input case '"+caseName+"'");
    }

    /*!
     * \brief Creates an Opm::parseContext object assuming that the parameters are ready.
     */
    static std::unique_ptr<Opm::ParseContext> createParseContext()
    {
        typedef std::pair<std::string, Opm::InputError::Action> ParseModePair;
        typedef std::vector<ParseModePair> ParseModePairs;
        ParseModePairs tmp;
        tmp.emplace_back(Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE);
        tmp.emplace_back(Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN);
        tmp.emplace_back(Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN);
        tmp.emplace_back(Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN);
        tmp.emplace_back(Opm::ParseContext::PARSE_EXTRA_RECORDS, Opm::InputError::WARN);

        std::unique_ptr<Opm::ParseContext> parseContext(new Opm::ParseContext(tmp));

        const std::string ignoredKeywords = EWOMS_GET_PARAM(TypeTag, std::string, IgnoreKeywords);
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

        if (EWOMS_GET_PARAM(TypeTag, bool, EclStrictParsing))
            parseContext->update(Opm::InputError::DELAYED_EXIT1);

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
    static void setExternalParseContext(std::unique_ptr<Opm::ParseContext> parseContext)
    { externalParseContext_ = std::move(parseContext); }

    /*!
     * \brief Set the Opm::ErrorGuard object which ought to be used for parsing the deck and creating the Opm::EclipseState object.
     */
    static void setExternalErrorGuard(std::unique_ptr<Opm::ErrorGuard> errorGuard)
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
    static void setExternalDeck(std::unique_ptr<Opm::Deck> deck)
    { externalDeck_ = std::move(deck); externalDeckSet_ = true; }
    /*!
     * \brief Set the Opm::EclipseState object which ought to be used when the simulator
     *        vanguard is instantiated.
     */
    static void setExternalEclState(std::unique_ptr<Opm::EclipseState> eclState)
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

        // Make proper case name.
        {
            if (fileName == "")
                throw std::runtime_error("No input deck file has been specified as a command line argument,"
                                        " or via '--ecl-deck-file-name=CASE.DATA'");

            fileName = canonicalDeckPath(fileName).string();

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
            parseContext_ = createParseContext();
        }

        readDeck(myRank, fileName, deck_, eclState_, eclSchedule_,
                 eclSummaryConfig_, std::move(errorGuard), python,
                 std::move(parseContext_), /* initFromRestart = */ false,
                 /* checkDeck = */ enableExperiments);

        this->summaryState_ = std::make_unique<Opm::SummaryState>( std::chrono::system_clock::from_time_t(this->eclSchedule_->getStartTime() ));
        this->actionState_ = std::make_unique<Opm::Action::State>() ;

        // Possibly override IOConfig setting for how often RESTART files should get
        // written to disk (every N report step)
        int outputInterval = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (outputInterval >= 0)
            schedule().restart().overrideRestartWriteInterval(outputInterval);
    }

    /*!
     * \brief Return a reference to the parsed ECL deck.
     */
    const Opm::Deck& deck() const
    { return *deck_; }

    Opm::Deck& deck()
    { return *deck_; }

    /*!
     * \brief Return a reference to the internalized ECL deck.
     */
    const Opm::EclipseState& eclState() const
    { return *eclState_; }

    Opm::EclipseState& eclState()
    { return *eclState_; }

    /*!
     * \brief Return a reference to the object that managages the ECL schedule.
     */
    const Opm::Schedule& schedule() const
    { return *eclSchedule_; }

    Opm::Schedule& schedule()
    { return *eclSchedule_; }

    /*!
     * \brief Set the schedule object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSchedule(std::unique_ptr<Opm::Schedule> schedule)
    { externalEclSchedule_ = std::move(schedule); }

    /*!
     * \brief Return a reference to the object that determines which quantities ought to
     *        be put into the ECL summary output.
     */
    const Opm::SummaryConfig& summaryConfig() const
    { return *eclSummaryConfig_; }

    /*!
     * \brief Set the summary configuration object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSummaryConfig(std::unique_ptr<Opm::SummaryConfig> summaryConfig)
    { externalEclSummaryConfig_ = std::move(summaryConfig); }


    /*!
    * \brief Returns the summary state
    *
    * The summary state is a small container object for
    * computed, ready to use summary values. The values will typically be used by
    * the UDQ, WTEST and ACTIONX calculations.
    */
    Opm::SummaryState& summaryState()
    { return *summaryState_; }

    const Opm::SummaryState& summaryState() const
    { return *summaryState_; }


    /*!
     * \brief Returns the action state
     *
     * The ActionState keeps track of how many times the various actions have run.
     */
    Opm::Action::State& actionState()
    { return *actionState_; }

    const Opm::Action::State& actionState() const
    { return *actionState_; }

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
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string& caseName() const
    { return caseName_; }

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
     * \brief Return the names of the wells which do not penetrate any cells on the local
     *        process.
     *
     * This is a kludge around the fact that for distributed grids, not all wells are
     * seen by all proccesses.
     */
    std::unordered_set<std::string> defunctWellNames() const
    { return std::unordered_set<std::string>(); }

    /*!
     * \brief Get the cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     */
    const std::vector<double>& cellCentroids() const
    {
        return centroids_;
    }
protected:
    void callImplementationInit()
    {
        asImp_().createGrids_();
        asImp_().filterConnections_();
        asImp_().updateOutputDir_();
        asImp_().finalizeInit_();

        if (enableExperiments) {
            if (asImp_().grid().size(0)) //grid not loadbalanced yet for ebos!
            {
                Opm::RelpermDiagnostics relpermDiagnostics;
                relpermDiagnostics.diagnosis(*eclState_, asImp_().grid());
            }
        }
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
        if (!Opm::filesystem::is_directory(outputDir)) {
            try {
                Opm::filesystem::create_directories(outputDir);
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

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::string caseName_;

    static Scalar externalSetupTime_;

    static std::unique_ptr<Opm::ParseContext> externalParseContext_;
    static std::unique_ptr<Opm::ErrorGuard> externalErrorGuard_;
    static std::unique_ptr<Opm::Deck> externalDeck_;
    static bool externalDeckSet_;
    static std::unique_ptr<Opm::EclipseState> externalEclState_;
    static std::unique_ptr<Opm::Schedule> externalEclSchedule_;
    static std::unique_ptr<Opm::SummaryConfig> externalEclSummaryConfig_;

    std::unique_ptr<Opm::SummaryState> summaryState_;
    std::unique_ptr<Opm::Action::State> actionState_;

    // these attributes point  either to the internal  or to the external version of the
    // parser objects.
    std::unique_ptr<Opm::ParseContext> parseContext_;
    std::unique_ptr<Opm::ErrorGuard> errorGuard_;
    std::unique_ptr<Opm::Deck> deck_;
    std::unique_ptr<Opm::EclipseState> eclState_;
    std::unique_ptr<Opm::Schedule> eclSchedule_;
    std::unique_ptr<Opm::SummaryConfig> eclSummaryConfig_;
    std::shared_ptr<Opm::Python> python = std::make_shared<Opm::Python>();

    Dune::EdgeWeightMethod edgeWeightsMethod_;
    bool ownersFirst_;

protected:
    /*! \brief The cell centroids after loadbalance was called.
     * Empty otherwise. Used by EclTransmissibilty.
     */
    std::vector<double> centroids_;
};

template <class TypeTag>
typename EclBaseVanguard<TypeTag>::Scalar EclBaseVanguard<TypeTag>::externalSetupTime_ = 0.0;

template <class TypeTag>
std::unique_ptr<Opm::ParseContext> EclBaseVanguard<TypeTag>::externalParseContext_ = nullptr;

template <class TypeTag>
std::unique_ptr<Opm::ErrorGuard> EclBaseVanguard<TypeTag>::externalErrorGuard_ = nullptr;

template <class TypeTag>
std::unique_ptr<Opm::Deck> EclBaseVanguard<TypeTag>::externalDeck_ = nullptr;

template <class TypeTag>
bool EclBaseVanguard<TypeTag>::externalDeckSet_ = false;

template <class TypeTag>
std::unique_ptr<Opm::EclipseState> EclBaseVanguard<TypeTag>::externalEclState_;

template <class TypeTag>
std::unique_ptr<Opm::Schedule> EclBaseVanguard<TypeTag>::externalEclSchedule_ = nullptr;

template <class TypeTag>
std::unique_ptr<Opm::SummaryConfig> EclBaseVanguard<TypeTag>::externalEclSummaryConfig_ = nullptr;

} // namespace Opm

#endif
