/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Ewoms::EclGridManager
 */
#ifndef EWOMS_ECL_GRID_MANAGER_HH
#define EWOMS_ECL_GRID_MANAGER_HH

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/basegridmanager.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/CpGrid.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Ewoms {
template <class TypeTag>
class EclProblem;

template <class TypeTag>
class EclGridManager;
} // namespace Ewoms

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(EclGridManager);

// declare the properties required by the for the ecl grid manager
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclipseDeckFileName);

SET_STRING_PROP(EclGridManager, EclipseDeckFileName, "data/ecl.DATA");

// set the Grid and GridManager properties
SET_TYPE_PROP(EclGridManager, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclGridManager, GridManager, Ewoms::EclGridManager<TypeTag>);
}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the ecl problem.
 */
template <class TypeTag>
class EclGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef std::unique_ptr<Grid> GridPointer;

public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclipseDeckFileName,
                             "The name of the file which contains the Eclipse deck to be simulated");
    }

    /*!
     * \brief Create the grid for the ecl problem
     */
    /*!
     * \brief Create the grid for the lens problem
     */
    EclGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclipseDeckFileName);
        boost::filesystem::path deckPath(fileName);
        caseName_ = boost::to_upper_copy(deckPath.stem().string());

        Opm::ParserPtr parser(new Opm::Parser());
        deck_ = parser->parseFile(deckPath.string());

        eclipseState_.reset(new Opm::EclipseState(deck_));

        grid_ = GridPointer(new Grid());
        grid_->processEclipseFormat(deck_,
                                    /*zTolerance=*/0,
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false);

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid.
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the grid.
     */
    const Grid& grid() const
    { return *grid_; }

    /*!
     * \brief Return a pointer to the parsed Eclipse deck
     */
    Opm::DeckConstPtr deck() const
    { return deck_; }

    /*!
     * \brief Return a pointer to the internalized Eclipse deck
     */
    Opm::EclipseStateConstPtr eclipseState() const
    { return eclipseState_; }

    /*!
     * \brief Return a pointer to the internalized schedule of the
     *        Eclipse deck
     */
    Opm::ScheduleConstPtr schedule() const
    { return eclipseState_->getSchedule(); }

    /*!
     * \brief Return a pointer to the EclipseGrid object
     *
     * The EclipseGrid class is used to internalize the cornerpoint
     * grid representation and, amongst others, can be used to write
     * EGRID files (which tends to be difficult with a plain
     * Dune::CpGrid)
     */
    Opm::EclipseGridConstPtr eclipseGrid() const
    { return eclipseState_->getEclipseGrid(); }

    /*!
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string &caseName() const
    { return caseName_; }

private:
    std::string caseName_;
    GridPointer grid_;
    Opm::DeckConstPtr deck_;
    Opm::EclipseStateConstPtr eclipseState_;
};

} // namespace Ewoms

#endif
