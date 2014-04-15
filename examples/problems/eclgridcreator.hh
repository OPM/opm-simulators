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
 * \copydoc Ewoms::EclGridCreator
 */
#ifndef EWOMS_ECL_GRID_CREATOR_HH
#define EWOMS_ECL_GRID_CREATOR_HH

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/basegridcreator.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/CpGrid.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Ewoms {
template <class TypeTag>
class EclProblem;

template <class TypeTag>
class EclGridCreator;
} // namespace Ewoms

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(EclGridCreator);

// declare the properties required by the for the ecl grid creator
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclipseDeckFileName);

SET_STRING_PROP(EclGridCreator, EclipseDeckFileName, "grids/ecl.DATA");

// set the Grid and GridCreator properties
SET_TYPE_PROP(EclGridCreator, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclGridCreator, GridCreator, Ewoms::EclGridCreator<TypeTag>);
}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the ecl problem.
 */
template <class TypeTag>
class EclGridCreator : public BaseGridCreator<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef std::shared_ptr<Grid> GridPointer;

private:
    static const int dim = Grid::dimension;

public:
    /*!
     * \brief Register all run-time parameters for the grid creator.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclipseDeckFileName,
                             "The name of the file which contains the Eclipse deck to be simulated");
    }

    /*!
     * \brief Create the grid for the ecl problem
     */
    static void makeGrid()
    {
        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclipseDeckFileName);

        Opm::ParserPtr parser(new Opm::Parser());
        deck_ = parser->parseFile(fileName);

        std::shared_ptr<Opm::RUNSPECSection> runspecSection(new Opm::RUNSPECSection(deck_) );
        std::shared_ptr<Opm::GRIDSection> gridSection(new Opm::GRIDSection(deck_) );
        eclipseGrid_.reset(new Opm::EclipseGrid(runspecSection, gridSection));

        grid_ = GridPointer(new Grid());
        grid_->processEclipseFormat(deck_,
                                    /*zTolerance=*/0,
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false);
    }

    /*!
     * \brief Return a reference to the grid.
     */
    static GridPointer gridPointer()
    { return grid_; }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); }

    /*!
     * \brief Destroy the grid
     *
     * This is required to guarantee that the grid is deleted before
     * MPI_Comm_free is called.
     */
    static void deleteGrid()
    {
        grid_.reset();
        deck_.reset();
    }

    /*!
     * \brief Return a pointer to the parsed Eclipse deck
     */
    static Opm::DeckConstPtr deck()
    { return deck_; }

    /*!
     * \brief Return a pointer to the EclipseGrid object
     *
     * The EclipseGrid class is used to internalize the cornerpoint
     * grid representation and, amongst others, can be used to write
     * EGRID files (which tends to be difficult with a plain
     * Dune::CpGrid)
     */
    static Opm::EclipseGridConstPtr eclipseGrid()
    { return eclipseGrid_; }

private:
    static GridPointer grid_;
    static Opm::DeckConstPtr deck_;
    static Opm::EclipseGridConstPtr eclipseGrid_;
};

template <class TypeTag>
typename EclGridCreator<TypeTag>::GridPointer EclGridCreator<TypeTag>::grid_;

template <class TypeTag>
Opm::DeckConstPtr EclGridCreator<TypeTag>::deck_;

template <class TypeTag>
Opm::EclipseGridConstPtr EclGridCreator<TypeTag>::eclipseGrid_;

} // namespace Ewoms

#endif
