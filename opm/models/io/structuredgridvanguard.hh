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
 * \copydoc Opm::StructuredGridVanguard
 */
#ifndef EWOMS_STRUCTURED_GRID_VANGUARD_HH
#define EWOMS_STRUCTURED_GRID_VANGUARD_HH

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <opm/models/io/basevanguard.hh>

#include <opm/models/utils/basicparameters.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#include <memory>
#include <sstream>

namespace Opm {

template <class TypeTag>
class StructuredGridVanguard;

} // namespace Opm

namespace Opm::Properties {

namespace TTag {

struct StructuredGridVanguard {};

} // namespace TTag

// GRIDDIM is only set by the finger problem
#ifndef GRIDDIM
static constexpr int dim = 2;
#else
static constexpr int dim = GRIDDIM;
#endif

// set the Grid and Vanguard properties
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::StructuredGridVanguard>
{ using type = Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::StructuredGridVanguard>
{ using type = Dune::YaspGrid< dim >; };
#endif

template<class TypeTag>
struct Vanguard<TypeTag, TTag::StructuredGridVanguard>
{ using type = Opm::StructuredGridVanguard<TypeTag>; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class StructuredGridVanguard : public BaseVanguard<TypeTag>
{
    using ParentType = BaseVanguard<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    using GridPointer = std::unique_ptr<Grid>;

    static constexpr int dim = Grid::dimension;

public:
    /*!
     * \brief Register all run-time parameters for the structured grid simulator vanguard.
     */
    static void registerParameters()
    {
        Parameters::Register<Parameters::GridGlobalRefinements>
            ("The number of global refinements of the grid "
             "executed after it was loaded");
        Parameters::Register<Parameters::DomainSizeX<Scalar>>
            ("The size of the domain in x direction");
        Parameters::Register<Parameters::CellsX>
            ("The number of intervalls in x direction");
        if constexpr (dim > 1) {
            Parameters::Register<Parameters::DomainSizeY<Scalar>>
                ("The size of the domain in y direction");
            Parameters::Register<Parameters::CellsY>
                ("The number of intervalls in y direction");
        }
        if constexpr (dim > 2) {
            Parameters::Register<Parameters::DomainSizeZ<Scalar>>
                ("The size of the domain in z direction");
            Parameters::Register<Parameters::CellsZ>
                ("The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the grid for the lens problem
     */
    explicit StructuredGridVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        Dune::FieldVector<int, dim> cellRes;

        using GridScalar = double;
        Dune::FieldVector<GridScalar, dim> upperRight;
        Dune::FieldVector<GridScalar, dim> lowerLeft( 0 );

        upperRight[0] = Parameters::Get<Parameters::DomainSizeX<Scalar>>();
        upperRight[1] = Parameters::Get<Parameters::DomainSizeY<Scalar>>();

        cellRes[0] = Parameters::Get<Parameters::CellsX>();
        cellRes[1] = Parameters::Get<Parameters::CellsY>();
        if constexpr (dim == 3) {
            upperRight[2] = Parameters::Get<Parameters::DomainSizeZ<Scalar>>();
            cellRes[2] = Parameters::Get<Parameters::CellsZ>();
        }

        std::stringstream dgffile;
        dgffile << "DGF" << '\n'
                << "INTERVAL" << '\n'
                << lowerLeft  << '\n'
                << upperRight << '\n'
                << cellRes    <<  '\n'
                << "#" << '\n'
                << "GridParameter" << '\n'
                << "overlap 1" << '\n'
                << "#" << '\n'
                << "Simplex" << '\n'
                << "#" << '\n';

        // use DGF parser to create a grid from interval block
        gridPtr_.reset(Dune::GridPtr<Grid>(dgffile).release());

        const int numRefinements = Parameters::Get<Parameters::GridGlobalRefinements>();
        gridPtr_->globalRefine(numRefinements);

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Return a constant reference to the grid object.
     */
    const Grid& grid() const
    { return *gridPtr_; }

private:
    GridPointer gridPtr_;
};

} // namespace Opm

#endif
