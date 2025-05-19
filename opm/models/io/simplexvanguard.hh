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
 * \copydoc Opm::SimplexGridVanguard
 */
#ifndef EWOMS_SIMPLEX_GRID_VANGUARD_HH
#define EWOMS_SIMPLEX_GRID_VANGUARD_HH

#include <dune/common/fvector.hh>

#include <dune/grid/utility/structuredgridfactory.hh>

#include <opm/models/io/basevanguard.hh>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <array>
#include <memory>

namespace Opm {

/*!
 * \brief Provides a simulator vanguard which a creates regular grid made of simplices.
 */
template <class TypeTag>
class SimplexGridVanguard
{
    using ParentType = BaseVanguard<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    using GridPointer = std::unique_ptr<Grid>;
    using CoordScalar = typename Grid::ctype;
    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld,
    };
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
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
        if (dimWorld > 1) {
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
     * \brief Create the Grid
     */
    explicit SimplexGridVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        std::array<unsigned, dim> cellRes{};
        GlobalPosition upperRight;
        GlobalPosition lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = Parameters::Get<Parameters::DomainSizeX<Scalar>>();
        cellRes[0] = Parameters::Get<Parameters::CellsX>();
        if constexpr (dim > 1) {
            lowerLeft[1] = 0.0;
            upperRight[1] = Parameters::Get<Parameters::DomainSizeY<Scalar>>();
            cellRes[1] = Parameters::Get<Parameters::CellsY>();
        }
        if constexpr (dim > 2) {
            lowerLeft[2] = 0.0;
            upperRight[2] = Parameters::Get<Parameters::DomainSizeZ<Scalar>>();
            cellRes[2] = Parameters::Get<Parameters::CellsZ>();
        }

        simplexGrid_ = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft,
                                                                            upperRight,
                                                                            cellRes);

        const unsigned numRefinments = Parameters::Get<Parameters::GridGlobalRefinements>();
        simplexGrid_->globalRefine(numRefinments);

        this->finalizeInit_();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return simplexGrid_; }

    /*!
     * \brief Returns a reference to the grid.
     */
    const Grid& grid() const
    { return *simplexGrid_; }

private:
    GridPointer simplexGrid_;
};

} // namespace Opm

#endif
