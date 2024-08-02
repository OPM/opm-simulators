// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil;
// c-basic-offset: 4 -*- vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or
  modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation,
  either version 2 of the License, or (at your option) any
  later version. OPM is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License
  for more details. You should have received a copy of the
  GNU General Public License along with OPM.  If not, see
  <http://www.gnu.org/licenses/>. Consult the COPYING file
  in the top-level source directory of this module for the
  precise wording of the license and the list of copyright
  holders.
*/
/*!
 * \file
 * \copydoc Opm::UnstructuredGridVanguard
 */
#ifndef EWOMS_UNSTRUCTURED_GRID_VANGUARD_HH
#define EWOMS_UNSTRUCTURED_GRID_VANGUARD_HH

#include <opm/models/io/basevanguard.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#ifdef HAVE_OPM_GRID
#include "opm/grid/UnstructuredGrid.h"
#endif

namespace Opm {
/*!
 * \brief Provides a simulator vanguard which creates a grid
 * by parsing an unstructured grid file
 */
template <class TypeTag>
class UnstructuredGridVanguard : public BaseVanguard<TypeTag> {
    using ParentType = BaseVanguard<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    using GridPointer = Dune::GridPtr< Grid >;

   public:
    /*!
     * \brief Register all run-time parameters for the
     * unstructured grid simulator vanguard.
     */
    static void registerParameters() {
        Parameters::registerParam<TypeTag, Parameters::GridGlobalRefinements>
            ("The number of global refinements of the grid "
             "executed after it was loaded");
        Parameters::registerParam<TypeTag, Parameters::GridFile>,
            ("The file name of the file to load");
    }

    /*!
     * \brief Load the grid from the file.
     */
    UnstructuredGridVanguard(Simulator& simulator) : ParentType(simulator){
#ifdef HAVE_OPM_GRID
        const std::string gridFileName = Parameters::get<TypeTag, Parameters::GridFile>();
        unsigned numRefinments = Parameters::get<TypeTag, Parameters::GridGlobalRefinements>();

        const char* c_str = gridFileName.c_str();

        UnstructuredGrid* grid = read_grid(c_str);
        if (grid == nullptr) {
            std::string msg =
                "RuntimeError: UnstructuredGridVanguard could not read grid file: " +
                gridFileName + ". Are you sure the filename is correct?";
            throw std::runtime_error(msg);
        }
        ugPtr_.reset(std::move( grid ));
        //GridPointer polygrid( new Grid(*ugPtr) );
        gridPtr_ = new Grid(*ugPtr_);//std::move(polygrid);
        if (numRefinments > 0) {
            gridPtr_->globalRefine(static_cast<int>(numRefinments));
        }
        this->finalizeInit_();
#endif
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid& grid() { return *gridPtr_; }

    /*!
     * \brief Return a constant reference to the grid
     * object.
     */
    const Grid& grid() const { return *gridPtr_; }

   private:
    GridPointer gridPtr_;
    typename Grid::UnstructuredGridPtr ugPtr_;
};

}  // namespace Opm

#endif
