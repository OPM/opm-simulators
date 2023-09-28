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
 * \copydoc Opm::EclPolyhedralGridVanguard
 */
#ifndef EWOMS_ECL_POLYHEDRAL_GRID_VANGUARD_HH
#define EWOMS_ECL_POLYHEDRAL_GRID_VANGUARD_HH

#include <ebos/eclbasevanguard.hh>
#include <ebos/ecltransmissibility.hh>

#include <opm/grid/polyhedralgrid.hh>

#include <opm/models/common/multiphasebaseproperties.hh>

#include <array>
#include <functional>
#include <string>
#include <tuple>
#include <unordered_set>

namespace Opm {
template <class TypeTag>
class EclPolyhedralGridVanguard;
}

namespace Opm::Properties {

namespace TTag {
struct EclPolyhedralGridVanguard {
    using InheritsFrom = std::tuple<EclBaseVanguard>;
};
}

// declare the properties
template<class TypeTag>
struct Vanguard<TypeTag, TTag::EclPolyhedralGridVanguard> {
    using type = Opm::EclPolyhedralGridVanguard<TypeTag>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::EclPolyhedralGridVanguard> {
    using type = Dune::PolyhedralGrid<3, 3>;
};
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::EclPolyhedralGridVanguard> {
    using type = GetPropType<TypeTag, Properties::Grid>;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::PolyhedralGrid as the simulation grid.
 */
template <class TypeTag>
class EclPolyhedralGridVanguard : public EclBaseVanguard<TypeTag>
{
    friend class EclBaseVanguard<TypeTag>;
    using ParentType = EclBaseVanguard<TypeTag>;

    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using EquilCartesianIndexMapper = Dune::CartesianIndexMapper<EquilGrid>;
    static constexpr int dimension = Grid::dimension;
    static constexpr int dimensionworld = Grid::dimensionworld;

private:
    using GridPointer = Grid*;
    using EquilGridPointer = EquilGrid*;

public:
    using TransmissibilityType = EclTransmissibility<Grid, GridView, ElementMapper,
                                                     CartesianIndexMapper, Scalar>;

    EclPolyhedralGridVanguard(Simulator& simulator)
        : EclBaseVanguard<TypeTag>(simulator),
          simulator_( simulator )
    {
        this->callImplementationInit();
        // add a copy in standard vector format to fullfill new interface
        const int* globalcellorg = this->grid().globalCell();
        int num_cells = this->gridView().size(0);
        globalcell_.resize(num_cells);
        for(int i=0; i < num_cells; ++i){
            globalcell_[i] = globalcellorg[i];
        }
    }

    ~EclPolyhedralGridVanguard()
    {
    }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    const Grid& grid() const
    { return *grid_; }

    /*!
     * \brief Returns a refefence to the grid which should be used by the EQUIL
     *        initialization code.
     *
     * The EQUIL keyword is used to specify the initial condition of the reservoir in
     * hydrostatic equilibrium. Since the code which does this is not accepting arbitrary
     * DUNE grids (the code is part of the opm-core module), this is not necessarily the
     * same as the grid which is used for the actual simulation.
     */
    const EquilGrid& equilGrid() const
    { return *grid_; }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    { /* do nothing: The EQUIL grid is the simulation grid! */ }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    { /* do nothing: PolyhedralGrid is not parallel! */
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     *
     * Since PolyhedralGrid is not parallel, that's always the same as
     * cartesianIndexMapper().
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    const std::vector<int>& globalCell()
    {
        return globalcell_;
    }

    unsigned int gridEquilIdxToGridIdx(unsigned int elemIndex) const {
         return elemIndex;
    }

    unsigned int gridIdxToEquilGridIdx(unsigned int elemIndex) const {
        return elemIndex;
    }

    /*!
     * \brief Free the memory occupied by the global transmissibility object.
     *
     * After writing the initial solution, this array should not be necessary anymore.
     */
    void releaseGlobalTransmissibilities()
    {
    }

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const TransmissibilityType& globalTransmissibility() const
    {
        return simulator_.problem().eclTransmissibilities();
    }

    /*!
     * \brief Get function to query cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     * It is a function return the centroid for the given element
     * index.
     */
    std::function<std::array<double,EclBaseVanguard<TypeTag>::dimensionworld>(int)>
    cellCentroids() const
    {
        return this->cellCentroids_(this->cartesianIndexMapper(), false);
    }

    std::vector<int> cellPartition() const
    {
        // not required for this type of grid yet (only from bdaBridge??)
        return {};
    }
protected:
    void createGrids_()
    {
        grid_ = std::make_unique<Grid>(this->eclState().getInputGrid(), this->eclState().fieldProps().porv(true));
        cartesianIndexMapper_ = std::make_unique<CartesianIndexMapper>(*grid_);
        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
    }

    void filterConnections_()
    {
        // not handling the removal of completions for this type of grid yet.
    }

    Simulator& simulator_;

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    //CartesianIndexMapperPointer cartesianIndexMapper_;

    std::unordered_set<std::string> defunctWellNames_;
    std::vector<int> globalcell_;
};

} // namespace Opm

#endif
