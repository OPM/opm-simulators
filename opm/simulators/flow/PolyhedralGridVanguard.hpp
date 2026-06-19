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
 * \copydoc Opm::PolyhedralGridVanguard
 */
#ifndef OPM_POLYHEDRAL_GRID_VANGUARD_HPP
#define OPM_POLYHEDRAL_GRID_VANGUARD_HPP

#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/polyhedralgrid/levelcartesianindexmapper.hh>

#include <opm/models/common/multiphasebaseproperties.hh>

#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <fmt/format.h>

#include <array>
#include <functional>
#include <string>
#include <tuple>
#include <unordered_set>

namespace Opm {
template <class TypeTag>
class PolyhedralGridVanguard;
}

namespace Opm::Properties {

namespace TTag {
struct PolyhedralGridVanguard {
    using InheritsFrom = std::tuple<FlowBaseVanguard>;
};
}

// declare the properties
template<class TypeTag>
struct Vanguard<TypeTag, TTag::PolyhedralGridVanguard> {
    using type = Opm::PolyhedralGridVanguard<TypeTag>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::PolyhedralGridVanguard> {
    using type = Dune::PolyhedralGrid<3, 3>;
};
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::PolyhedralGridVanguard> {
    using type = GetPropType<TypeTag, Properties::Grid>;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::PolyhedralGrid as the simulation grid.
 */
template <class TypeTag>
class PolyhedralGridVanguard : public FlowBaseVanguard<TypeTag>
{
    friend class FlowBaseVanguard<TypeTag>;
    using ParentType = FlowBaseVanguard<TypeTag>;

    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using LevelCartesianIndexMapper = Opm::LevelCartesianIndexMapper<Grid>;
    using EquilCartesianIndexMapper = Dune::CartesianIndexMapper<EquilGrid>;
    static constexpr int dimension = Grid::dimension;
    static constexpr int dimensionworld = Grid::dimensionworld;

private:
    using GridPointer = Grid*;
    using EquilGridPointer = EquilGrid*;

public:
    using TransmissibilityType = Transmissibility<Grid, GridView, ElementMapper,
                                                  CartesianIndexMapper, Scalar>;

    explicit PolyhedralGridVanguard(Simulator& simulator)
        : FlowBaseVanguard<TypeTag>(simulator)
        , simulator_(simulator)
    {
        this->callImplementationInit();
        // add a copy in standard vector format to fullfill new interface
        const int* globalcellorg = this->grid().globalCellPtr();
        int num_cells = this->gridView().size(0);
        globalcell_.resize(num_cells);
        for(int i=0; i < num_cells; ++i){
            // For grids without global cell numbering, globalcellorg is nullptr
            // and we just use the local cell index as global cell index.
            if (globalcellorg) {
                globalcell_[i] = globalcellorg[i];
            } else {
                globalcell_[i] = i;
            }
        }
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

    void addLgrs()
    { /* do nothing: PolyhedralGrid with LGRs not supported yet! */
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the level logically Cartesian index.
     *        No refinement is supported for AluGrid so it coincides with CartesianIndexMapper.
     */
    const LevelCartesianIndexMapper levelCartesianIndexMapper() const
    { return LevelCartesianIndexMapper(*cartesianIndexMapper_); }

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

    bool gridFromFile() const
    {
        return cellCentroidsFromGrid_;
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
    std::function<std::array<double,FlowBaseVanguard<TypeTag>::dimensionworld>(int)>
    cellCentroids() const
    {
        if (!cellCentroidsFromGrid_) {
            return this->cellCentroids_(this->cartesianIndexMapper(), false);
        } else {
            using UnstructuredGridType = Grid::UnstructuredGridType;
            UnstructuredGridType ugPtr( *this->grid_ );
            double* centroids = ugPtr.cell_centroids;
            int num_cells = this->grid_->size(0);
            return [centroids, num_cells](int i) -> std::array<double, FlowBaseVanguard<TypeTag>::dimensionworld> {
                if (i < 0 || i >= num_cells) return std::array<double, FlowBaseVanguard<TypeTag>::dimensionworld>{};
                return {centroids[i * 3], centroids[i * 3 + 1], centroids[i * 3 + 2]};
            };
        }
    }

    std::vector<int> cellPartition() const
    {
        // not required for this type of grid yet (only from bdaBridge??)
        return {};
    }

protected:
    void createGrids_()
    {
        // Read the grid from file if specified, otherwise use the grid from the ECL file.
        std::string gridFileName = Parameters::Get<Parameters::UnstructuredGridFileName>();
        if (gridFileName.empty()) {
            this->grid_ = std::make_unique<Grid>
                (this->eclState().getInputGrid(),
                this->eclState().fieldProps().porv(true),
                this->edgeConformal());
        } else {
            const auto& input_grid = this->eclState().getInputGrid();
            if (!input_grid.allActive()) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Cannot use unstructured grid file '{}': the Eclipse input grid "
                                      "contains inactive cells ({} active out of {} total). "
                                      "Unstructured grid files require all cells to be active.",
                                      gridFileName,
                                      input_grid.getNumActive(),
                                      input_grid.getCartesianSize()));
            }
            std::cout << "Using unstructured grid from file: " << gridFileName << std::endl;
            this->grid_ = std::make_unique<Grid>(gridFileName.empty() ? "test.txt" : gridFileName);
            cellCentroidsFromGrid_ = true;
            const int numCellsFile = static_cast<int>(this->grid_->size(0));
            const int numCellsEcl  = static_cast<int>(input_grid.getNumActive());
            if (numCellsFile != numCellsEcl) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Cell count mismatch: unstructured grid file '{}' has {} cells, "
                                      "but the Eclipse input grid has {} active cells.",
                                      gridFileName, numCellsFile, numCellsEcl));
            }
        }
        this->cartesianIndexMapper_ =
            std::make_unique<CartesianIndexMapper>(*this->grid_);

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
    bool cellCentroidsFromGrid_ = false;
};

} // namespace Opm

#endif // OPM_POLYHEDRAL_GRID_VANGUARD_HPP
