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

#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/common/UnstructuredGridPartitioning.hpp>
#include <opm/grid/polyhedralgrid/levelcartesianindexmapper.hh>

#include <opm/models/common/multiphasebaseproperties.hh>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>

#if HAVE_MPI
#include <opm/simulators/flow/partitionCells.hpp>
#include <opm/simulators/utils/MPISerializer.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/utils/PropsDataHandle.hpp>
#endif

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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
        this->updateGlobalCell_();
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
    { return equilGrid_ ? *equilGrid_ : *grid_; }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    {
        // In a serial run the EQUIL grid is the simulation grid, which must stay.
        equilGrid_.reset();
        equilCartesianIndexMapper_.reset();
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
#if HAVE_MPI
        if (this->grid_->comm().size() > 1) {
            this->doLoadBalance_();

            this->updateGridView_();
            // The mapping of the global grid, made when it was created, is not valid anymore.
            this->cartesianToCompressed_.clear();
            this->updateCartesianToCompressedMapping_();
            this->updateCellThickness_();
            this->distributeFieldProps_();
            this->updateCellDepths_();
            this->updateGlobalCell_();
        }
#endif
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
     * In a serial run that's the same as cartesianIndexMapper(). In a
     * parallel run it is the mapper of the global grid, which only exists on
     * rank 0.
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const
    { return equilCartesianIndexMapper_ ? *equilCartesianIndexMapper_ : *cartesianIndexMapper_; }

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
        globalTrans_.reset();
    }

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const TransmissibilityType& globalTransmissibility() const
    {
        // In a parallel run, the global transmissibilities of rank 0 are
        // computed before the grid is distributed.
        return globalTrans_ ? *globalTrans_ : simulator_.problem().eclTransmissibilities();
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
        if (eclCentroids_.empty()) {
            // A serial run, or rank 0 before the grid has been distributed: use the input grid.
            return this->cellCentroids_(this->cartesianIndexMapper(), false);
        }
        // The processes other than rank 0 do not have the input grid, so the
        // centers have been distributed with the grid.
        return [this](int elemIdx) { return eclCentroids_[elemIdx]; };
    }

    std::vector<int> cellPartition() const
    {
        // not required for this type of grid yet (only from bdaBridge??)
        return {};
    }

protected:
    void createGrids_()
    {
        const auto& comm = FlowGenericVanguard::comm();
        const bool parallel = comm.size() > 1;
        const bool isRoot = comm.rank() == 0;

        // In a parallel run only the root process has the input grid.
        // The other processes get an empty grid, which receives its cells in loadBalance().
        std::array<int, 3> cartDims{};
        if (parallel) {
            if (isRoot) {
                const auto& inputGrid = this->eclState().getInputGrid();
                cartDims = {static_cast<int>(inputGrid.getNX()),
                            static_cast<int>(inputGrid.getNY()),
                            static_cast<int>(inputGrid.getNZ())};
            }
            comm.broadcast(cartDims.data(), 3, 0);
        }

        if (!parallel || isRoot) {
            this->grid_ = createGridFromInput_();
            if (parallel) {
                // The simulation grid is distributed, so the output needs its own copy of the global grid.
                this->equilGrid_ = createGridFromInput_();
                this->equilCartesianIndexMapper_ =
                    std::make_unique<CartesianIndexMapper>(*this->equilGrid_);
            }

            // The grid can have fewer cells than the deck has active cells, for instance if cells are
            // removed because of MINPV. The properties of the cells are given for the active cells,
            // so they must follow the cells of the grid, as they do for the CpGrid. Only the process
            // that has the input grid, which is the one that has the properties as well, can do this.
            auto& eclState = this->eclState();
            eclState.reset_actnum(UgGridHelpers::createACTNUM(*this->grid_));
            const int* globalCell = this->grid_->globalCellPtr();
            eclState.set_active_indices(std::vector<int>(globalCell, globalCell + this->grid_->size(0)));
        }
        else {
            typename Grid::UnstructuredGridPtr empty(allocate_grid(3, 0, 0, 0, 0, 0));
            std::copy(cartDims.begin(), cartDims.end(), empty->cartdims);
            this->grid_ = std::make_unique<Grid>(std::move(empty));
        }

        this->cartesianIndexMapper_ =
            std::make_unique<CartesianIndexMapper>(*this->grid_);

        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
    }

    std::unique_ptr<Grid> createGridFromInput_() const
    {
        return std::make_unique<Grid>(this->eclState().getInputGrid(),
                                      this->eclState().fieldProps().porv(true),
                                      this->edgeConformal());
    }


    void updateGlobalCell_()
    {
        // A grid without global_cell, such as the empty grid of a process that
        // has no cells yet, uses the compressed index.
        const int* globalCell = this->grid().globalCellPtr();
        globalcell_.resize(this->gridView().size(0));
        for (int i = 0; i < static_cast<int>(globalcell_.size()); ++i) {
            globalcell_[i] = globalCell ? globalCell[i] : i;
        }
    }

#if HAVE_MPI
    // Distribute the simulation grid and the field properties.
    void doLoadBalance_()
    {
        auto* eclState = dynamic_cast<ParallelEclipseState*>(&this->eclState());
        if (eclState == nullptr) {
            const std::string message = "Parallel simulator setup is incorrect as "
                                        "it does not use ParallelEclipseState";
            OpmLog::error(message);
            throw std::invalid_argument { message };
        }

        const auto& comm = this->grid_->comm();
        const bool isRoot = comm.rank() == 0;

        // NNCs that depend on the grid processing exist on rank 0 only.
        if (eclState->aquifer().hasNumericalAquifer()) {
            auto nncInput = eclState->getInputNNC();
            Parallel::MpiSerializer ser(comm);
            ser.broadcast(Parallel::RootRank{0}, nncInput);
            if (!isRoot) {
                eclState->setInputNNC(nncInput);
            }
        }
        bool hasPinchNnc = eclState->hasPinchNNC();
        comm.broadcast(&hasPinchNnc, 1, 0);
        if (hasPinchNnc) {
            auto pinchNnc = eclState->getPinchNNC();
            Parallel::MpiSerializer ser(comm);
            ser.broadcast(Parallel::RootRank{0}, pinchNnc);
            if (!isRoot) {
                eclState->setPinchNNC(std::move(pinchNnc));
            }
        }

        // The transmissibilities of the whole grid, which the output needs, are computed before it is distributed.
        if (isRoot && this->enableEclOutput()) {
            this->allocTrans_();
        }

        std::vector<int> parts;
        if (isRoot) {
            parts = this->computePartition_();
        }

        PropsDataHandle<Grid> handle { *this->grid_, *eclState };
        this->grid_->loadBalance(handle, parts, this->ownersFirst(),
                                 /* addCornerCells = */ false, this->numOverlap());

        // The cell centers of the input grid, which only rank 0 has, are needed on all processes.
        eclCentroids_.assign(this->grid_->size(0), {});
        CentroidHandle centroidHandle { *this->grid_, isRoot ? &eclState->getInputGrid() : nullptr, eclCentroids_ };
        this->grid_->scatterData(centroidHandle);

        this->setParallelWells_();
    }

    // Moves the cell centers of the input grid from rank 0 to the cells of the distributed grid.
    class CentroidHandle
    {
    public:
        using DataType = double;

        CentroidHandle(const Grid& grid,
                       const EclipseGrid* inputGrid,
                       std::vector<std::array<double, 3>>& centroids)
            : grid_(grid), inputGrid_(inputGrid), centroids_(centroids)
        {}

        bool contains(int, int codim) const { return codim == 0; }
        bool fixedSize(int, int) const { return true; }

        template <class Entity>
        std::size_t size(const Entity&) const { return 3; }

        // Called for the cells of the global grid, on rank 0.
        template <class Buffer, class Entity>
        void gather(Buffer& buffer, const Entity& entity) const
        {
            const int cell = grid_.leafIndexSet().index(entity);
            const int* globalCell = grid_.globalCellPtr();
            for (const double value : inputGrid_->getCellCenter(globalCell ? globalCell[cell] : cell)) {
                buffer.write(value);
            }
        }

        // Called for the cells of the distributed grid.
        template <class Buffer, class Entity>
        void scatter(Buffer& buffer, const Entity& entity, std::size_t)
        {
            for (double& value : centroids_[grid_.leafIndexSet().index(entity)]) {
                buffer.read(value);
            }
        }

    private:
        const Grid& grid_;
        const EclipseGrid* inputGrid_;
        std::vector<std::array<double, 3>>& centroids_;
    };

    void allocTrans_()
    {
        globalTrans_ = std::make_unique<TransmissibilityType>
            (this->eclState(),
             this->gridView(),
             this->cartesianIndexMapper(),
             this->grid(),
             this->cellCentroids(),
             getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal ||
             getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::SequentialImplicitThermal,
             getPropValue<TypeTag, Properties::EnableDiffusion>(),
             getPropValue<TypeTag, Properties::EnableDispersion>());
        globalTrans_->update(false, TransmissibilityType::TransUpdateQuantities::Trans);
    }

    // The process that owns each cell of the global grid, only on rank 0.
    std::vector<int> computePartition_() const
    {
        const int numCells = this->grid_->size(0);
        const int numParts = this->grid_->comm().size();

        std::vector<int> parts;
        if (const auto& file = this->externalPartitionFile(); !file.empty() && file != "none") {
            parts = partitionCellsFromFile(file, numCells).first;
            if (*std::max_element(parts.begin(), parts.end()) >= numParts) {
                throw std::invalid_argument { "The partition file has more parts than there are processes." };
            }
        }
        else if (this->partitionMethod() == Dune::PartitionMethod::simple || !Dune::unstructured::metisAvailable()) {
            if (this->partitionMethod() != Dune::PartitionMethod::simple) {
                OpmLog::warning("The grid library was built without METIS, so the grid is partitioned in "
                                "blocks of consecutive cells.");
            }
            parts = partitionCellsSimple(numCells, numParts).first;
        }
        else {
            // Zoltan is only available for the CpGrid, so the other methods use METIS.
            if (this->partitionMethod() != Dune::PartitionMethod::metis) {
                OpmLog::info("The partition method is not available for PolyhedralGrid, using METIS.");
            }
            parts = Dune::unstructured::partitionWithMetis(static_cast<const UnstructuredGrid&>(*this->grid_),
                                                           numParts, this->imbalanceTol());
        }
        this->keepWellsTogether_(parts);

        std::vector<int> size(numParts, 0);
        for (const int part : parts) {
            ++size[part];
        }
        if (std::ranges::any_of(size, [](const int n) { return n == 0; })) {
            throw std::runtime_error { "At least one process has no cells. Use fewer processes, or "
                                       "a different partition method." };
        }
        return parts;
    }

    // The cells that a well is connected to must be on one process, unless distributed wells are allowed.
    void keepWellsTogether_(std::vector<int>& parts) const
    {
        if (this->enableDistributedWells()) {
            return;
        }
        const auto* globalCell = this->grid_->globalCellPtr();
        std::unordered_map<int, int> cartesianToCompressed;
        for (int cell = 0; cell < static_cast<int>(parts.size()); ++cell) {
            cartesianToCompressed.emplace(globalCell ? globalCell[cell] : cell, cell);
        }

        // The compressed indices of the cells of each well.
        std::vector<std::vector<int>> wellCells;
        const auto& possibleFuture = this->schedule().getPossibleFutureConnections();
        for (const auto& well : this->schedule().getWellsatEnd()) {
            std::set<int> cartesian;
            for (const auto& connection : well.getConnections()) {
                cartesian.insert(connection.global_index());
            }
            if (const auto it = possibleFuture.find(well.name()); it != possibleFuture.end()) {
                cartesian.insert(it->second.begin(), it->second.end());
            }
            std::vector<int> cells;
            for (const int index : cartesian) {
                if (const auto it = cartesianToCompressed.find(index); it != cartesianToCompressed.end()) {
                    cells.push_back(it->second);
                }
            }
            if (cells.size() > 1) {
                wellCells.push_back(std::move(cells));
            }
        }

        // Wells may share cells, so repeat until nothing moves.
        for (bool changed = true; changed;) {
            changed = false;
            for (const auto& cells : wellCells) {
                const int target = parts[cells.front()];
                for (const int cell : cells) {
                    if (parts[cell] != target) {
                        parts[cell] = target;
                        changed = true;
                    }
                }
            }
        }
    }

    // Which wells have connections on this process.
    void setParallelWells_()
    {
        std::unordered_set<int> cellsOnProcess;
        const auto* globalCell = this->grid_->globalCellPtr();
        for (int cell = 0; cell < this->grid_->size(0); ++cell) {
            cellsOnProcess.insert(globalCell ? globalCell[cell] : cell);
        }

        this->parallelWells_.clear();
        for (const auto& well : this->schedule().getWellsatEnd()) {
            const auto& connections = well.getConnections();
            const bool local = std::any_of(connections.begin(), connections.end(),
                                           [&](const auto& connection)
                                           { return cellsOnProcess.count(connection.global_index()) > 0; });
            this->parallelWells_.emplace_back(well.name(), local);
        }
        std::ranges::sort(this->parallelWells_);
    }

    void distributeFieldProps_()
    {
        auto* eclState = dynamic_cast<ParallelEclipseState*>(&this->eclState());
        // The field properties keep a pointer to the mapper.
        levelCartesianIndexMapper_ = std::make_unique<LevelCartesianIndexMapper>(*this->cartesianIndexMapper_);
        eclState->resetCartesianMapper(levelCartesianIndexMapper_.get());
        eclState->switchToDistributedProps();
    }
#endif

    void filterConnections_()
    {
        // not handling the removal of completions for this type of grid yet.
    }

    Simulator& simulator_;

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<LevelCartesianIndexMapper> levelCartesianIndexMapper_;
    // The global grid on rank 0 of a parallel run, which the output uses.
    std::unique_ptr<Grid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;
    // The transmissibilities of the global grid on rank 0, computed before the grid is distributed.
    std::unique_ptr<TransmissibilityType> globalTrans_;
    // The cell centers of the input grid for the cells of the distributed grid, in a parallel run.
    std::vector<std::array<double, 3>> eclCentroids_;
    //CartesianIndexMapperPointer cartesianIndexMapper_;

    std::unordered_set<std::string> defunctWellNames_;
    std::vector<int> globalcell_;
};

} // namespace Opm

#endif // OPM_POLYHEDRAL_GRID_VANGUARD_HPP
