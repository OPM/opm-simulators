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
 * \copydoc Opm::EclAluGridVanguard
 */
#ifndef EWOMS_ECL_ALU_GRID_VANGUARD_HH
#define EWOMS_ECL_ALU_GRID_VANGUARD_HH

#include <dune/alugrid/common/fromtogridfactory.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>

#include <ebos/alucartesianindexmapper.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/ecltransmissibility.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>

#include <opm/simulators/utils/ParallelEclipseState.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <tuple>
#include <vector>

namespace Opm {
template <class TypeTag>
class EclAluGridVanguard;

} // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct EclAluGridVanguard {
    using InheritsFrom = std::tuple<EclBaseVanguard>;
};
}

// declare the properties
template<class TypeTag>
struct Vanguard<TypeTag, TTag::EclAluGridVanguard> {
    using type = Opm::EclAluGridVanguard<TypeTag>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::EclAluGridVanguard> {
#if HAVE_MPI
    using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>; 
#else    
    using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI     
};
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::EclAluGridVanguard> {
    using type = Dune::CpGrid;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::ALUGrid as the simulation grid.
 */
template <class TypeTag>
class EclAluGridVanguard : public EclBaseVanguard<TypeTag>
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
    using TransmissibilityType = EclTransmissibility<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>;
    using Factory = Dune::FromToGridFactory<Grid>;

    static constexpr int dimension = Grid::dimension;
    static constexpr int dimensionworld = Grid::dimensionworld;

    EclAluGridVanguard(Simulator& simulator)
        : EclBaseVanguard<TypeTag>(simulator)
    { 
      this->mpiRank = EclGenericVanguard::comm().rank();
      this->callImplementationInit();
    }

    ~EclAluGridVanguard() = default;

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
    { return *equilGrid_; }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    {
        delete equilCartesianIndexMapper_;
        equilCartesianIndexMapper_ = nullptr;

        delete equilGrid_;
        equilGrid_ = nullptr;
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
        auto gridView = grid().leafGridView();
        auto dataHandle = cartesianIndexMapper_->dataHandle(gridView);
        grid().loadBalance(*dataHandle);

        // communicate non-interior cells values
        grid().communicate(*dataHandle,
                           Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication );

       if (grid().size(0))
       {
           globalTrans_ = std::make_unique<TransmissibilityType>(this->eclState(),
                                                                 this->gridView(),
                                                                 this->cartesianIndexMapper(),
                                                                 this->grid(),
                                                                 this->cellCentroids(),
                                                                 getPropValue<TypeTag,
                                                                 Properties::EnableEnergy>(),
                                                                 getPropValue<TypeTag,
                                                                 Properties::EnableDiffusion>(),
                                                                 getPropValue<TypeTag,
                                                                 Properties::EnableDispersion>());
            // Re-ordering  for ALUGrid
            globalTrans_->update(false, [&](unsigned int i) { return gridEquilIdxToGridIdx(i);});
        }
        
    }

    template<class DataHandle>
    void scatterData(DataHandle& /*handle*/) const
    {
    // not existing for this type of grid yet
    }

    template<class DataHandle>
    void gatherData(DataHandle& /*handle*/) const
    {
    // not existing for this type of grid yet
    }

    template<class DataHandle, class InterfaceType, class CommunicationDirection>
    void communicate (DataHandle& /*data*/, InterfaceType /*iftype*/,
                      CommunicationDirection /*dir*/) const
    {
    // not existing for this type of grid yet
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

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     */
    const EquilCartesianIndexMapper& equilCartesianIndexMapper() const
    { return *equilCartesianIndexMapper_; }

    /*!
     * \brief Get function to query cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     * It is a function return the centroid for the given element
     * index.
     */
    std::function<std::array<double,dimensionworld>(int)>
    cellCentroids() const
    {
        return this->cellCentroids_(this->cartesianIndexMapper(), false);
    }

    const TransmissibilityType& globalTransmissibility() const
    {
        assert( globalTrans_ != nullptr );
        return *globalTrans_;
    }

    void releaseGlobalTransmissibility()
    {
        globalTrans_.reset();
    }

    const std::vector<int>& globalCell()
    {
        return cartesianCellId_;
    }

    std::vector<int> cellPartition() const
    {      
        // not required for this type of grid yet
        return {};
    }
    
    unsigned int gridEquilIdxToGridIdx(unsigned int elemIndex) const {
        return equilGridToGrid_[elemIndex];
    }

    unsigned int gridIdxToEquilGridIdx(unsigned int elemIndex) const {
        return ordering_[elemIndex];
    }

protected:
    void createGrids_()
    {
        // we use separate grid objects: one for the calculation of the initial condition
        // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
        // cannot cope with arbitrary Dune grids and is also allergic to distributed
        // grids.

        /////
        // create the EQUIL grid
        /////
        const EclipseGrid* input_grid = nullptr;
        std::vector<double> global_porv;
        // At this stage the ParallelEclipseState instance is still in global
        // view; on rank 0 we have undistributed data for the entire grid, on
        // the other ranks the EclipseState is empty.
        if (mpiRank == 0) {
            // Processing grid
            input_grid = &this->eclState().getInputGrid();
            global_porv = this->eclState().fieldProps().porv(true);
            OpmLog::info("\nProcessing grid");
        }

#if HAVE_MPI
        this->equilGrid_ = std::make_unique<Dune::CpGrid>(EclGenericVanguard::comm());
#else
        this->equilGrid_ = std::make_unique<Dune::CpGrid>();
#endif
        // Note: removed_cells is guaranteed to be empty on ranks other than 0.
        auto removed_cells =
            this->equilGrid_->processEclipseFormat(input_grid,
                                                   &this->eclState(),
                                                   /*isPeriodic=*/false,
                                                   /*flipNormals=*/false,
                                                   /*clipZ=*/false);

        cartesianCellId_ = this->equilGrid_->globalCell();

        for (unsigned i = 0; i < dimension; ++i)
            cartesianDimension_[i] = this->equilGrid_->logicalCartesianSize()[i];

        equilCartesianIndexMapper_ = std::make_unique<EquilCartesianIndexMapper>(*equilGrid_);

        /////
        // create the simulation grid
        /////

        factory_ = std::make_unique<Factory>();
        grid_ = factory_->convert(*equilGrid_, cartesianCellId_, ordering_);
        OpmLog::warning("Space Filling Curve Ordering is not yet supported: DISABLE_ALUGRID_SFC_ORDERING is enabled");
        equilGridToGrid_.resize(ordering_.size());
        for (std::size_t index = 0; index < ordering_.size(); ++index) {
            equilGridToGrid_[ordering_[index]] = index;
        }

        cartesianIndexMapper_ = std::make_unique<CartesianIndexMapper>(*grid_, cartesianDimension_, cartesianCellId_);
        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
        this->updateCellThickness_();
    }

    void filterConnections_()
    {
        // not handling the removal of completions for this type of grid yet.
    }

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<EquilGrid> equilGrid_;
    std::vector<int> cartesianCellId_;
    std::vector<unsigned int> ordering_;
    std::vector<unsigned int> equilGridToGrid_;
    std::array<int,dimension> cartesianDimension_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<EquilCartesianIndexMapper> equilCartesianIndexMapper_;
    std::unique_ptr<Factory> factory_;
    std::unique_ptr<TransmissibilityType> globalTrans_;
    int mpiRank;
};

} // namespace Opm

#endif
