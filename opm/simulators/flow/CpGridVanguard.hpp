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
 * \copydoc Opm::CpGridVanguard
 */
#ifndef OPM_CPGRID_VANGUARD_HPP
#define OPM_CPGRID_VANGUARD_HPP

#include <opm/common/TimingMacros.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/simulators/flow/FemCpGridCompat.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/GenericCpGridVanguard.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>

#include <array>
#include <functional>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Opm {
template <class TypeTag>
class CpGridVanguard;
}

namespace Opm::Properties {

namespace TTag {
struct CpGridVanguard {
    using InheritsFrom = std::tuple<FlowBaseVanguard>;
};
}

// declare the properties
template<class TypeTag>
struct Vanguard<TypeTag, TTag::CpGridVanguard> {
    using type = CpGridVanguard<TypeTag>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::CpGridVanguard> {
    using type = Dune::CpGrid;
};
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::CpGridVanguard> {
    using type = GetPropType<TypeTag, Properties::Grid>;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::CpGrid as the simulation grid.
 */
template <class TypeTag>
class CpGridVanguard : public FlowBaseVanguard<TypeTag>
                     , public GenericCpGridVanguard<GetPropType<TypeTag, Properties::ElementMapper>,
                                                    GetPropType<TypeTag, Properties::GridView>,
                                                    GetPropType<TypeTag, Properties::Scalar>>
{
    friend class FlowBaseVanguard<TypeTag>;
    using ParentType = FlowBaseVanguard<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using TransmissibilityType = Transmissibility<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>;
    static constexpr int dimensionworld = Grid::dimensionworld;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static constexpr bool waterEnabled = Indices::waterEnabled;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;
private:
    using Element = typename GridView::template Codim<0>::Entity;

public:
    explicit CpGridVanguard(Simulator& simulator)
        : FlowBaseVanguard<TypeTag>(simulator)
    {
        this->checkConsistency();
        this->callImplementationInit();
    }

    int compressedIndexForInteriorLGR(const std::string& lgr_tag, const Connection& conn) const override
    {
        const std::array<int,3> lgr_ijk = {conn.getI(), conn.getJ(), conn.getK()};
        const auto& lgr_level = this->grid().getLgrNameToLevel().at(lgr_tag);
        if (ParentType::lgrMappers_.has_value() == false) {
            ParentType::lgrMappers_.emplace(this->grid().mapLocalCartesianIndexSetsToLeafIndexSet());
        }
        const auto& lgr_dim = this->grid().currentData()[lgr_level]->logicalCartesianSize();
        const auto lgr_cartesian_index = (lgr_ijk[2]*lgr_dim[0]*lgr_dim[1]) + (lgr_ijk[1]*lgr_dim[0]) + (lgr_ijk[0]);
        return ParentType::lgrMappers_.value()[lgr_level].at(lgr_cartesian_index);
    }   
    /*!
     * Checking consistency of simulator
     */
    void checkConsistency()
    {
        const auto& runspec = this->eclState().runspec();
        const auto& config = this->eclState().getSimulationConfig();
        const auto& phases = runspec.phases();

        // check for correct module setup
        if (config.isThermal()) {
            if (getPropValue<TypeTag, Properties::EnableEnergy>() == false) {
                throw std::runtime_error("Input specifies energy while simulator has disabled it, try xxx_energy");
            }
        } else {
            if (getPropValue<TypeTag, Properties::EnableEnergy>() == true) {
                throw std::runtime_error("Input specifies no energy while simulator has energy, try run without _energy");
            }
        }

        if (config.isDiffusive()) {
            if (getPropValue<TypeTag, Properties::EnableDiffusion>() == false) {
                throw std::runtime_error("Input specifies diffusion while simulator has disabled it, try xxx_diffusion");
            }
        }

        if (runspec.micp()) {
            if (getPropValue<TypeTag, Properties::EnableBioeffects>() == false) {
                throw std::runtime_error("Input specifies MICP while simulator has it disabled");
            }
        }

        if (runspec.biof()) {
            if (getPropValue<TypeTag, Properties::EnableBioeffects>() == false) {
                throw std::runtime_error("Input specifies Biofilm while simulator has it disabled");
            }
        }

        if (phases.active(Phase::BRINE)) {
            if (getPropValue<TypeTag, Properties::EnableBrine>() == false) {
                throw std::runtime_error("Input specifies Brine while simulator has it disabled");
            }
        }

        if (phases.active(Phase::POLYMER)) {
            if (getPropValue<TypeTag, Properties::EnablePolymer>() == false) {
                throw std::runtime_error("Input specifies Polymer while simulator has it disabled");
            }
        }

        // checking for correct phases is more difficult TODO!
        if (phases.active(Phase::ZFRACTION)) {
            if (getPropValue<TypeTag, Properties::EnableExtbo>() == false) {
                throw std::runtime_error("Input specifies ExBo while simulator has it disabled");
            }
        }
        if (phases.active(Phase::FOAM)) {
            if (getPropValue<TypeTag, Properties::EnableFoam>() == false) {
                throw std::runtime_error("Input specifies Foam while simulator has it disabled");
            }
        }

        if (phases.active(Phase::SOLVENT)) {
            if (getPropValue<TypeTag, Properties::EnableSolvent>() == false) {
                throw std::runtime_error("Input specifies Solvent while simulator has it disabled");
            }
        }
        if(phases.active(Phase::WATER)){
            if(waterEnabled == false){
                throw std::runtime_error("Input specifies water while simulator has it disabled");
            }
        }
        if(phases.active(Phase::GAS)){
            if(gasEnabled == false){
                throw std::runtime_error("Input specifies gas while simulator has it disabled");
            }
        }
        if(phases.active(Phase::OIL)){
            if(oilEnabled == false){
                throw std::runtime_error("Input specifies oil while simulator has it disabled");
            }
        }

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

    const TransmissibilityType& globalTransmissibility() const
    {
        assert( globalTrans_ != nullptr );
        return *globalTrans_;
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
#if HAVE_MPI
        if (const auto& extPFile = this->externalPartitionFile();
            !extPFile.empty() && (extPFile != "none"))
        {
            this->setExternalLoadBalancer(details::MPIPartitionFromFile { extPFile });
        }

        this->doLoadBalance_(this->edgeWeightsMethod(), this->ownersFirst(),
                             this->addCorners(), this->numOverlap(),
                             this->partitionMethod(), this->serialPartitioning(),
                             this->enableDistributedWells(),
                             this->allow_splitting_inactive_wells_,
                             this->imbalanceTol(),
                             this->gridView(), this->schedule(),
                             this->eclState(), this->parallelWells_,
                             this->numJacobiBlocks(), this->enableEclOutput());
#endif

        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
        this->updateCellThickness_();

#if HAVE_MPI
        this->distributeFieldProps_(this->eclState());
#endif
    }

    /*!
     * \brief Add LGRs and update Leaf Grid View in the simulation grid.
     */
    void addLgrs()
    {
        // Check if input file contains Lgrs. Add them, if any.
        // In a parallel run, this adds the LGRs on the distributed simulation grid.
        if (const auto& lgrs = this->eclState().getLgrs(); lgrs.size() > 0) {
            OpmLog::info("\nAdding LGRs to the grid and updating its leaf grid view");
            this->addLgrsUpdateLeafView(lgrs, lgrs.size(), *this->grid_);

            this->updateGridView_();
            this->updateCellDepths_();
            this->updateCellThickness_();

            if (this->grid_->comm().size()>1) {
                // Add LGRs and update the leaf grid view in the global (undistributed) simulation grid.
                // Purpose: To enable synchronization of cell ids in 'serial mode',
                //          we rely on the "parent-to-children" cell id mapping.
                OpmLog::info("\nAdding LGRs to the global view and updating its leaf grid view");
                this->grid_->switchToGlobalView();
                this->addLgrsUpdateLeafView(lgrs, lgrs.size(), *this->grid_);
                this->grid_->switchToDistributedView();
                this->grid_->syncDistributedGlobalCellIds();
            }
        }
    }

    unsigned int gridEquilIdxToGridIdx(unsigned int elemIndex) const {
        return elemIndex;
    }

    unsigned int gridIdxToEquilGridIdx(unsigned int elemIndex) const {
        return elemIndex;
    }
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
        return this->cellCentroids_(this->cartesianIndexMapper(), true);
    }

    const std::vector<int>& globalCell()
    {
        return this->grid().globalCell();
    }

protected:
    void createGrids_()
    {
        this->doCreateGrids_(this->edgeConformal(), this->eclState());
    }

    void allocTrans() override
    {
        OPM_TIMEBLOCK(allocateTrans);
        globalTrans_.reset(new TransmissibilityType(this->eclState(),
                                                    this->gridView(),
                                                    this->cartesianIndexMapper(),
                                                    this->grid(),
                                                    this->cellCentroids(),
                                                    getPropValue<TypeTag, Properties::EnableEnergy>(),
                                                    getPropValue<TypeTag, Properties::EnableDiffusion>(),
                                                    getPropValue<TypeTag, Properties::EnableDispersion>()));
        globalTrans_->update(false, TransmissibilityType::TransUpdateQuantities::Trans);
    }

    double getTransmissibility(unsigned I, unsigned J) const override
    {
       return globalTrans_->transmissibility(I,J);
    }

#if HAVE_MPI
    const std::string& zoltanParams() const override
    {
        return this->zoltanParams_;
    }

    double zoltanPhgEdgeSizeThreshold() const override
    {
        return this->zoltanPhgEdgeSizeThreshold_;
    }

    const std::string& metisParams() const override
    {
        return this->metisParams_;
    }
#endif

    // removing some connection located in inactive grid cells
    void filterConnections_()
    {
        this->doFilterConnections_(this->schedule());
    }

    // \Note: this globalTrans_ is used for domain decomposition and INIT file output.
    // It only contains trans_ due to permeability and does not contain thermalHalfTrans_,
    // diffusivity_ abd dispersivity_. The main reason is to reduce the memory usage for rank 0
    // during parallel running.
    std::unique_ptr<TransmissibilityType> globalTrans_;
};

} // namespace Opm

#endif // OPM_CPGRID_VANGUARD_HPP
