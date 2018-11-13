/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/grid/utility/extractPvtTableIndex.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/common/ErrorMacros.hpp>

namespace Opm
{
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const Opm::Deck& deck,
                                                     const Opm::EclipseState& eclState,
                                                     std::shared_ptr<MaterialLawManager> materialLawManager,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        init(deck, eclState, materialLawManager, grid.number_of_cells, grid.global_cell, grid.cartdims,
             init_rock);
    }

#ifdef HAVE_OPM_GRID
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const Opm::Deck& deck,
                                                     const Opm::EclipseState& eclState,
                                                     const Dune::CpGrid& grid,
                                                     const bool init_rock )
    {
        auto materialLawManager = std::make_shared<MaterialLawManager>();
        unsigned number_of_cells = grid.size(0);
        std::vector<int> compressedToCartesianIdx(number_of_cells);
        for (unsigned cellIdx = 0; cellIdx < number_of_cells; ++cellIdx) {
            compressedToCartesianIdx[cellIdx] = grid.globalCell()[cellIdx];
        }
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);
        init(deck, eclState, materialLawManager, grid.numCells(), static_cast<const int*>(&grid.globalCell()[0]),
             static_cast<const int*>(&grid.logicalCartesianSize()[0]),
             init_rock);
    }
#endif

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const Opm::Deck& deck,
                                                     const Opm::EclipseState& eclState,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        auto materialLawManager = std::make_shared<MaterialLawManager>();
        std::vector<int> compressedToCartesianIdx(grid.number_of_cells);
        for (int cellIdx = 0; cellIdx < grid.number_of_cells; ++cellIdx) {
            if (grid.global_cell) {
                compressedToCartesianIdx[cellIdx] = grid.global_cell[cellIdx];
            }
            else {
                compressedToCartesianIdx[cellIdx] = cellIdx;
            }
        }
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);
        init(deck, eclState, materialLawManager, grid.number_of_cells, grid.global_cell, grid.cartdims,
             init_rock);
    }

#ifdef HAVE_OPM_GRID
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const Opm::Deck& deck,
                                                     const Opm::EclipseState& eclState,
                                                     std::shared_ptr<MaterialLawManager> materialLawManager,
                                                     const Dune::CpGrid& grid,
                                                     const bool init_rock )
    {
        init(deck, eclState, materialLawManager, grid.numCells(), static_cast<const int*>(&grid.globalCell()[0]),
             static_cast<const int*>(&grid.logicalCartesianSize()[0]),
             init_rock);
    }
#endif

/// Constructor for properties on a subgrid
BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const BlackoilPropsAdFromDeck& props,
                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                 const int number_of_cells)
    : rock_(number_of_cells), satprops_(new SaturationPropsFromDeck())
{
    const int original_size = props.cellPvtRegionIdx_.size();
    if (number_of_cells > original_size) {
        OPM_THROW(std::runtime_error, "The number of cells is larger than the one of the original grid!");
    }
    if (number_of_cells < 0) {
        OPM_THROW(std::runtime_error, "The number of cells is has to be larger than 0.");
    }

    materialLawManager_ = materialLawManager;

    // Copy properties that do not depend on the postion within the grid.
    phase_usage_      = props.phase_usage_;
    vap1_             = props.vap1_;
    vap2_             = props.vap2_;
    vap_satmax_guard_ = props.vap_satmax_guard_;
    // For data that is dependant on the subgrid we simply allocate space
    // and initialize with obviously bogus numbers.
    cellPvtRegionIdx_.resize(number_of_cells, std::numeric_limits<int>::min());
    satprops_->init(phase_usage_, materialLawManager_);
}

    /// Initializes the properties.
    void BlackoilPropsAdFromDeck::init(const Opm::Deck& deck,
                                       const Opm::EclipseState& eclState,
                                       std::shared_ptr<MaterialLawManager> materialLawManager,
                                       int number_of_cells,
                                       const int* global_cell,
                                       const int* cart_dims,
                                       const bool init_rock)
    {
        materialLawManager_ = materialLawManager;

        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, eclState, number_of_cells, global_cell);

        if (init_rock){
            rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }

        phase_usage_ = phaseUsageFromDeck(deck);

        if (!FluidSystem::isInitialized()) {
            // make sure that we don't initialize the fluid system twice
            FluidSystem::initFromDeck(deck, eclState);
        }

        // Oil vaporization controls (kw VAPPARS)
        vap1_ = vap2_ = 0.0;
        if (deck.hasKeyword("VAPPARS") && deck.hasKeyword("VAPOIL") && deck.hasKeyword("DISGAS")) {
            vap1_ = deck.getKeyword("VAPPARS").getRecord(0).getItem(0).get< double >(0);
            vap2_ = deck.getKeyword("VAPPARS").getRecord(0).getItem(1).get< double >(0);
        } else if (deck.hasKeyword("VAPPARS")) {
            OPM_THROW(std::runtime_error, "Input has VAPPARS, but missing VAPOIL and/or DISGAS\n");
        }
        satOilMax_.resize(number_of_cells, 0.0);

        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        satprops_.reset(ptr);
        ptr->init(deck, materialLawManager_);

        if (phase_usage_.num_phases != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck() - "
                      "Inconsistent number of phases in pvt data (" << phase_usage_.num_phases
                      << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
        vap_satmax_guard_ = 0.01;
    }

    ////////////////////////////
    //      Rock interface    //
    ////////////////////////////

    /// \return   D, the number of spatial dimensions.
    int BlackoilPropsAdFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int BlackoilPropsAdFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropsAdFromDeck::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropsAdFromDeck::permeability() const
    {
        return rock_.permeability();
    }


    ////////////////////////////
    //      Fluid interface   //
    ////////////////////////////

    /// \return   Number of active phases (also the number of components).
    int BlackoilPropsAdFromDeck::numPhases() const
    {
        return phase_usage_.num_phases;
    }

    /// \return   Object describing the active phases.
    PhaseUsage BlackoilPropsAdFromDeck::phaseUsage() const
    {
        return phase_usage_;
    }

    /// Saturation update for hysteresis behavior.
    /// \param[in]  cells       Array of n cell indices to be associated with the saturation values.
    void BlackoilPropsAdFromDeck::updateSatHyst(const std::vector<double>& saturation,
                                                const std::vector<int>& cells)
    {
        const int n = cells.size();
        satprops_->updateSatHyst(n, cells.data(), saturation.data());
    }

    /// Set gas-oil hysteresis parameters
    /// \param[in]  pcswmdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::pcSwMdc(...))
    /// \param[in]  krnswdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::krnSwMdc(...))
    void BlackoilPropsAdFromDeck::setGasOilHystParams(const std::vector<double>& pcswmdc,
                             const std::vector<double>& krnswdc,
                             const std::vector<int>& cells)
    {
        const size_t n = cells.size();
        assert(pcswmdc.size() == n);
        assert(krnswdc.size() == n);
        satprops_->setGasOilHystParams(n, cells.data(), pcswmdc.data(), krnswdc.data());
    }

    /// Get gas-oil hysteresis parameters
    /// \param[in]  pcswmdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::pcSwMdc(...))
    /// \param[in]  krnswdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::krnSwMdc(...))
    void BlackoilPropsAdFromDeck::getGasOilHystParams(std::vector<double>& pcswmdc,
                             std::vector<double>& krnswdc,
                             const std::vector<int>& cells) const
    {
        const size_t n = cells.size();
        pcswmdc.resize(n);
        krnswdc.resize(n);
        satprops_->getGasOilHystParams(n, cells.data(), pcswmdc.data(), krnswdc.data());
    }

    /// Set oil-water hysteresis parameters
    /// \param[in]  pcswmdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::pcSwMdc(...))
    /// \param[in]  krnswdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::krnSwMdc(...))
    void BlackoilPropsAdFromDeck::setOilWaterHystParams(const std::vector<double>& pcswmdc,
                               const std::vector<double>& krnswdc,
                               const std::vector<int>& cells)
    {
        const size_t n = cells.size();
        assert(pcswmdc.size() == n);
        assert(krnswdc.size() == n);
        satprops_->setOilWaterHystParams(n, cells.data(), pcswmdc.data(), krnswdc.data());
    }

    /// Get oil-water hysteresis parameters
    /// \param[in]  pcswmdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::pcSwMdc(...))
    /// \param[in]  krnswdc  Vector of hysteresis parameters (@see EclHysteresisTwoPhaseLawParams::krnSwMdc(...))
    void BlackoilPropsAdFromDeck::getOilWaterHystParams(std::vector<double>& pcswmdc,
                               std::vector<double>& krnswdc,
                               const std::vector<int>& cells) const
    {
        const size_t n = cells.size();
        pcswmdc.resize(n);
        krnswdc.resize(n);
        satprops_->getOilWaterHystParams(n, cells.data(), pcswmdc.data(), krnswdc.data());
    }

    /// Update for max oil saturation.
    void BlackoilPropsAdFromDeck::updateSatOilMax(const std::vector<double>& saturation)
    {
        const int n = satOilMax_.size();
        const int np = phase_usage_.num_phases;
        const int posOil = phase_usage_.phase_pos[Oil];
        const double* s = saturation.data();
        for (int i=0; i<n; ++i) {
            if (satOilMax_[i] < s[np*i+posOil]) {
                satOilMax_[i] = s[np*i+posOil];
            }
        }
    }

    // Get max oil saturation vector
    const std::vector<double>& BlackoilPropsAdFromDeck::satOilMax() const
    {
        return satOilMax_;
    }

    // Set max oil saturation vector
    void BlackoilPropsAdFromDeck::setSatOilMax(const std::vector<double>& max_sat) {
        assert(satOilMax_.size() == max_sat.size());
        satOilMax_ = max_sat;
    }

    /// Set capillary pressure scaling according to pressure diff. and initial water saturation.
    /// \param[in]  saturation Array of n*numPhases saturation values.
    /// \param[in]  pc         Array of n*numPhases capillary pressure values.
    void BlackoilPropsAdFromDeck::setSwatInitScaling(const std::vector<double>& saturation,
                                           const std::vector<double>& pc)
    {
        const int nc = rock_.numCells();
        const int numActivePhases = numPhases();
        for (int i = 0; i < nc; ++i) {
            double pcow = pc[numActivePhases*i + phase_usage_.phase_pos[Water]];
            double swat = saturation[numActivePhases*i + phase_usage_.phase_pos[Water]];
            satprops_->swatInitScaling(i, pcow, swat);
        }
    }

} // namespace Opm
