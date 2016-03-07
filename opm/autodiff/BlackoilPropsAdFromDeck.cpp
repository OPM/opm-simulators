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
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/extractPvtTableIndex.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/common/ErrorMacros.hpp>

namespace Opm
{
    // Making these typedef to make the code more readable.
    typedef BlackoilPropsAdFromDeck::ADB ADB;
    typedef BlackoilPropsAdFromDeck::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
                                                     std::shared_ptr<MaterialLawManager> materialLawManager,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        init(deck, eclState, materialLawManager, grid.number_of_cells, grid.global_cell, grid.cartdims,
             init_rock);
    }

#ifdef HAVE_DUNE_CORNERPOINT
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
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
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
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

#ifdef HAVE_DUNE_CORNERPOINT
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
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
    oilPvt_           = props.oilPvt_;
    gasPvt_           = props.gasPvt_;
    waterPvt_         = props.waterPvt_;
    phase_usage_      = props.phase_usage_;
    surfaceDensity_   = props.surfaceDensity_;
    vap1_             = props.vap1_;
    vap2_             = props.vap2_;
    vap_satmax_guard_ = props.vap_satmax_guard_;
    // For data that is dependant on the subgrid we simply allocate space
    // and initialize with obviously bogus numbers.
    cellPvtRegionIdx_.resize(number_of_cells, std::numeric_limits<int>::min());
    satprops_->init(phase_usage_, materialLawManager_);
}

    /// Initializes the properties.
    void BlackoilPropsAdFromDeck::init(Opm::DeckConstPtr deck,
                                       Opm::EclipseStateConstPtr eclState,
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

        gasPvt_ = std::make_shared<GasPvt>();
        oilPvt_ = std::make_shared<OilPvt>();
        waterPvt_ = std::make_shared<WaterPvt>();

        gasPvt_->initFromDeck(deck, eclState);
        oilPvt_->initFromDeck(deck, eclState);
        waterPvt_->initFromDeck(deck, eclState);

        // Surface densities. Accounting for different orders in eclipse and our code.
        const auto& densityKeyword = deck->getKeyword("DENSITY");
        int numRegions = densityKeyword.size();
        auto tables = eclState->getTableManager();

        surfaceDensity_.resize(numRegions);
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            if (phase_usage_.phase_used[Liquid]) {
                surfaceDensity_[regionIdx][phase_usage_.phase_pos[Liquid]]
                    = densityKeyword.getRecord(regionIdx).getItem("OIL").getSIDouble(0);
            }
            if (phase_usage_.phase_used[Aqua]) {
                surfaceDensity_[regionIdx][phase_usage_.phase_pos[Aqua]]
                    = densityKeyword.getRecord(regionIdx).getItem("WATER").getSIDouble(0);
            }
            if (phase_usage_.phase_used[Vapour]) {
                surfaceDensity_[regionIdx][phase_usage_.phase_pos[Vapour]]
                    = densityKeyword.getRecord(regionIdx).getItem("GAS").getSIDouble(0);
            }
        }

        // Oil vaporization controls (kw VAPPARS)
        vap1_ = vap2_ = 0.0;
        if (deck->hasKeyword("VAPPARS") && deck->hasKeyword("VAPOIL") && deck->hasKeyword("DISGAS")) {
            vap1_ = deck->getKeyword("VAPPARS").getRecord(0).getItem(0).get< double >(0);
            vap2_ = deck->getKeyword("VAPPARS").getRecord(0).getItem(1).get< double >(0);
            satOilMax_.resize(number_of_cells, 0.0);
        } else if (deck->hasKeyword("VAPPARS")) {
            OPM_THROW(std::runtime_error, "Input has VAPPARS, but missing VAPOIL and/or DISGAS\n");
        }

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

    // ------ Density ------

    /// Densities of stock components at surface conditions.
    /// \param[in] phaseIdx
    /// \param[in] cells  Array of n cell indices to be associated with the pressure values.
    /// \return Array of n density values for phase given by phaseIdx.
    V BlackoilPropsAdFromDeck::surfaceDensity(const int phaseIdx, const Cells& cells) const
    {
        assert( !(phaseIdx > numPhases()));
        const int n = cells.size();
        V rhos = V::Zero(n);
        for (int cellIdx = 0; cellIdx < n; ++cellIdx) {
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            const auto* rho = &surfaceDensity_[pvtRegionIdx][0];
            rhos[cellIdx] = rho[phaseIdx];
        }
        return rhos;
    }


    // ------ Viscosity ------


    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muWat(const ADB& pw,
                                       const ADB& T,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not active.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V mu(n);
        V dmudp(n);

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = pw.value()[i];
            TLad.value = T.value()[i];

            const LadEval& muLad = waterPvt_->viscosity(pvtRegionIdx, TLad, pLad);

            mu[i] = muLad.value;
            dmudp[i] = muLad.derivatives[0];
        }

        if (pw.derivative().empty()) {
            return ADB::constant(std::move(mu));
        } else {
            ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
            const int num_blocks = pw.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                fastSparseProduct(dmudp_diag, pw.derivative()[block], jacs[block]);
            }
            return ADB::function(std::move(mu), std::move(jacs));
        }
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muOil(const ADB& po,
                                       const ADB& T,
                                       const ADB& rs,
                                       const std::vector<PhasePresence>& cond,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        enum PressureRsEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureRsEvalTag, /*size=*/2> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RsLad = 0.0;

        pLad.derivatives[0] = 1.0;
        RsLad.derivatives[1] = 1.0;

        LadEval muLad;
        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = po.value()[i];
            TLad.value = T.value()[i];

            if (cond[i].hasFreeGas()) {
                muLad = oilPvt_->saturatedViscosity(pvtRegionIdx, TLad, pLad);
            }
            else {
                RsLad.value = rs.value()[i];
                muLad = oilPvt_->viscosity(pvtRegionIdx, TLad, pLad, RsLad);
            }

            mu[i] = muLad.value;
            dmudp[i] = muLad.derivatives[0];
            dmudr[i] = muLad.derivatives[1];
        }

        ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
        ADB::M dmudr_diag(dmudr.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dmudp_diag, po.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dmudr_diag, rs.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(mu), std::move(jacs));
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muGas(const ADB& pg,
                                       const ADB& T,
                                       const ADB& rv,
                                       const std::vector<PhasePresence>& cond,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.value().size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        enum PressureRvEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureRvEvalTag, /*size=*/2> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RvLad = 0.0;
        LadEval muLad;

        pLad.derivatives[0] = 1.0;
        RvLad.derivatives[1] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = pg.value()[i];
            TLad.value = T.value()[i];

            if (cond[i].hasFreeOil()) {
                muLad = gasPvt_->saturatedViscosity(pvtRegionIdx, TLad, pLad);
            }
            else {
                RvLad.value = rv.value()[i];
                muLad = gasPvt_->viscosity(pvtRegionIdx, TLad, pLad, RvLad);
            }

            mu[i] = muLad.value;
            dmudp[i] = muLad.derivatives[0];
            dmudr[i] = muLad.derivatives[1];
        }

        ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
        ADB::M dmudr_diag(dmudr.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dmudp_diag, pg.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dmudr_diag, rv.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(mu), std::move(jacs));
    }


    // ------ Formation volume factor (b) ------


    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bWat(const ADB& pw,
                                      const ADB& T,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call bWat(): water phase not active.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;

        pLad.derivatives[0] = 1.0;
        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = pw.value()[i];
            TLad.value = T.value()[i];

            const LadEval& bLad = waterPvt_->inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);

            b[i] = bLad.value;
            dbdp[i] = bLad.derivatives[0];
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, pw.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(b), std::move(jacs));
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bOil(const ADB& po,
                                      const ADB& T,
                                      const ADB& rs,
                                      const std::vector<PhasePresence>& cond,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call bOil(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        enum PressureRsEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureRsEvalTag, /*size=*/2> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RsLad = 0.0;
        LadEval bLad;

        pLad.derivatives[0] = 1.0;
        RsLad.derivatives[1] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = po.value()[i];
            TLad.value = T.value()[i];

            if (cond[i].hasFreeGas()) {
                bLad = oilPvt_->saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
            }
            else {
                RsLad.value = rs.value()[i];
                bLad = oilPvt_->inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RsLad);
            }

            b[i] = bLad.value;
            dbdp[i] = bLad.derivatives[0];
            dbdr[i] = bLad.derivatives[1];
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        ADB::M dbdr_diag(dbdr.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, po.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dbdr_diag, rs.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(b), std::move(jacs));
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bGas(const ADB& pg,
                                      const ADB& T,
                                      const ADB& rv,
                                      const std::vector<PhasePresence>& cond,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call bGas(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        enum PressureRvEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureRvEvalTag, /*size=*/2> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RvLad = 0.0;
        LadEval bLad;

        pLad.derivatives[0] = 1.0;
        RvLad.derivatives[1] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = pg.value()[i];
            TLad.value = T.value()[i];

            if (cond[i].hasFreeOil()) {
                bLad = gasPvt_->saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
            }
            else {
                RvLad.value = rv.value()[i];
                bLad = gasPvt_->inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RvLad);
            }

            b[i] = bLad.value;
            dbdp[i] = bLad.derivatives[0];
            dbdr[i] = bLad.derivatives[1];
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        ADB::M dbdr_diag(dbdr.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, pg.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dbdr_diag, rv.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(b), std::move(jacs));
    }



    // ------ Rs bubble point curve ------

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAdFromDeck::rsSat(const ADB& po,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call rsSat(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V rbub(n);
        V drbubdp(n);

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 293.15; // temperature is not supported by this API!

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = po.value()[i];

            const LadEval& RsLad = oilPvt_->saturatedGasDissolutionFactor(pvtRegionIdx, TLad, pLad);

            rbub[i] = RsLad.value;
            drbubdp[i] = RsLad.derivatives[0];
        }

        ADB::M drbubdp_diag(drbubdp.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(drbubdp_diag, po.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(rbub), std::move(jacs));
    }

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAdFromDeck::rsSat(const ADB& po,
                                       const ADB& so,
                                       const Cells& cells) const
    {
        ADB rs = rsSat(po, cells);
        applyVap(rs, so, cells, vap2_);
        return rs;
    }

    // ------ Rv condensation curve ------

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n condensation point values for Rv.
    ADB BlackoilPropsAdFromDeck::rvSat(const ADB& pg,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call rvSat(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.size() == n);
        V rv(n);
        V drvdp(n);

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 293.15; // temperature is not supported by this API!

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pLad.value = pg.value()[i];

            const LadEval& RvLad = gasPvt_->saturatedOilVaporizationFactor(pvtRegionIdx, TLad, pLad);

            rv[i] = RvLad.value;
            drvdp[i] = RvLad.derivatives[0];
        }

        ADB::M drvdp_diag(drvdp.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(drvdp_diag, pg.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(rv), std::move(jacs));
    }

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n condensation point values for Rv.
    ADB BlackoilPropsAdFromDeck::rvSat(const ADB& po,
                                       const ADB& so,
                                       const Cells& cells) const
    {
        ADB rv = rvSat(po, cells);
        applyVap(rv, so, cells, vap1_);
        return rv;
    }

    // ------ Relative permeability ------

    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  sg     Array of n gas saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 3 elements, each an array of n relperm values,
    ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
    std::vector<ADB> BlackoilPropsAdFromDeck::relperm(const ADB& sw,
                                                      const ADB& so,
                                                      const ADB& sg,
                                                      const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Water]) = sw.value();
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Oil]) = so.value();
        } else {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::relperm() assumes oil phase is active.");
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Gas]) = sg.value();
        }
        Block kr(n, np);
        Block dkr(n, np*np);
        satprops_->relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(3);
        typedef const ADB* ADBPtr;
        ADBPtr s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (phase_usage_.phase_used[phase1]) {
                const int phase1_pos = phase_usage_.phase_pos[phase1];
                std::vector<ADB::M> jacs(num_blocks);
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!phase_usage_.phase_used[phase2]) {
                        continue;
                    }
                    const int phase2_pos = phase_usage_.phase_pos[phase2];
                    // Assemble dkr1/ds2.
                    const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()

                    ADB::M dkr1_ds2_diag(dkr.col(column).matrix().asDiagonal());
                    for (int block = 0; block < num_blocks; ++block) {
                        ADB::M temp;
                        fastSparseProduct(dkr1_ds2_diag, s[phase2]->derivative()[block], temp);
                        jacs[block] += temp;
                    }
                }
                ADB::V val = kr.col(phase1_pos);
                relperms.emplace_back(ADB::function(std::move(val), std::move(jacs)));
            } else {
                relperms.emplace_back(ADB::null());
            }
        }
        return relperms;
    }

    std::vector<ADB> BlackoilPropsAdFromDeck::capPress(const ADB& sw,
                                                       const ADB& so,
                                                       const ADB& sg,
                                                       const Cells& cells) const
    {
        const int nCells = cells.size();
        const int nActivePhases = numPhases();
        const int nBlocks = so.numBlocks();

        Block activeSat(nCells, nActivePhases);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Water]) = sw.value();
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Oil]) = so.value();
        } else {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::relperm() assumes oil phase is active.");
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Gas]) = sg.value();
        }

        Block pc(nCells, nActivePhases);
        Block dpc(nCells, nActivePhases*nActivePhases);
        satprops_->capPress(nCells, activeSat.data(), cells.data(), pc.data(), dpc.data());

        std::vector<ADB> adbCapPressures;
        adbCapPressures.reserve(3);
        const ADB* s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (phase_usage_.phase_used[phase1]) {
                const int phase1_pos = phase_usage_.phase_pos[phase1];
                std::vector<ADB::M> jacs(nBlocks);
                for (int block = 0; block < nBlocks; ++block) {
                    jacs[block] = ADB::M(nCells, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!phase_usage_.phase_used[phase2])
                        continue;
                    const int phase2_pos = phase_usage_.phase_pos[phase2];
                    // Assemble dpc1/ds2.
                    const int column = phase1_pos + nActivePhases*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                    ADB::M dpc1_ds2_diag(dpc.col(column).matrix().asDiagonal());
                    for (int block = 0; block < nBlocks; ++block) {
                        ADB::M temp;
                        fastSparseProduct(dpc1_ds2_diag, s[phase2]->derivative()[block], temp);
                        jacs[block] += temp;
                    }
                }
                ADB::V val = pc.col(phase1_pos);
                adbCapPressures.emplace_back(ADB::function(std::move(val), std::move(jacs)));
            } else {
                adbCapPressures.emplace_back(ADB::null());
            }
        }
        return adbCapPressures;
    }

    /// Saturation update for hysteresis behavior.
    /// \param[in]  cells       Array of n cell indices to be associated with the saturation values.
    void BlackoilPropsAdFromDeck::updateSatHyst(const std::vector<double>& saturation,
                                                const std::vector<int>& cells)
    {
        const int n = cells.size();
        satprops_->updateSatHyst(n, cells.data(), saturation.data());
    }

    /// Update for max oil saturation.
    void BlackoilPropsAdFromDeck::updateSatOilMax(const std::vector<double>& saturation)
    {
        if (!satOilMax_.empty()) {
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


    /// Apply correction to rs/rv according to kw VAPPARS
    /// \param[in/out] r     Array of n rs/rv values.
    /// \param[in]     so    Array of n oil saturation values.
    /// \param[in]     cells Array of n cell indices to be associated with the r and so values.
    /// \param[in]     vap   Correction parameter.
    void BlackoilPropsAdFromDeck::applyVap(V& r,
                                           const V& so,
                                           const std::vector<int>& cells,
                                           const double vap) const
    {
        if (!satOilMax_.empty() && vap > 0.0) {
            const int n = cells.size();
            V factor = V::Ones(n, 1);
            const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
            for (int i=0; i<n; ++i) {
                if (satOilMax_[cells[i]] > vap_satmax_guard_ && so[i] < satOilMax_[cells[i]]) {
                    // guard against too small saturation values.
                    const double so_i= std::max(so[i],eps_sqrt);
                    factor[i] = std::pow(so_i/satOilMax_[cells[i]], vap);
                }
            }
            r = factor*r;
        }
    }

    /// Apply correction to rs/rv according to kw VAPPARS
    /// \param[in/out] r     Array of n rs/rv values.
    /// \param[in]     so    Array of n oil saturation values.
    /// \param[in]     cells Array of n cell indices to be associated with the r and so values.
    /// \param[in]     vap   Correction parameter.
    void BlackoilPropsAdFromDeck::applyVap(ADB& r,
                                           const ADB& so,
                                           const std::vector<int>& cells,
                                           const double vap) const
    {
        if (!satOilMax_.empty() && vap > 0.0) {
            const int n = cells.size();
            V factor = V::Ones(n, 1);
            const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
            V dfactor_dso = V::Zero(n, 1);
            for (int i=0; i<n; ++i) {
                if (satOilMax_[cells[i]] > vap_satmax_guard_ && so.value()[i] < satOilMax_[cells[i]]) {
                    // guard against too small saturation values.
                    const double so_i= std::max(so.value()[i],eps_sqrt);
                    factor[i] = std::pow(so_i/satOilMax_[cells[i]], vap);
                    dfactor_dso[i] = vap*std::pow(so_i/satOilMax_[cells[i]], vap-1.0)/satOilMax_[cells[i]];
                }
            }
            ADB::M dfactor_dso_diag(dfactor_dso.matrix().asDiagonal());
            const int num_blocks = so.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = dfactor_dso_diag * so.derivative()[block];
            }
            r = ADB::function(std::move(factor), std::move(jacs))*r;
        }
    }

    /// Obtain the scaled critical oil in gas saturation values.
    /// \param[in]  cells  Array of cell indices.
    /// \return Array of critical oil in gas saturaion values.
    V BlackoilPropsAdFromDeck::scaledCriticalOilinGasSaturations(const Cells& cells) const {


        assert(phase_usage_.phase_used[Gas]);
        assert(phase_usage_.phase_used[Oil]);

        const int n = cells.size();
        V sogcr = V::Zero(n);
        const MaterialLawManager& materialLawManager = satprops_->materialLawManager();
        for (int i = 0; i < n; ++i) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cells[i]);
                sogcr[i] = scaledDrainageInfo.Sogcr;
        }
        return sogcr;
    }

    /// Obtain the scaled critical gas saturation values.
    /// \param[in]  cells  Array of cell indices.
    /// \return Array of scaled critical gas saturaion values.
    V BlackoilPropsAdFromDeck::scaledCriticalGasSaturations(const Cells& cells) const {

        assert(phase_usage_.phase_used[Gas]);

        const int n = cells.size();
        V sgcr = V::Zero(n);
        const MaterialLawManager& materialLawManager = satprops_->materialLawManager();
        for (int i = 0; i < n; ++i) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cells[i]);
                sgcr[i] = scaledDrainageInfo.Sgcr;
        }
        return sgcr;
    }


} // namespace Opm

