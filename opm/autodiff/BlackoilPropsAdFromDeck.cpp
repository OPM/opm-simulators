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
#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/core/props/pvt/PvtConstCompr.hpp>
#include <opm/core/props/pvt/PvtDead.hpp>
#include <opm/core/props/pvt/PvtDeadSpline.hpp>
#include <opm/core/props/pvt/PvtLiveOil.hpp>
#include <opm/core/props/pvt/PvtLiveGas.hpp>
#include <opm/core/props/pvt/ThermalWaterPvtWrapper.hpp>
#include <opm/core/props/pvt/ThermalOilPvtWrapper.hpp>
#include <opm/core/props/pvt/ThermalGasPvtWrapper.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

namespace Opm
{
    // Making these typedef to make the code more readable.
    typedef BlackoilPropsAdFromDeck::ADB ADB;
    typedef BlackoilPropsAdFromDeck::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        init(deck, eclState, grid.number_of_cells, grid.global_cell, grid.cartdims, 
             grid.cell_centroids, grid.dimensions, init_rock);
    }

#ifdef HAVE_DUNE_CORNERPOINT
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr deck,
                                                     Opm::EclipseStateConstPtr eclState,
                                                     const Dune::CpGrid& grid,
                                                     const bool init_rock )
    {
        init(deck, eclState, grid.numCells(), static_cast<const int*>(&grid.globalCell()[0]),
             static_cast<const int*>(&grid.logicalCartesianSize()[0]),
             grid.beginCellCentroids(), Dune::CpGrid::dimension, init_rock);
    }
#endif

/// Constructor for properties on a subgrid
BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const BlackoilPropsAdFromDeck& props,
                                                 const int number_of_cells)
{
    const int original_size = props.cellPvtRegionIdx_.size();
    if (number_of_cells > original_size) {
        OPM_THROW(std::runtime_error, "The number of cells is larger than the one of the original grid!");
    }
    if (number_of_cells < 0) {
        OPM_THROW(std::runtime_error, "The number of cells is has to be larger than 0.");
    }
    // Copy properties that do not depend on the postion within the grid.
    rock_             = props.rock_;
    satprops_         = props.satprops_;
    phase_usage_      = props.phase_usage_;
    props_            = props.props_;
    densities_        = props.densities_;
    vap1_             = props.vap1_;
    vap2_             = props.vap2_;
    vap_satmax_guard_ = props.vap_satmax_guard_;
    // For data that is dependant on the subgrid we simply allocate space
    // and initialize with obviously bogus numbers.
    cellPvtRegionIdx_.resize(number_of_cells, std::numeric_limits<int>::min());
}

    /// Initializes the properties.
    template <class CentroidIterator>
    void BlackoilPropsAdFromDeck::init(Opm::DeckConstPtr deck,
                                       Opm::EclipseStateConstPtr eclState,
                                       int number_of_cells,
                                       const int* global_cell,
                                       const int* cart_dims,
                                       const CentroidIterator& begin_cell_centroids,
                                       int dimension,
                                       const bool init_rock)
    {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, deck, number_of_cells, global_cell);

        if (init_rock){
            rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }

        phase_usage_ = phaseUsageFromDeck(deck);



        // Surface densities. Accounting for different orders in eclipse and our code.
        Opm::DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");
        int numRegions = densityKeyword->size();

        densities_.resize(numRegions);
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            if (phase_usage_.phase_used[Liquid]) {
                densities_[regionIdx][phase_usage_.phase_pos[Liquid]]
                    = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Aqua]) {
                densities_[regionIdx][phase_usage_.phase_pos[Aqua]]
                    = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Vapour]) {
                densities_[regionIdx][phase_usage_.phase_pos[Vapour]]
                    = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            }
        }

        const int numSamples = 0;

        // Resize the property objects container
        props_.resize(phase_usage_.num_phases);


        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            // if water is used, we require the presence of the "PVTW"
            // keyword for now...
            std::shared_ptr<PvtConstCompr> pvtw(new PvtConstCompr);
            pvtw->initFromWater(deck->getKeyword("PVTW"));
            props_[phase_usage_.phase_pos[Aqua]] = pvtw;

            // handle temperature dependence of the oil phase
            if (!eclState->getWatvisctTables().empty() || deck->hasKeyword("WATDENT")) {
                // deal with temperature dependent properties
                std::shared_ptr<ThermalWaterPvtWrapper> waterNiPvt(new ThermalWaterPvtWrapper);
                waterNiPvt->initFromDeck(props_[phase_usage_.phase_pos[Aqua]], deck, eclState);

                props_[phase_usage_.phase_pos[Aqua]] = waterNiPvt;
            }
        }
        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            // for oil, we support the "PVDO", "PVTO" and "PVCDO"
            // keywords...
            const auto& pvdoTables = eclState->getPvdoTables();
            const auto& pvtoTables = eclState->getPvtoTables();
            if (!pvdoTables.empty()) {
                if (numSamples > 0) {
                    auto splinePvdo = std::shared_ptr<PvtDeadSpline>(new PvtDeadSpline);
                    splinePvdo->initFromOil(pvdoTables, numSamples);
                    props_[phase_usage_.phase_pos[Liquid]] = splinePvdo;
                } else {
                    auto pvdo = std::shared_ptr<PvtDead>(new PvtDead);
                    pvdo->initFromOil(pvdoTables);
                    props_[phase_usage_.phase_pos[Liquid]] = pvdo;
                }
            } else if (!pvtoTables.empty()) {
                std::shared_ptr<PvtLiveOil> pvto(new PvtLiveOil(pvtoTables));
                props_[phase_usage_.phase_pos[Liquid]] = pvto;
            } else if (deck->hasKeyword("PVCDO")) {
                std::shared_ptr<PvtConstCompr> pvcdo(new PvtConstCompr);
                pvcdo->initFromOil(deck->getKeyword("PVCDO"));
                props_[phase_usage_.phase_pos[Liquid]] = pvcdo;
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDO, PVCDO or PVTO\n");
            }

            // handle temperature dependence of the oil phase
            if (!eclState->getOilvisctTables().empty() || deck->hasKeyword("THERMEX1")) {
                std::shared_ptr<ThermalOilPvtWrapper> oilNiPvt(new ThermalOilPvtWrapper);
                oilNiPvt->initFromDeck(props_[phase_usage_.phase_pos[Liquid]], deck, eclState);

                props_[phase_usage_.phase_pos[Liquid]] = oilNiPvt;
            }
        }
        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            // gas can be specified using the "PVDG" or "PVTG" keywords...
            const auto& pvdgTables = eclState->getPvdgTables();
            const auto& pvtgTables = eclState->getPvtgTables();
            if (!pvdgTables.empty()) {
                if (numSamples > 0) {
                    std::shared_ptr<PvtDeadSpline> splinePvt(new PvtDeadSpline);
                    splinePvt->initFromGas(pvdgTables, numSamples);
                    props_[phase_usage_.phase_pos[Vapour]] = splinePvt;
                } else {
                    std::shared_ptr<PvtDead> deadPvt(new PvtDead);
                    deadPvt->initFromGas(pvdgTables);
                    props_[phase_usage_.phase_pos[Vapour]] = deadPvt;
                }
            } else if (!pvtgTables.empty()) {
                props_[phase_usage_.phase_pos[Vapour]].reset(new PvtLiveGas(pvtgTables));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDG or PVTG\n");
            }

            // handle temperature dependence of the gas phase
            if (!eclState->getGasvisctTables().empty() || deck->hasKeyword("TREF")) {
                std::shared_ptr<ThermalGasPvtWrapper> gasNiPvt(new ThermalGasPvtWrapper);
                gasNiPvt->initFromDeck(props_[phase_usage_.phase_pos[Vapour]], deck, eclState);

                props_[phase_usage_.phase_pos[Vapour]] = gasNiPvt;
            }
        }
        // Oil vaporization controls (kw VAPPARS)
        vap1_ = vap2_ = 0.0;
        if (deck->hasKeyword("VAPPARS") && deck->hasKeyword("VAPOIL") && deck->hasKeyword("DISGAS")) {
            vap1_ = deck->getKeyword("VAPPARS")->getRecord(0)->getItem(0)->getRawDouble(0);
            vap2_ = deck->getKeyword("VAPPARS")->getRecord(0)->getItem(1)->getRawDouble(0);
            satOilMax_.resize(number_of_cells, 0.0);
        } else if (deck->hasKeyword("VAPPARS")) {
            OPM_THROW(std::runtime_error, "Input has VAPPARS, but missing VAPOIL and/or DISGAS\n");
        }

        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        satprops_.reset(ptr);
        ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids, dimension);

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
    /// \return Array of 3 density values.
    const double* BlackoilPropsAdFromDeck::surfaceDensity(const int cellIdx) const
    {
        int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
        return &densities_[pvtRegionIdx][0];
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
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(pw.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->mu(n, pvt_region_.data(), pw.value().data(), T.value().data(), rs,
                                                  mu.data(), dmudp.data(), dmudr.data());
        if (pw.derivative().empty()) {
            return ADB::constant(std::move(mu));
        } else {
            ADB::M dmudp_diag = spdiag(dmudp);
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
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(po.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Oil]]->mu(n, pvt_region_.data(), po.value().data(), T.value().data(), rs.value().data(),
                                                &cond[0], mu.data(), dmudp.data(), dmudr.data());

        ADB::M dmudp_diag = spdiag(dmudp);
        ADB::M dmudr_diag = spdiag(dmudr);
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
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(pg.value().size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Gas]]->mu(n, pvt_region_.data(), pg.value().data(), T.value().data(), rv.value().data(),&cond[0],
                                                  mu.data(), dmudp.data(), dmudr.data());

        ADB::M dmudp_diag = spdiag(dmudp);
        ADB::M dmudr_diag = spdiag(dmudr);
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
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(pw.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->b(n, pvt_region_.data(), pw.value().data(), T.value().data(), rs,
                                                 b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
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
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(po.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Oil]]->b(n, pvt_region_.data(), po.value().data(), T.value().data(), rs.value().data(),
                                               &cond[0], b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        ADB::M dbdr_diag = spdiag(dbdr);
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
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Gas]]->b(n, pvt_region_.data(), pg.value().data(), T.value().data(), rv.value().data(), &cond[0],
                                               b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        ADB::M dbdr_diag = spdiag(dbdr);
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
            OPM_THROW(std::runtime_error, "Cannot call rsMax(): oil phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(po.size() == n);
        V rbub(n);
        V drbubdp(n);
        props_[phase_usage_.phase_pos[Oil]]->rsSat(n, pvt_region_.data(), po.value().data(), rbub.data(), drbubdp.data());
        ADB::M drbubdp_diag = spdiag(drbubdp);
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
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n condensation point values for Rv.
    ADB BlackoilPropsAdFromDeck::rvSat(const ADB& po,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call rvMax(): gas phase not present.");
        }
        const int n = cells.size();
        mapPvtRegions(cells);
        assert(po.size() == n);
        V rv(n);
        V drvdp(n);
        props_[phase_usage_.phase_pos[Gas]]->rvSat(n, pvt_region_.data(), po.value().data(), rv.data(), drvdp.data());
        ADB::M drvdp_diag = spdiag(drvdp);
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(drvdp_diag, po.derivative()[block], jacs[block]);
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
                    ADB::M dkr1_ds2_diag = spdiag(dkr.col(column));
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
        const int numCells = cells.size();
        const int numActivePhases = numPhases();
        const int numBlocks = so.numBlocks();

        Block activeSat(numCells, numActivePhases);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.value().size() == numCells);
            activeSat.col(phase_usage_.phase_pos[Water]) = sw.value();
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.value().size() == numCells);
            activeSat.col(phase_usage_.phase_pos[Oil]) = so.value();
        } else {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::relperm() assumes oil phase is active.");
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.value().size() == numCells);
            activeSat.col(phase_usage_.phase_pos[Gas]) = sg.value();
        }

        Block pc(numCells, numActivePhases);
        Block dpc(numCells, numActivePhases*numActivePhases);
        satprops_->capPress(numCells, activeSat.data(), cells.data(), pc.data(), dpc.data());

        std::vector<ADB> adbCapPressures;
        adbCapPressures.reserve(3);
        const ADB* s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (phase_usage_.phase_used[phase1]) {
                const int phase1_pos = phase_usage_.phase_pos[phase1];
                std::vector<ADB::M> jacs(numBlocks);
                for (int block = 0; block < numBlocks; ++block) {
                    jacs[block] = ADB::M(numCells, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!phase_usage_.phase_used[phase2])
                        continue;
                    const int phase2_pos = phase_usage_.phase_pos[phase2];
                    // Assemble dpc1/ds2.
                    const int column = phase1_pos + numActivePhases*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                    ADB::M dpc1_ds2_diag = spdiag(dpc.col(column));
                    for (int block = 0; block < numBlocks; ++block) {
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
            ADB::M dfactor_dso_diag = spdiag(dfactor_dso);
            const int num_blocks = so.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = dfactor_dso_diag * so.derivative()[block];
            }
            r = ADB::function(std::move(factor), std::move(jacs))*r;
        }
    }





    // Fills pvt_region_ with cellPvtRegionIdx_[cells].
    void BlackoilPropsAdFromDeck::mapPvtRegions(const std::vector<int>& cells) const
    {
        const int n = cells.size();
        pvt_region_.resize(n);
        for (int ii = 0; ii < n; ++ii) {
            pvt_region_[ii] = cellPvtRegionIdx_[cells[ii]];
        }
    }


} // namespace Opm

