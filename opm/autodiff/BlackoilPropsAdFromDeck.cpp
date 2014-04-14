/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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
#include <opm/core/props/pvt/SinglePvtInterface.hpp>
#include <opm/core/props/pvt/SinglePvtConstCompr.hpp>
#include <opm/core/props/pvt/SinglePvtDead.hpp>
#include <opm/core/props/pvt/SinglePvtDeadSpline.hpp>
#include <opm/core/props/pvt/SinglePvtLiveOil.hpp>
#include <opm/core/props/pvt/SinglePvtLiveGas.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Utility/PvdoTable.hpp>
#include <opm/parser/eclipse/Utility/PvtgTable.hpp>
#include <opm/parser/eclipse/Utility/PvcdoTable.hpp>

namespace Opm
{
    // Making these typedef to make the code more readable.
    typedef BlackoilPropsAdFromDeck::ADB ADB;
    typedef BlackoilPropsAdFromDeck::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;
    enum { Aqua = BlackoilPhases::Aqua,
           Liquid = BlackoilPhases::Liquid,
           Vapour = BlackoilPhases::Vapour };

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const EclipseGridParser& deck,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        init(deck, grid.number_of_cells, grid.global_cell, grid.cartdims, 
             grid.cell_centroids, grid.dimensions, init_rock);
    }

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr newParserDeck,
                                                     const UnstructuredGrid& grid,
                                                     const bool init_rock)
    {
        init(newParserDeck, grid.number_of_cells, grid.global_cell, grid.cartdims, 
             grid.cell_centroids, grid.dimensions, init_rock);
    }

#ifdef HAVE_DUNE_CORNERPOINT
        
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(const EclipseGridParser& deck,
                                                     const Dune::CpGrid& grid,
                                                     const bool init_rock )
    {
        init(deck, grid.numCells(), static_cast<const int*>(&grid.globalCell()[0]),
             static_cast<const int*>(&grid.logicalCartesianSize()[0]),
             grid.beginCellCentroids(), Dune::CpGrid::dimension, init_rock);
    }

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr newParserDeck,
                                                     const Dune::CpGrid& grid,
                                                     const bool init_rock )
    {
        init(newParserDeck, grid.numCells(), static_cast<const int*>(&grid.globalCell()[0]),
             static_cast<const int*>(&grid.logicalCartesianSize()[0]),
             grid.beginCellCentroids(), Dune::CpGrid::dimension, init_rock);
    }
#endif

    /// Constructor wrapping an opm-core black oil interface.
    template<class T>
    BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck(Opm::DeckConstPtr newParserDeck,
                                                     int number_of_cells,
                                                     const int* global_cell,
                                                     const int* cart_dims,
                                                     T begin_cell_centroids,
                                                     int dimensions,
                                                     const bool init_rock)
    {
        init(newParserDeck, number_of_cells, global_cell, cart_dims, begin_cell_centroids,
             dimensions, init_rock);
    }

    template<class T>
    void BlackoilPropsAdFromDeck::init(const EclipseGridParser& deck,
                                       int number_of_cells,
                                       const int* global_cell,
                                       const int* cart_dims,
                                       T begin_cell_centroids,
                                       int dimensions,
                                       const bool init_rock)
    {
        if (init_rock){
            rock_.init(deck, number_of_cells, global_cell, cart_dims);
        }
        const int samples = 0;
        const int region_number = 0;

        phase_usage_ = phaseUsageFromDeck(deck);

        // Surface densities. Accounting for different orders in eclipse and our code.
        if (deck.hasField("DENSITY")) {
            const std::vector<double>& d = deck.getDENSITY().densities_[region_number];
            enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
            if (phase_usage_.phase_used[Aqua]) {
                densities_[phase_usage_.phase_pos[Aqua]]   = d[ECL_water];
            }
            if (phase_usage_.phase_used[Vapour]) {
                densities_[phase_usage_.phase_pos[Vapour]] = d[ECL_gas];
            }
            if (phase_usage_.phase_used[Liquid]) {
                densities_[phase_usage_.phase_pos[Liquid]] = d[ECL_oil];
            }
        } else {
            OPM_THROW(std::runtime_error, "Input is missing DENSITY\n");
        }

        // Set the properties.
        props_.resize(phase_usage_.num_phases);
        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            if (deck.hasField("PVTW")) {
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(deck.getPVTW().pvtw_));
            } else {
                // Eclipse 100 default.
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(0.5*Opm::prefix::centi*Opm::unit::Poise));
            }
        }
        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            if (deck.hasField("PVDO")) {
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDeadSpline(deck.getPVDO().pvdo_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDead(deck.getPVDO().pvdo_));
                }
            } else if (deck.hasField("PVTO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtLiveOil(deck.getPVTO().pvto_));
            } else if (deck.hasField("PVCDO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtConstCompr(deck.getPVCDO().pvcdo_));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDO, PVTO or PVCDO\n");
            }
        }
        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            if (deck.hasField("PVDG")) {
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDeadSpline(deck.getPVDG().pvdg_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDead(deck.getPVDG().pvdg_));
                }
             } else if (deck.hasField("PVTG")) {
                 props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtLiveGas(deck.getPVTG().pvtg_));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDG or PVTG\n");
            }
        }

        SaturationPropsFromDeck<SatFuncGwsegNonuniform>* ptr
            = new SaturationPropsFromDeck<SatFuncGwsegNonuniform>();
        satprops_.reset(ptr);
        ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimensions, -1);

        if (phase_usage_.num_phases != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck() - "
                  "Inconsistent number of phases in pvt data (" << phase_usage_.num_phases
                  << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class T>
    void BlackoilPropsAdFromDeck::init(Opm::DeckConstPtr newParserDeck,
                                       int number_of_cells,
                                       const int* global_cell,
                                       const int* cart_dims,
                                       T begin_cell_centroids,
                                       int dimensions,
                                       const bool init_rock)
    {
        if (init_rock){
            rock_.init(newParserDeck, number_of_cells, global_cell, cart_dims);
        }

        const int samples = 0;
        const int region_number = 0;

        phase_usage_ = phaseUsageFromDeck(newParserDeck);

        // Surface densities. Accounting for different orders in eclipse and our code.
        if (newParserDeck->hasKeyword("DENSITY")) {
            const auto keyword = newParserDeck->getKeyword("DENSITY");
            const auto record = keyword->getRecord(region_number);
            enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
            if (phase_usage_.phase_used[Aqua]) {
                densities_[phase_usage_.phase_pos[Aqua]]   = record->getItem("WATER")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Vapour]) {
                densities_[phase_usage_.phase_pos[Vapour]] = record->getItem("GAS")->getSIDouble(0);
            }
            if (phase_usage_.phase_used[Liquid]) {
                densities_[phase_usage_.phase_pos[Liquid]] = record->getItem("OIL")->getSIDouble(0);
            }
        } else {
            OPM_THROW(std::runtime_error, "Input is missing DENSITY\n");
        }

        // Set the properties.
        props_.resize(phase_usage_.num_phases);
        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            if (newParserDeck->hasKeyword("PVTW")) {
                Opm::PvtwTable pvtwTable(newParserDeck->getKeyword("PVTW"), region_number);
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(pvtwTable));
            } else {
                // Eclipse 100 default.
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(0.5*Opm::prefix::centi*Opm::unit::Poise));
            }
        }

        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            if (newParserDeck->hasKeyword("PVDO")) {
                Opm::PvdoTable pvdoTable(newParserDeck->getKeyword("PVDO"), region_number);
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDeadSpline(pvdoTable, samples));
                } else {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDead(pvdoTable));
                }
            }
            else if (newParserDeck->hasKeyword("PVTO")) {
                Opm::PvtoTable pvtoTable(newParserDeck->getKeyword("PVTO"), /*tableIdx=*/0);
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtLiveOil(pvtoTable));
            } else if (newParserDeck->hasKeyword("PVCDO")) {
                Opm::PvcdoTable pvdcoTable(newParserDeck->getKeyword("PVCDO"), region_number);
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtConstCompr(pvdcoTable));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDO, PVTO or PVCDO\n");
            }
        }

        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            if (newParserDeck->hasKeyword("PVDG")) {
                Opm::PvdoTable pvdgTable(newParserDeck->getKeyword("PVDG"), region_number);
                if (samples > 0) {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDeadSpline(pvdgTable, samples));
                } else {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDead(pvdgTable));
                }
            } else if (newParserDeck->hasKeyword("PVTG")) {
                Opm::PvtgTable pvtgTable(newParserDeck->getKeyword("PVTG"), /*tableIdx=*/0);
                props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtLiveGas(pvtgTable));
            } else {
                OPM_THROW(std::runtime_error, "Input is missing PVDG or PVTG\n");
            }
        }

        SaturationPropsFromDeck<SatFuncGwsegNonuniform>* ptr
            = new SaturationPropsFromDeck<SatFuncGwsegNonuniform>();
        satprops_.reset(ptr);
        ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids, dimensions, -1);

        if (phase_usage_.num_phases != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck() - "
                  "Inconsistent number of phases in pvt data (" << phase_usage_.num_phases
                  << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
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
    const double* BlackoilPropsAdFromDeck::surfaceDensity() const
    {
        return densities_;
    }


    // ------ Viscosity ------

    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAdFromDeck::muWat(const V& pw,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        assert(pw.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->mu(n, pw.data(), rs,
                                                  mu.data(), dmudp.data(), dmudr.data());
        return mu;
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAdFromDeck::muOil(const V& po,
                                     const V& rs,
                                     const std::vector<PhasePresence>& cond,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Oil]]->mu(n, po.data(), rs.data(), &cond[0],
                                                mu.data(), dmudp.data(), dmudr.data());
        return mu;
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAdFromDeck::muGas(const V& pg,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Gas]]->mu(n, pg.data(), rs,
                                                mu.data(), dmudp.data(), dmudr.data());
        return mu;
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAdFromDeck::muGas(const V& pg,
                                     const V& rv,
                                     const std::vector<PhasePresence>& cond,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Gas]]->mu(n, pg.data(), rv.data(),&cond[0],
                                                mu.data(), dmudp.data(), dmudr.data());
        return mu;
    }

    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muWat(const ADB& pw,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        assert(pw.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->mu(n, pw.value().data(), rs,
                                                  mu.data(), dmudp.data(), dmudr.data());
        ADB::M dmudp_diag = spdiag(dmudp);
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmudp_diag * pw.derivative()[block];
        }
        return ADB::function(mu, jacs);
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muOil(const ADB& po,
                                       const ADB& rs,
                                       const std::vector<PhasePresence>& cond,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Oil]]->mu(n, po.value().data(), rs.value().data(),
                                                &cond[0], mu.data(), dmudp.data(), dmudr.data());

        ADB::M dmudp_diag = spdiag(dmudp);
        ADB::M dmudr_diag = spdiag(dmudr);
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmudp_diag * po.derivative()[block] + dmudr_diag * rs.derivative()[block];
        }
        return ADB::function(mu, jacs);
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muGas(const ADB& pg,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.value().size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);
        const double* rv = 0;

        props_[phase_usage_.phase_pos[Gas]]->mu(n, pg.value().data(), rv,
                                                  mu.data(), dmudp.data(), dmudr.data());

        ADB::M dmudp_diag = spdiag(dmudp);
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmudp_diag * pg.derivative()[block];
        }
        return ADB::function(mu, jacs);
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeck::muGas(const ADB& pg,
                                       const ADB& rv,
                                       const std::vector<PhasePresence>& cond,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.value().size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        props_[phase_usage_.phase_pos[Gas]]->mu(n, pg.value().data(), rv.value().data(),&cond[0],
                                                  mu.data(), dmudp.data(), dmudr.data());

        ADB::M dmudp_diag = spdiag(dmudp);
        ADB::M dmudr_diag = spdiag(dmudr);
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmudp_diag * pg.derivative()[block] + dmudr_diag * rv.derivative()[block];
        }
        return ADB::function(mu, jacs);
    }


    // ------ Formation volume factor (b) ------

    // These methods all call the matrix() method, after which the variable
    // (also) called 'matrix' contains, in each row, the A = RB^{-1} matrix for
    // a cell. For three-phase black oil:
    //  A = [  bw       0       0
    //          0       bo      0
    //          0      b0*rs   bw ]
    // Where b = B^{-1}.
    // Therefore, we extract the correct diagonal element, and are done.
    // When we need the derivatives (w.r.t. p, since we don't do w.r.t. rs),
    // we also get the following derivative matrix:
    //  A = [  dbw       0       0
    //          0       dbo      0
    //          0      db0*rs   dbw ]
    // Again, we just extract a diagonal element.

    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAdFromDeck::bWat(const V& pw,
                                    const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call bWat(): water phase not present.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->b(n, pw.data(), rs,
                                                 b.data(), dbdp.data(), dbdr.data());

        return b;
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAdFromDeck::bOil(const V& po,
                                    const V& rs,
                                    const std::vector<PhasePresence>& cond,
                                    const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call bOil(): oil phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Oil]]->b(n, po.data(), rs.data(), &cond[0],
                                               b.data(), dbdp.data(), dbdr.data());

        return b;
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAdFromDeck::bGas(const V& pg,
                                    const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call bGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Gas]]->b(n, pg.data(), rs,
                                               b.data(), dbdp.data(), dbdr.data());

        return b;
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAdFromDeck::bGas(const V& pg,
           const V& rv,
           const std::vector<PhasePresence>& cond,
           const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Gas]]->b(n, pg.data(), rv.data(), &cond[0],
                                               b.data(), dbdp.data(), dbdr.data());

        return b;
    }

    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bWat(const ADB& pw,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);
        const double* rs = 0;

        props_[phase_usage_.phase_pos[Water]]->b(n, pw.value().data(), rs,
                                                 b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dbdp_diag * pw.derivative()[block];
        }
        return ADB::function(b, jacs);
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bOil(const ADB& po,
                                      const ADB& rs,
                                      const std::vector<PhasePresence>& cond,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Oil]]->b(n, po.value().data(), rs.value().data(),
                                               &cond[0], b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        ADB::M dbdr_diag = spdiag(dbdr);
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dbdp_diag * po.derivative()[block] + dbdr_diag * rs.derivative()[block];
        }
        return ADB::function(b, jacs);
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bGas(const ADB& pg,
                                      const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);
        const double* rv = 0;

        props_[phase_usage_.phase_pos[Gas]]->b(n, pg.value().data(), rv,
                                               b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dbdp_diag * pg.derivative()[block];
        }
        return ADB::function(b, jacs);
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeck::bGas(const ADB& pg,
           const ADB& rv,
           const std::vector<PhasePresence>& cond,
           const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        props_[phase_usage_.phase_pos[Gas]]->b(n, pg.value().data(), rv.value().data(), &cond[0],
                                               b.data(), dbdp.data(), dbdr.data());

        ADB::M dbdp_diag = spdiag(dbdp);
        ADB::M dmudr_diag = spdiag(dbdr);
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dbdp_diag * pg.derivative()[block] + dmudr_diag * rv.derivative()[block];;
        }
        return ADB::function(b, jacs);
    }



    // ------ Rs bubble point curve ------

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    V BlackoilPropsAdFromDeck::rsSat(const V& po,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call rsMax(): oil phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V rbub(n);
        V drbubdp(n);
        props_[Oil]->rsSat(n, po.data(), rbub.data(), drbubdp.data());
        return rbub;
    }

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
        assert(po.size() == n);
        V rbub(n);
        V drbubdp(n);
        props_[Oil]->rsSat(n, po.value().data(), rbub.data(), drbubdp.data());
        ADB::M drbubdp_diag = spdiag(drbubdp);
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = drbubdp_diag * po.derivative()[block];
        }
        return ADB::function(rbub, jacs);
    }

    // ------ Condensation curve ------

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    V BlackoilPropsAdFromDeck::rvSat(const V& po,
                                     const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call rvMax(): gas phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V rv(n);
        V drvdp(n);
        props_[Gas]->rvSat(n, po.data(), rv.data(), drvdp.data());
        return rv;
    }

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAdFromDeck::rvSat(const ADB& po,
                                       const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call rvMax(): gas phase not present.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V rv(n);
        V drvdp(n);
        props_[Gas]->rvSat(n, po.value().data(), rv.data(), drvdp.data());
        ADB::M drvdp_diag = spdiag(drvdp);
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = drvdp_diag * po.derivative()[block];
        }
        return ADB::function(rv, jacs);
    }

    // ------ Relative permeability ------

    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  sg     Array of n gas saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 3 elements, each an array of n relperm values,
    ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
    std::vector<V> BlackoilPropsAdFromDeck::relperm(const V& sw,
                                                    const V& so,
                                                    const V& sg,
                                                    const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.size() == n);
            s_all.col(phase_usage_.phase_pos[Water]) = sw;
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.size() == n);
            s_all.col(phase_usage_.phase_pos[Oil]) = so;
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.size() == n);
            s_all.col(phase_usage_.phase_pos[Gas]) = sg;
        }
        Block kr(n, np);
        satprops_->relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        std::vector<V> relperms;
        relperms.reserve(3);
        for (int phase = 0; phase < 3; ++phase) {
            if (phase_usage_.phase_used[phase]) {
                relperms.emplace_back(kr.col(phase_usage_.phase_pos[phase]));
            } else {
                relperms.emplace_back();
            }
        }
        return relperms;
    }

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
                        jacs[block] += dkr1_ds2_diag * s[phase2]->derivative()[block];
                    }
                }
                relperms.emplace_back(ADB::function(kr.col(phase1_pos), jacs));
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
                        jacs[block] += dpc1_ds2_diag * s[phase2]->derivative()[block];
                    }
                }
                adbCapPressures.emplace_back(ADB::function(pc.col(phase1_pos), jacs));
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

} // namespace Opm

