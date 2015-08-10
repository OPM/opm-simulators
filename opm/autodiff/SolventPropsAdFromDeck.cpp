/*
  Copyright 2015 IRIS
  
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

#include "SolventPropsAdFromDeck.hpp"
#include <opm/autodiff/AutoDiffHelpers.hpp>


namespace Opm
{
// Making these typedef to make the code more readable.
typedef SolventPropsAdFromDeck::ADB ADB;
typedef SolventPropsAdFromDeck::V V;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

SolventPropsAdFromDeck::SolventPropsAdFromDeck(DeckConstPtr deck,
                                                     EclipseStateConstPtr eclState,
                                                     const int number_of_cells,
                                                     const int* global_cell)
{
    // retrieve the cell specific PVT table index from the deck
    // and using the grid...
    extractPvtTableIndex(cellPvtRegionIdx_, deck, number_of_cells, global_cell);

    // surface densities
    if (deck->hasKeyword("SDENSITY")) {
        Opm::DeckKeywordConstPtr densityKeyword = deck->getKeyword("SDENSITY");
        int numRegions = densityKeyword->size();
        solvent_surface_densities_.resize(numRegions);
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            solvent_surface_densities_[regionIdx]
                    = densityKeyword->getRecord(regionIdx)->getItem("SOLVENT_DENSITY")->getSIDouble(0);
        }
    } else {
        OPM_THROW(std::runtime_error, "SDENSITY must be specified in SOLVENT runs\n");
    }

    // pvt
    const auto& pvdsTables = eclState->getPvdsTables();
    if (!pvdsTables.empty()) {

        int numRegions = pvdsTables.size();
        // resize the attributes of the object
        b_.resize(numRegions);
        viscosity_.resize(numRegions);
        inverseBmu_.resize(numRegions);

        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const Opm::PvdsTable& pvdsTable = pvdsTables[regionIdx];

            // Copy data
            const std::vector<double>& press = pvdsTable.getPressureColumn();
            const std::vector<double>& b = pvdsTable.getFormationFactorColumn();
            const std::vector<double>& visc = pvdsTable.getViscosityColumn();

            const int sz = b.size();
            std::vector<double> inverseB(sz);
            for (int i = 0; i < sz; ++i) {
                inverseB[i] = 1.0 / b[i];
            }

            std::vector<double> inverseBmu(sz);
            for (int i = 0; i < sz; ++i) {
                inverseBmu[i] = 1.0 / (b[i] * visc[i]);
            }

            b_[regionIdx] = NonuniformTableLinear<double>(press, inverseB);
            viscosity_[regionIdx] = NonuniformTableLinear<double>(press, visc);
            inverseBmu_[regionIdx] = NonuniformTableLinear<double>(press, inverseBmu);
        }
    } else {
        OPM_THROW(std::runtime_error, "PVDS must be specified in SOLVENT runs\n");
    }

    const auto& ssfnTables = eclState->getSsfnTables();
    // relative permeabilty multiplier
    if (!ssfnTables.empty()) {

        int numRegions = pvdsTables.size();

        if(numRegions > 1) {
            OPM_THROW(std::runtime_error, "Only single table saturation function supported for SSFN");
        }
        // resize the attributes of the object
        krg_.resize(numRegions);
        krs_.resize(numRegions);
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            const Opm::SsfnTable& ssfnTable = ssfnTables[regionIdx];

            // Copy data
            const std::vector<double>& solventFraction = ssfnTable.getSolventFractionColumn();
            const std::vector<double>& krg = ssfnTable.getGasRelPermMultiplierColumn();
            const std::vector<double>& krs = ssfnTable.getSolventRelPermMultiplierColumn();

            krg_[regionIdx] = NonuniformTableLinear<double>(solventFraction, krg);
            krs_[regionIdx] = NonuniformTableLinear<double>(solventFraction, krs);
        }

    } else {
        OPM_THROW(std::runtime_error, "SSFN must be specified in SOLVENT runs\n");
    }
}

ADB SolventPropsAdFromDeck::muSolvent(const ADB& pg,
                                 const Cells& cells) const
{
    const int n = cells.size();
    assert(pg.value().size() == n);
    V mu(n);
    V dmudp(n);
    for (int i = 0; i < n; ++i) {
        const double& pg_i = pg.value()[i];
        int regionIdx = cellPvtRegionIdx_[i];
        double tempInvB = b_[regionIdx](pg_i);
        double tempInvBmu = inverseBmu_[regionIdx](pg_i);
        mu[i] = tempInvB / tempInvBmu;
        dmudp[i] = (tempInvBmu * b_[regionIdx].derivative(pg_i)
                         - tempInvB * inverseBmu_[regionIdx].derivative(pg_i)) / (tempInvBmu * tempInvBmu);
    }

    ADB::M dmudp_diag = spdiag(dmudp);
    const int num_blocks = pg.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        fastSparseProduct(dmudp_diag, pg.derivative()[block], jacs[block]);
    }
    return ADB::function(std::move(mu), std::move(jacs));
}

ADB SolventPropsAdFromDeck::bSolvent(const ADB& pg,
                                const Cells& cells) const
{
    const int n = cells.size();
    assert(pg.size() == n);

    V b(n);
    V dbdp(n);
    for (int i = 0; i < n; ++i) {
        const double& pg_i = pg.value()[i];
        int regionIdx = cellPvtRegionIdx_[i];
        b[i] = b_[regionIdx](pg_i);
        dbdp[i] = b_[regionIdx].derivative(pg_i);
    }

    ADB::M dbdp_diag = spdiag(dbdp);
    const int num_blocks = pg.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        fastSparseProduct(dbdp_diag, pg.derivative()[block], jacs[block]);
    }
    return ADB::function(std::move(b), std::move(jacs));
}

ADB SolventPropsAdFromDeck::gasRelPermMultiplier(const ADB& solventFraction,
                                 const Cells& cells) const
{
    const int n = cells.size();
    assert(solventFraction.value().size() == n);
    V krg(n);
    V dkrgdsf(n);
    for (int i = 0; i < n; ++i) {
        const double& solventFraction_i = solventFraction.value()[i];
        int regionIdx = 0; // TODO add mapping from cells to sat function table
        krg[i] = krg_[regionIdx](solventFraction_i);
        dkrgdsf[i] = krg_[regionIdx].derivative(solventFraction_i);
    }

    ADB::M dkrgdsf_diag = spdiag(dkrgdsf);
    const int num_blocks = solventFraction.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        fastSparseProduct(dkrgdsf_diag, solventFraction.derivative()[block], jacs[block]);
    }
    return ADB::function(std::move(krg), std::move(jacs));
}

ADB SolventPropsAdFromDeck::solventRelPermMultiplier(const ADB& solventFraction,
                                 const Cells& cells) const
{
    const int n = cells.size();
    assert(solventFraction.value().size() == n);
    V krs(n);
    V dkrsdsf(n);
    for (int i = 0; i < n; ++i) {
        const double& solventFraction_i = solventFraction.value()[i];
        int regionIdx = 0; // TODO add mapping from cells to sat function table
        krs[i] = krs_[regionIdx](solventFraction_i);
        dkrsdsf[i] = krs_[regionIdx].derivative(solventFraction_i);
    }

    ADB::M dkrsdsf_diag = spdiag(dkrsdsf);
    const int num_blocks = solventFraction.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        fastSparseProduct(dkrsdsf_diag, solventFraction.derivative()[block], jacs[block]);
    }
    return ADB::function(std::move(krs), std::move(jacs));
}

V SolventPropsAdFromDeck::solventSurfaceDensity(const Cells& cells) const {
    const int n = cells.size();
    V density(n);
    for (int i = 0; i < n; ++i) {
        int regionIdx = cellPvtRegionIdx_[i];
        density[i] = solvent_surface_densities_[regionIdx];
    }
    return density;
}

} //namespace OPM
