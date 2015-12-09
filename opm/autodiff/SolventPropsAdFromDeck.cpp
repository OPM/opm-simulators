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

#include <opm/autodiff/SolventPropsAdFromDeck.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <opm/parser/eclipse/EclipseState/Tables/PvdsTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SsfnTable.hpp>


namespace Opm
{
// Making these typedef to make the code more readable.
typedef SolventPropsAdFromDeck::ADB ADB;
typedef Eigen::SparseMatrix<double> S;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> D;
typedef SolventPropsAdFromDeck::V V;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

SolventPropsAdFromDeck::SolventPropsAdFromDeck(DeckConstPtr deck,
                                                     EclipseStateConstPtr eclState,
                                                     const int number_of_cells,
                                                     const int* global_cell)
{
    if (deck->hasKeyword("SOLVENT")) {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, eclState, number_of_cells, global_cell);

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

        auto tables = eclState->getTableManager();
        // pvt
        const TableContainer& pvdsTables = tables->getPvdsTables();
        if (!pvdsTables.empty()) {

            int numRegions = pvdsTables.size();
            // resize the attributes of the object
            b_.resize(numRegions);
            viscosity_.resize(numRegions);
            inverseBmu_.resize(numRegions);

            for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const Opm::PvdsTable& pvdsTable = pvdsTables.getTable<PvdsTable>(regionIdx);

                const auto& press = pvdsTable.getPressureColumn();
                const auto& b = pvdsTable.getFormationFactorColumn();
                const auto& visc = pvdsTable.getViscosityColumn();

                const int sz = b.size();
                std::vector<double> inverseBmu(sz);
                std::vector<double> inverseB(sz);
                for (int i = 0; i < sz; ++i) {
                    inverseB[i] = 1.0 / b[i];
                    inverseBmu[i] = 1.0 / (b[i] * visc[i]);
                }


                b_[regionIdx] = NonuniformTableLinear<double>(press, inverseB);
                viscosity_[regionIdx] = NonuniformTableLinear<double>(press, visc);
                inverseBmu_[regionIdx] = NonuniformTableLinear<double>(press, inverseBmu);
            }
        } else {
            OPM_THROW(std::runtime_error, "PVDS must be specified in SOLVENT runs\n");
        }

        const TableContainer& ssfnTables = tables->getSsfnTables();
        // relative permeabilty multiplier
        if (!ssfnTables.empty()) {

            int numRegions = ssfnTables.size();

            if(numRegions > 1) {
                OPM_THROW(std::runtime_error, "Only single table saturation function supported for SSFN");
            }
            // resize the attributes of the object
            krg_.resize(numRegions);
            krs_.resize(numRegions);
            for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const Opm::SsfnTable& ssfnTable = ssfnTables.getTable<SsfnTable>(regionIdx);

                // Copy data
                const auto& solventFraction = ssfnTable.getSolventFractionColumn();
                const auto& krg = ssfnTable.getGasRelPermMultiplierColumn();
                const auto& krs = ssfnTable.getSolventRelPermMultiplierColumn();

                krg_[regionIdx] = NonuniformTableLinear<double>(solventFraction, krg);
                krs_[regionIdx] = NonuniformTableLinear<double>(solventFraction, krs);
            }

        } else {
            OPM_THROW(std::runtime_error, "SSFN must be specified in SOLVENT runs\n");
        }


        if (deck->hasKeyword("MISCIBLE") ) {
            // misicible hydrocabon relative permeability wrt water
            const TableContainer& sof2Tables = tables->getSof2Tables();
            if (!sof2Tables.empty()) {

                int numRegions = sof2Tables.size();

                if(numRegions > 1) {
                    OPM_THROW(std::runtime_error, "Only single table saturation function supported for SOF2");
                }
                // resize the attributes of the object
                krn_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::Sof2Table& sof2Table = sof2Tables.getTable<Sof2Table>(regionIdx);

                    // Copy data
                    // Sn = So + Sg + Ss;
                    const std::vector<double>& sn = sof2Table.getSoColumn();
                    const std::vector<double>& krn = sof2Table.getKroColumn();

                    for (size_t i = 0; i < sn.size(); ++i) {
                        std::cout << sn[i] << " " << krn[i] <<std::endl;
                    }


                    krn_[regionIdx] = NonuniformTableLinear<double>(sn, krn);
                }

            } else {
                OPM_THROW(std::runtime_error, "SOF2 must be specified in MISCIBLE (SOLVENT) runs\n");
            }

            const TableContainer& miscTables = tables->getMiscTables();
            if (!miscTables.empty()) {

                int numRegions = miscTables.size();

                if(numRegions > 1) {
                    OPM_THROW(std::runtime_error, "Only single table miscibility function supported for MISC");
                }
                // resize the attributes of the object
                misc_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::MiscTable& miscTable = miscTables.getTable<MiscTable>(regionIdx);

                    // Copy data
                    // solventFraction = Ss / (Ss + Sg);
                    const std::vector<double>& solventFraction = miscTable.getSolventFractionColumn();
                    const std::vector<double>& misc = miscTable.getMiscibilityColumn();

                    misc_[regionIdx] = NonuniformTableLinear<double>(solventFraction, misc);

                }
            } else {
                OPM_THROW(std::runtime_error, "MISC must be specified in MISCIBLE (SOLVENT) runs\n");
            }

            // miscible relative permeability multipleiers
            const TableContainer& msfnTables = tables->getMsfnTables();
            if (!msfnTables.empty()) {

                int numRegions = msfnTables.size();

                if(numRegions > 1) {
                    OPM_THROW(std::runtime_error, "Only single table saturation function supported for MSFN");
                }
                // resize the attributes of the object
                mkrsg_.resize(numRegions);
                mkro_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::MsfnTable& msfnTable = msfnTables.getTable<MsfnTable>(regionIdx);

                    // Copy data
                    // Ssg = Ss + Sg;
                    const std::vector<double>& Ssg = msfnTable.getGasPhaseFractionColumn();
                    const std::vector<double>& krsg = msfnTable.getGasSolventRelpermMultiplierColumn();
                    const std::vector<double>& kro = msfnTable.getOilRelpermMultiplierColumn();

                    mkrsg_[regionIdx] = NonuniformTableLinear<double>(Ssg, krsg);
                    mkro_[regionIdx] = NonuniformTableLinear<double>(Ssg, kro);

                }
            }

            const TableContainer& sorwmisTables = tables->getSorwmisTables();
            if (!sorwmisTables.empty()) {

                int numRegions = sorwmisTables.size();

                if(numRegions > 1) {
                    OPM_THROW(std::runtime_error, "Only single table miscibility function supported for SORWMIS");
                }
                // resize the attributes of the object
                sorwmis_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::SorwmisTable& sorwmisTable = sorwmisTables.getTable<SorwmisTable>(regionIdx);

                    // Copy data
                    const std::vector<double>& sw = sorwmisTable.getWaterSaturationColumn();
                    const std::vector<double>& sorwmis = sorwmisTable.getMiscibleResidualOilColumn();

                    sorwmis_[regionIdx] = NonuniformTableLinear<double>(sw, sorwmis);
                }
            }

            const TableContainer& sgcwmisTables = tables->getSgcwmisTables();
            if (!sgcwmisTables.empty()) {

                int numRegions = sgcwmisTables.size();

                if(numRegions > 1) {
                    OPM_THROW(std::runtime_error, "Only single table miscibility function supported for SGCWMIS");
                }
                // resize the attributes of the object
                sgcwmis_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::SgcwmisTable& sgcwmisTable = sgcwmisTables.getTable<SgcwmisTable>(regionIdx);

                    // Copy data
                    const std::vector<double>& sw = sgcwmisTable.getWaterSaturationColumn();
                    const std::vector<double>& sgcwmis = sgcwmisTable.getMiscibleResidualGasColumn();

                    sgcwmis_[regionIdx] = NonuniformTableLinear<double>(sw, sgcwmis);
                }
            }

        }
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

    ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
    const int num_blocks = pg.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        jacs[block] = dmudp_diag * pg.derivative()[block];
    }
    return ADB::function(std::move(mu), std::move(jacs));
}

ADB SolventPropsAdFromDeck::bSolvent(const ADB& pg,
                                const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeAD(pg, cells, b_);

}

ADB SolventPropsAdFromDeck::gasRelPermMultiplier(const ADB& solventFraction,
                                 const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeAD(solventFraction, cells, krg_);

}

ADB SolventPropsAdFromDeck::solventRelPermMultiplier(const ADB& solventFraction,
                                 const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeAD(solventFraction, cells, krs_);
}


ADB SolventPropsAdFromDeck::misicibleHydrocarbonWaterRelPerm(const ADB& Sn,
                                 const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeAD(Sn, cells, krn_);
}

ADB SolventPropsAdFromDeck::miscibleSolventGasRelPermMultiplier(const ADB& Ssg,
                                 const Cells& cells) const
{
    if (mkrsg_.size() > 0) {
        return SolventPropsAdFromDeck::makeAD(Ssg, cells, mkrsg_);
    }
    // trivial function if not specified
    return Ssg;
}

ADB SolventPropsAdFromDeck::miscibleOilRelPermMultiplier(const ADB& So,
                                 const Cells& cells) const
{
    if (mkro_.size() > 0) {
        return SolventPropsAdFromDeck::makeAD(So, cells, mkro_);
    }
    // trivial function if not specified
    return So;
}

ADB SolventPropsAdFromDeck::miscibilityFunction(const ADB& solventFraction,
                                 const Cells& cells) const
{

    return SolventPropsAdFromDeck::makeAD(solventFraction, cells, misc_);
}


ADB SolventPropsAdFromDeck::miscibleCriticalGasSaturationFunction (const ADB& Sw,
                                           const Cells& cells) const {
    if (sgcwmis_.size()>0) {
        return SolventPropsAdFromDeck::makeAD(Sw, cells, sgcwmis_);
    }
    // return zeros if not specified
    return ADB::constant(V::Zero(Sw.size()));
}


ADB SolventPropsAdFromDeck::miscibleResidualOilSaturationFunction (const ADB& Sw,
                                           const Cells& cells) const {
    if (sorwmis_.size()>0) {
        return SolventPropsAdFromDeck::makeAD(Sw, cells, sorwmis_);
    }
    // return zeros if not specified
    return ADB::constant(V::Zero(Sw.size()));
}

ADB SolventPropsAdFromDeck::makeAD(const ADB& X_AD, const Cells& cells, std::vector<NonuniformTableLinear<double>> table) const {
    const int n = cells.size();
    assert(Sn.value().size() == n);
    V x(n);
    V dx(n);
    for (int i = 0; i < n; ++i) {
        const double& X_i = X_AD.value()[i];
        int regionIdx = 0; // TODO add mapping from cells to sat function table
        x[i] = table[regionIdx](X_i);
        dx[i] = table[regionIdx].derivative(X_i);
    }

    ADB::M dx_diag(dx.matrix().asDiagonal());
    const int num_blocks = X_AD.numBlocks();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        fastSparseProduct(dx_diag, X_AD.derivative()[block], jacs[block]);
    }
    return ADB::function(std::move(x), std::move(jacs));
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
