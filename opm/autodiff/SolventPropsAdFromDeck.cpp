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
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>


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
        extractTableIndex("SATNUM", eclState, number_of_cells, global_cell, cellSatNumRegionIdx_);

        // surface densities
        if (deck->hasKeyword("SDENSITY")) {
            const auto& densityKeyword = deck->getKeyword("SDENSITY");
            int numRegions = densityKeyword.size();
            solvent_surface_densities_.resize(numRegions);
            for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                solvent_surface_densities_[regionIdx]
                        = densityKeyword.getRecord(regionIdx).getItem("SOLVENT_DENSITY").getSIDouble(0);
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


            // retrieve the cell specific Misc table index from the deck
            // and using the grid...
            extractTableIndex("MISCNUM", eclState, number_of_cells, global_cell, cellMiscRegionIdx_);

            // misicible hydrocabon relative permeability wrt water
            const TableContainer& sof2Tables = tables->getSof2Tables();
            if (!sof2Tables.empty()) {

                int numRegions = sof2Tables.size();

                // resize the attributes of the object
                krn_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::Sof2Table& sof2Table = sof2Tables.getTable<Sof2Table>(regionIdx);

                    // Copy data
                    // Sn = So + Sg + Ss;
                    const auto& sn = sof2Table.getSoColumn();
                    const auto& krn = sof2Table.getKroColumn();

                    krn_[regionIdx] = NonuniformTableLinear<double>(sn, krn);
                }

            } else {
                OPM_THROW(std::runtime_error, "SOF2 must be specified in MISCIBLE (SOLVENT) runs\n");
            }

            const TableContainer& miscTables = tables->getMiscTables();
            if (!miscTables.empty()) {

                int numRegions = miscTables.size();

                // resize the attributes of the object
                misc_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::MiscTable& miscTable = miscTables.getTable<MiscTable>(regionIdx);

                    // Copy data
                    // solventFraction = Ss / (Ss + Sg);
                    const auto& solventFraction = miscTable.getSolventFractionColumn();
                    const auto& misc = miscTable.getMiscibilityColumn();

                    misc_[regionIdx] = NonuniformTableLinear<double>(solventFraction, misc);

                }
            } else {
                OPM_THROW(std::runtime_error, "MISC must be specified in MISCIBLE (SOLVENT) runs\n");
            }

            // miscible relative permeability multipleiers
            const TableContainer& msfnTables = tables->getMsfnTables();
            if (!msfnTables.empty()) {

                int numRegions = msfnTables.size();

                // resize the attributes of the object
                mkrsg_.resize(numRegions);
                mkro_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::MsfnTable& msfnTable = msfnTables.getTable<MsfnTable>(regionIdx);

                    // Copy data
                    // Ssg = Ss + Sg;
                    const auto& Ssg = msfnTable.getGasPhaseFractionColumn();
                    const auto& krsg = msfnTable.getGasSolventRelpermMultiplierColumn();
                    const auto& kro = msfnTable.getOilRelpermMultiplierColumn();

                    mkrsg_[regionIdx] = NonuniformTableLinear<double>(Ssg, krsg);
                    mkro_[regionIdx] = NonuniformTableLinear<double>(Ssg, kro);

                }
            }

            const TableContainer& sorwmisTables = tables->getSorwmisTables();
            if (!sorwmisTables.empty()) {

                int numRegions = sorwmisTables.size();

                // resize the attributes of the object
                sorwmis_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::SorwmisTable& sorwmisTable = sorwmisTables.getTable<SorwmisTable>(regionIdx);

                    // Copy data
                    const auto& sw = sorwmisTable.getWaterSaturationColumn();
                    const auto& sorwmis = sorwmisTable.getMiscibleResidualOilColumn();

                    sorwmis_[regionIdx] = NonuniformTableLinear<double>(sw, sorwmis);
                }
            }

            const TableContainer& sgcwmisTables = tables->getSgcwmisTables();
            if (!sgcwmisTables.empty()) {

                int numRegions = sgcwmisTables.size();

                // resize the attributes of the object
                sgcwmis_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const Opm::SgcwmisTable& sgcwmisTable = sgcwmisTables.getTable<SgcwmisTable>(regionIdx);

                    // Copy data
                    const auto& sw = sgcwmisTable.getWaterSaturationColumn();
                    const auto& sgcwmis = sgcwmisTable.getMiscibleResidualGasColumn();

                    sgcwmis_[regionIdx] = NonuniformTableLinear<double>(sw, sgcwmis);
                }
            }

            if (deck->hasKeyword("TLMIXPAR")) {
                const int numRegions = deck->getKeyword("TLMIXPAR").size();

                // resize the attributes of the object
                mix_param_viscosity_.resize(numRegions);
                mix_param_density_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const auto& tlmixparRecord = deck->getKeyword("TLMIXPAR").getRecord(regionIdx);
                    const auto& mix_params_viscosity = tlmixparRecord.getItem("TL_VISCOSITY_PARAMETER").getSIDoubleData();
                    mix_param_viscosity_[regionIdx] = mix_params_viscosity[0];
                    const auto& mix_params_density = tlmixparRecord.getItem("TL_DENSITY_PARAMETER").getSIDoubleData();
                    const int numDensityItems = mix_params_density.size();
                    if (numDensityItems == 0) {
                        mix_param_density_[regionIdx] = mix_param_viscosity_[regionIdx];
                    } else if (numDensityItems == 1) {
                        mix_param_density_[regionIdx] = mix_params_density[0];
                    } else {
                        OPM_THROW(std::runtime_error, "Only one value can be entered for the TL parameter pr MISC region.");
                    }
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
        int regionIdx = cellPvtRegionIdx_[cells[i]];
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
    return SolventPropsAdFromDeck::makeADBfromTables(pg, cells, cellPvtRegionIdx_, b_);

}

ADB SolventPropsAdFromDeck::gasRelPermMultiplier(const ADB& solventFraction,
                                                 const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeADBfromTables(solventFraction, cells, cellSatNumRegionIdx_, krg_);

}

ADB SolventPropsAdFromDeck::solventRelPermMultiplier(const ADB& solventFraction,
                                                     const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeADBfromTables(solventFraction, cells, cellSatNumRegionIdx_, krs_);
}


ADB SolventPropsAdFromDeck::misicibleHydrocarbonWaterRelPerm(const ADB& Sn,
                                                             const Cells& cells) const
{
    return SolventPropsAdFromDeck::makeADBfromTables(Sn, cells, cellSatNumRegionIdx_, krn_);
}

ADB SolventPropsAdFromDeck::miscibleSolventGasRelPermMultiplier(const ADB& Ssg,
                                                                const Cells& cells) const
{
    if (mkrsg_.size() > 0) {
        return SolventPropsAdFromDeck::makeADBfromTables(Ssg, cells, cellSatNumRegionIdx_, mkrsg_);
    }
    // trivial function if not specified
    return Ssg;
}

ADB SolventPropsAdFromDeck::miscibleOilRelPermMultiplier(const ADB& So,
                                                         const Cells& cells) const
{
    if (mkro_.size() > 0) {
        return SolventPropsAdFromDeck::makeADBfromTables(So, cells, cellSatNumRegionIdx_, mkro_);
    }
    // trivial function if not specified
    return So;
}

ADB SolventPropsAdFromDeck::miscibilityFunction(const ADB& solventFraction,
                                                const Cells& cells) const
{

    return SolventPropsAdFromDeck::makeADBfromTables(solventFraction, cells, cellMiscRegionIdx_, misc_);
}


ADB SolventPropsAdFromDeck::miscibleCriticalGasSaturationFunction (const ADB& Sw,
                                                                   const Cells& cells) const {
    if (sgcwmis_.size()>0) {
        return SolventPropsAdFromDeck::makeADBfromTables(Sw, cells, cellMiscRegionIdx_, sgcwmis_);
    }
    // return zeros if not specified
    return ADB::constant(V::Zero(Sw.size()));
}


ADB SolventPropsAdFromDeck::miscibleResidualOilSaturationFunction (const ADB& Sw,
                                                                   const Cells& cells) const {
    if (sorwmis_.size()>0) {
        return SolventPropsAdFromDeck::makeADBfromTables(Sw, cells, cellMiscRegionIdx_, sorwmis_);
    }
    // return zeros if not specified
    return ADB::constant(V::Zero(Sw.size()));
}

ADB SolventPropsAdFromDeck::makeADBfromTables(const ADB& X_AD,
                                              const Cells& cells,
                                              const std::vector<int>& regionIdx,
                                              const std::vector<NonuniformTableLinear<double>>& tables) const {
    const int n = cells.size();
    assert(X_AD.value().size() == n);
    V x(n);
    V dx(n);
    for (int i = 0; i < n; ++i) {
        const double& X_i = X_AD.value()[i];
        x[i] = tables[regionIdx[cells[i]]](X_i);
        dx[i] = tables[regionIdx[cells[i]]].derivative(X_i);
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
        int regionIdx = cellPvtRegionIdx_[cells[i]];
        density[i] = solvent_surface_densities_[regionIdx];
    }
    return density;
}

V SolventPropsAdFromDeck::mixingParameterViscosity(const Cells& cells) const {
    const int n = cells.size();
    if (mix_param_viscosity_.size() > 0) {
        V mix_param(n);
        for (int i = 0; i < n; ++i) {
            int regionIdx = cellMiscRegionIdx_[cells[i]];
            mix_param[i] = mix_param_viscosity_[regionIdx];
        }
        return mix_param;
    }
    // return zeros if not specified
    return V::Zero(n);
}

V SolventPropsAdFromDeck::mixingParameterDensity(const Cells& cells) const {
    const int n = cells.size();
    if (mix_param_viscosity_.size() > 0) {
        V mix_param(n);
        for (int i = 0; i < n; ++i) {
            int regionIdx = cellMiscRegionIdx_[cells[i]];
            mix_param[i] = mix_param_density_[regionIdx];
        }
        return mix_param;
    }
    // return zeros if not specified
    return V::Zero(n);
}

void SolventPropsAdFromDeck::extractTableIndex(const std::string& keyword,
                                               Opm::EclipseStateConstPtr eclState,
                                               size_t numCompressed,
                                               const int* compressedToCartesianCellIdx,
                                               std::vector<int>& tableIdx) const {
    //Get the Region data
    const std::vector<int>& regionData = eclState->getIntGridProperty(keyword)->getData();
    // Convert this into an array of compressed cells
    // Eclipse uses Fortran-style indices which start at 1
    // instead of 0, we subtract 1.
    tableIdx.resize(numCompressed);
    for (size_t cellIdx = 0; cellIdx < numCompressed; ++ cellIdx) {
        size_t cartesianCellIdx = compressedToCartesianCellIdx ? compressedToCartesianCellIdx[cellIdx]:cellIdx;
        assert(cartesianCellIdx < regionData.size());
        tableIdx[cellIdx] = regionData[cartesianCellIdx] - 1;
    }
}


} //namespace OPM
