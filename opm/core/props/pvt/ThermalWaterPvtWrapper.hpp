/*
  Copyright 2015 Andreas Lauser

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

#ifndef OPM_THERMAL_WATER_PVT_WRAPPER_HPP
#define OPM_THERMAL_WATER_PVT_WRAPPER_HPP

#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/WatvisctTable.hpp>

#include <vector>

namespace Opm
{
    /// Class which wraps another (i.e., isothermal) PVT object into one which adds
    /// temperature dependence of water
    class ThermalWaterPvtWrapper : public PvtInterface
    {
    public:
        ThermalWaterPvtWrapper()
        {}


        /// set the tables which specify the temperature dependence of the water viscosity
        void initFromDeck(std::shared_ptr<const PvtInterface> isothermalPvt,
                          Opm::DeckConstPtr deck,
                          Opm::EclipseStateConstPtr eclipseState)
        {
            isothermalPvt_ = isothermalPvt;
            watvisctTables_ = 0;

            // stuff which we need to get from the PVTW keyword
            const auto& pvtwKeyword = deck->getKeyword("PVTW");
            int numRegions = pvtwKeyword.size();
            pvtwRefPress_.resize(numRegions);
            pvtwRefB_.resize(numRegions);
            pvtwCompressibility_.resize(numRegions);
            pvtwViscosity_.resize(numRegions);
            pvtwViscosibility_.resize(numRegions);
            for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                const auto& pvtwRecord = pvtwKeyword.getRecord(regionIdx);
                pvtwRefPress_[regionIdx] = pvtwRecord.getItem("P_REF").getSIDouble(0);
                pvtwRefB_[regionIdx] = pvtwRecord.getItem("WATER_VOL_FACTOR").getSIDouble(0);
                pvtwViscosity_[regionIdx] = pvtwRecord.getItem("WATER_VISCOSITY").getSIDouble(0);
                pvtwViscosibility_[regionIdx] = pvtwRecord.getItem("WATER_VISCOSIBILITY").getSIDouble(0);
            }

            // quantities required for the temperature dependence of the viscosity
            // (basically we expect well-behaved VISCREF and WATVISCT keywords.)
            if (deck->hasKeyword("VISCREF")) {
                auto tables = eclipseState->getTableManager();
                watvisctTables_ = &tables->getWatvisctTables();
                const auto& viscrefKeyword = deck->getKeyword("VISCREF");

                assert(int(watvisctTables_->size()) == numRegions);
                assert(int(viscrefKeyword.size()) == numRegions);

                viscrefPress_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                    const auto& viscrefRecord = viscrefKeyword.getRecord(regionIdx);

                    viscrefPress_[regionIdx] = viscrefRecord.getItem("REFERENCE_PRESSURE").getSIDouble(0);
                }
            }

            // quantities required for the temperature dependence of the density
            if (deck->hasKeyword("WATDENT")) {
                const auto& watdentKeyword = deck->getKeyword("WATDENT");

                assert(int(watdentKeyword.size()) == numRegions);

                watdentRefTemp_.resize(numRegions);
                watdentCT1_.resize(numRegions);
                watdentCT2_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    const auto& watdentRecord = watdentKeyword.getRecord(regionIdx);

                    watdentRefTemp_[regionIdx] = watdentRecord.getItem("REFERENCE_TEMPERATURE").getSIDouble(0);
                    watdentCT1_[regionIdx] = watdentRecord.getItem("EXPANSION_COEFF_LINEAR").getSIDouble(0);
                    watdentCT2_[regionIdx] = watdentRecord.getItem("EXPANSION_COEFF_QUADRATIC").getSIDouble(0);
                }
            }
        }

        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* T,
                        const double* z,
                        double* output_mu) const
        {
            if (watvisctTables_)
                // TODO: temperature dependence for viscosity depending on z
                OPM_THROW(std::runtime_error,
                          "temperature dependent viscosity as a function of z "
                          "is not yet implemented!");

            // compute the isothermal viscosity
            isothermalPvt_->mu(n, pvtRegionIdx, p, T, z, output_mu);
        }

        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* T,
                        const double* r,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const
        {
            // compute the isothermal viscosity and its derivatives
            isothermalPvt_->mu(n, pvtRegionIdx, p, T, r, output_mu, output_dmudp, output_dmudr);

            if (!watvisctTables_)
                // isothermal case
                return;

            // temperature dependence
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);

                // calculate the viscosity of the isothermal keyword for the reference
                // pressure given by the VISCREF keyword.
                double x = -pvtwViscosibility_[tableIdx]*(viscrefPress_[tableIdx] - pvtwRefPress_[tableIdx]);
                double muRef = pvtwViscosity_[tableIdx]/(1.0 + x + 0.5*x*x);

                // compute the viscosity deviation due to temperature
                double alpha;
                {
                    const WatvisctTable& watVisctTable = watvisctTables_->getTable<WatvisctTable>(tableIdx);
                    double muWatvisct = watVisctTable.evaluate("Viscosity", T[i]);
                    alpha = muWatvisct/muRef;
                }

                output_mu[i] *= alpha;
                output_dmudp[i] *= alpha;
                output_dmudr[i] *= alpha;
                // TODO (?): derivative of viscosity w.r.t. temperature.
            }
        }

        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* T,
                        const double* r,
                        const PhasePresence* cond,
                        double* output_mu,
                        double* output_dmudp,
                        double* output_dmudr) const
        {
            // compute the isothermal viscosity and its derivatives
            isothermalPvt_->mu(n, pvtRegionIdx, p, T, r, cond, output_mu, output_dmudp, output_dmudr);

            if (!watvisctTables_)
                // isothermal case
                return;

            // temperature dependence
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);

                // calculate the viscosity of the isothermal keyword for the reference
                // pressure given by the VISCREF keyword.
                double x = -pvtwViscosibility_[tableIdx]*(viscrefPress_[tableIdx] - pvtwRefPress_[tableIdx]);
                double muRef = pvtwViscosity_[tableIdx]/(1.0 + x + 0.5*x*x);

                // compute the viscosity deviation due to temperature
                double alpha;
                {
                    const WatvisctTable& watVisctTable = watvisctTables_->getTable<WatvisctTable>(tableIdx);
                    double muWatvisct = watVisctTable.evaluate("Viscosity", T[i]);
                    alpha = muWatvisct/muRef;
                }
                output_mu[i] *= alpha;
                output_dmudp[i] *= alpha;
                output_dmudr[i] *= alpha;
                // TODO (?): derivative of viscosity w.r.t. temperature.
            }
        }

        virtual void B(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* T,
                       const double* z,
                       double* output_B) const
        {
            if (watdentRefTemp_.empty()) {
                // isothermal case
                isothermalPvt_->B(n, pvtRegionIdx, p, T, z, output_B);
                return;
            }

            // This changes how the water density depends on pressure compared to what's
            // used for the PVTW keyword, but it seems to be what Eclipse does. For
            // details, see the documentation for the WATDENT keyword in the Eclipse RM.
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);
                double BwRef = pvtwRefB_[tableIdx];
                double TRef = watdentRefTemp_[tableIdx];
                double X = pvtwCompressibility_[tableIdx]*(p[i] - pvtwRefPress_[tableIdx]);
                double cT1 = watdentCT1_[tableIdx];
                double cT2 = watdentCT2_[tableIdx];
                double Y = T[i] - TRef;
                double Bw = BwRef*(1 - X)*(1 + cT1*Y + cT2*Y*Y);
                output_B[i] = Bw;
            }
        }

        virtual void dBdp(const int n,
                          const int* pvtRegionIdx,
                          const double* p,
                          const double* T,
                          const double* z,
                          double* output_B,
                          double* output_dBdp) const
        {
            if (watdentRefTemp_.empty()) {
                // isothermal case
                isothermalPvt_->dBdp(n, pvtRegionIdx, p, T, z, output_B, output_dBdp);
                return;
            }

            // This changes how the water density depends on pressure. This is awkward,
            // but it seems to be what Eclipse does. See the documentation for the
            // WATDENT keyword in the Eclipse RM
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);
                double BwRef = pvtwRefB_[tableIdx];
                double TRef = watdentRefTemp_[tableIdx];
                double X = pvtwCompressibility_[tableIdx]*(p[i] - pvtwRefPress_[tableIdx]);
                double cT1 = watdentCT1_[tableIdx];
                double cT2 = watdentCT2_[tableIdx];
                double Y = T[i] - TRef;
                double Bw = BwRef*(1 - X)*(1 + cT1*Y + cT2*Y*Y);
                output_B[i] = Bw;
            }

            std::fill(output_dBdp, output_dBdp + n, 0.0);
        }

        virtual void b(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* T,
                       const double* r,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const
        {
            if (watdentRefTemp_.empty()) {
                // isothermal case
                isothermalPvt_->b(n, pvtRegionIdx, p, T, r, output_b, output_dbdp, output_dbdr);
                return;
            }

            // This changes how the water density depends on pressure. This is awkward,
            // but it seems to be what Eclipse does. See the documentation for the
            // WATDENT keyword in the Eclipse RM
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);
                double BwRef = pvtwRefB_[tableIdx];
                double TRef = watdentRefTemp_[tableIdx];
                double X = pvtwCompressibility_[tableIdx]*(p[i] - pvtwRefPress_[tableIdx]);
                double cT1 = watdentCT1_[tableIdx];
                double cT2 = watdentCT2_[tableIdx];
                double Y = T[i] - TRef;
                double Bw = BwRef*(1 - X)*(1 + cT1*Y + cT2*Y*Y);
                output_b[i] = 1.0/Bw;
            }

            std::fill(output_dbdp, output_dbdp + n, 0.0);
            std::fill(output_dbdr, output_dbdr + n, 0.0);
        }

        virtual void b(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* T,
                       const double* r,
                       const PhasePresence* cond,
                       double* output_b,
                       double* output_dbdp,
                       double* output_dbdr) const
        {
            if (watdentRefTemp_.empty()) {
                // isothermal case
                isothermalPvt_->b(n, pvtRegionIdx, p, T, r, cond, output_b, output_dbdp, output_dbdr);
                return;
            }

            // This changes pressure dependence of the water density, but it seems to be
            // what Eclipse does. See the documentation for the WATDENT keyword in the
            // Eclipse RM
            for (int i = 0; i < n; ++i) {
                int tableIdx = getTableIndex_(pvtRegionIdx, i);
                double BwRef = pvtwRefB_[tableIdx];
                double TRef = watdentRefTemp_[tableIdx];
                double X = pvtwCompressibility_[tableIdx]*(p[i] - pvtwRefPress_[tableIdx]);
                double cT1 = watdentCT1_[tableIdx];
                double cT2 = watdentCT2_[tableIdx];
                double Y = T[i] - TRef;
                double Bw = BwRef*(1 - X)*(1 + cT1*Y + cT2*Y*Y);
                output_b[i] = 1.0/Bw;
            }

            std::fill(output_dbdp, output_dbdp + n, 0.0);
            std::fill(output_dbdr, output_dbdr + n, 0.0);

        }

        virtual void rsSat(const int n,
                           const int* pvtRegionIdx,
                           const double* p,
                           double* output_rsSat,
                           double* output_drsSatdp) const
        {
            isothermalPvt_->rsSat(n, pvtRegionIdx, p, output_rsSat, output_drsSatdp);
        }

        virtual void rvSat(const int n,
                           const int* pvtRegionIdx,
                           const double* p,
                           double* output_rvSat,
                           double* output_drvSatdp) const
        {
            isothermalPvt_->rvSat(n, pvtRegionIdx, p, output_rvSat, output_drvSatdp);
        }

        virtual void R(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* z,
                       double* output_R) const
        {
            isothermalPvt_->R(n, pvtRegionIdx, p, z, output_R);
        }

        virtual void dRdp(const int n,
                          const int* pvtRegionIdx,
                          const double* p,
                          const double* z,
                          double* output_R,
                          double* output_dRdp) const
        {
            isothermalPvt_->dRdp(n, pvtRegionIdx, p, z, output_R, output_dRdp);
        }

    private:
        int getTableIndex_(const int* pvtTableIdx, int cellIdx) const
        {
            if (!pvtTableIdx)
                return 0;
            return pvtTableIdx[cellIdx];
        }

        // the PVT propertied for the isothermal case
        std::shared_ptr<const PvtInterface> isothermalPvt_;

        // The PVT properties needed for temperature dependence. We need to store one
        // value per PVT region.
        std::vector<double> viscrefPress_;

        std::vector<double> watdentRefTemp_;
        std::vector<double> watdentCT1_;
        std::vector<double> watdentCT2_;

        std::vector<double> pvtwRefPress_;
        std::vector<double> pvtwRefB_;
        std::vector<double> pvtwCompressibility_;
        std::vector<double> pvtwViscosity_;
        std::vector<double> pvtwViscosibility_;

        const TableContainer* watvisctTables_;
    };

}

#endif

