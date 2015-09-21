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
#ifndef OPM_THERMAL_GAS_PVT_WRAPPER_HPP
#define OPM_THERMAL_GAS_PVT_WRAPPER_HPP

#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <vector>

namespace Opm
{
    /// Class which wraps another (i.e., isothermal) PVT object into one which adds
    /// temperature dependence of gas
    class ThermalGasPvtWrapper : public PvtInterface
    {
    public:
        ThermalGasPvtWrapper()
        {}


        /// extract the quantities needed specify the temperature dependence of the gas
        /// viscosity and density from the deck
        void initFromDeck(std::shared_ptr<const PvtInterface> isothermalPvt,
                          Opm::DeckConstPtr deck,
                          Opm::EclipseStateConstPtr eclipseState)
        {
            isothermalPvt_ = isothermalPvt;
            auto tables = eclipseState->getTableManager();
            int numRegions;
            if (deck->hasKeyword("PVTG"))
                numRegions = tables->getPvtgTables().size();
            else if (deck->hasKeyword("PVDG"))
                numRegions = tables->getPvdgTables().size();
            else
                OPM_THROW(std::runtime_error, "Gas phase was not initialized using a known way");

            // viscosity
            if (deck->hasKeyword("GASVISCT")) {
                gasvisctTables_ = &tables->getGasvisctTables();
                assert(int(gasvisctTables_->size()) == numRegions);
                static_cast<void>(numRegions); //Silence compiler warning

                gasCompIdx_ = deck->getKeyword("GCOMPIDX")->getRecord(0)->getItem("GAS_COMPONENT_INDEX")->getInt(0) - 1;
                gasvisctColumnName_ = "Viscosity"+std::to_string(static_cast<long long>(gasCompIdx_));
            }

            // density
            if (deck->hasKeyword("TREF")) {
                tref_ = deck->getKeyword("TREF")->getRecord(0)->getItem("TEMPERATURE")->getSIDouble(0);
            }
        }

        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* T,
                        const double* z,
                        double* output_mu) const
        {
            if (gasvisctTables_)
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
            if (gasvisctTables_ != 0) {
                for (int i = 0; i < n; ++i) {
                    // temperature dependence of the gas phase. this assumes that the gas
                    // component index has been set properly, and it also looses the
                    // pressure dependence of gas. (This does not make much sense, but it
                    // seems to be what the documentation for the GASVISCT keyword in the
                    // RM says.)


                    int regionIdx = getPvtRegionIndex_(pvtRegionIdx, i);
                    double muGasvisct;
                    {
                        const GasvisctTable& gasvisctTable = gasvisctTables_->getTable<GasvisctTable>(regionIdx);
                        muGasvisct = gasvisctTable.evaluate(gasvisctColumnName_, T[i]);
                    }

                    output_mu[i] = muGasvisct;
                    output_dmudp[i] = 0.0;
                    output_dmudr[i] = 0.0;

                    // TODO (?): derivative of gas viscosity w.r.t. temperature.
                }
            }
            else {
                // compute the isothermal viscosity and its derivatives
                isothermalPvt_->mu(n, pvtRegionIdx, p, T, r, output_mu, output_dmudp, output_dmudr);
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
            if (gasvisctTables_ != 0) {
                for (int i = 0; i < n; ++i) {
                    // temperature dependence of the gas phase. this assumes that the gas
                    // component index has been set properly, and it also looses the
                    // pressure dependence of gas. (This does not make much sense, but it
                    // seems to be what the documentation for the GASVISCT keyword in the
                    // RM says.)
                    int regionIdx = getPvtRegionIndex_(pvtRegionIdx, i);
                    double muGasvisct;
                    {
                        const GasvisctTable& gasvisctTable = gasvisctTables_->getTable<GasvisctTable>(regionIdx);
                        muGasvisct = gasvisctTable.evaluate(gasvisctColumnName_, T[i]);
                    }

                    output_mu[i] = muGasvisct;
                    output_dmudp[i] = 0.0;
                    output_dmudr[i] = 0.0;

                    // TODO (?): derivative of gas viscosity w.r.t. temperature.
                }
            }
            else {
                // compute the isothermal viscosity and its derivatives
                isothermalPvt_->mu(n, pvtRegionIdx, p, T, r, cond, output_mu, output_dmudp, output_dmudr);
            }
        }

        virtual void B(const int n,
                       const int* pvtRegionIdx,
                       const double* p,
                       const double* T,
                       const double* z,
                       double* output_B) const
        {
            // isothermal case
            isothermalPvt_->B(n, pvtRegionIdx, p, T, z, output_B);

            if (tref_ > 0.0) {
                // the Eclipse TD/RM do not explicitly specify the relation of the gas
                // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
                // implies that the temperature dependence of the gas phase is rho(T, p) =
                // rho(tref_, p)/tref_*T ...
                for (int i = 0; i < n; ++i) {
                    double alpha = tref_/T[i];
                    output_B[i] *= alpha;
                }
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
            isothermalPvt_->dBdp(n, pvtRegionIdx, p, T, z, output_B, output_dBdp);

            if (tref_ > 0.0) {
                // the Eclipse TD/RM do not explicitly specify the relation of the gas
                // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
                // implies that the temperature dependence of the gas phase is rho(T, p) =
                // rho(tref_, p)/tref_*T ...
                for (int i = 0; i < n; ++i) {
                    double alpha = tref_/T[i];
                    output_B[i] *= alpha;
                    output_dBdp[i] *= alpha;
                }
            }
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
            isothermalPvt_->b(n, pvtRegionIdx, p, T, r, output_b, output_dbdp, output_dbdr);

            if (tref_ > 0.0) {
                // the Eclipse TD/RM do not explicitly specify the relation of the gas
                // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
                // implies that the temperature dependence of the gas phase is rho(T, p) =
                // rho(tref_, p)/tref_*T ...
                for (int i = 0; i < n; ++i) {
                    double alpha = T[i]/tref_;
                    output_b[i] *= alpha;
                    output_dbdp[i] *= alpha;
                    output_dbdr[i] *= alpha;
                }
            }
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
            isothermalPvt_->b(n, pvtRegionIdx, p, T, r, cond, output_b, output_dbdp, output_dbdr);

            if (tref_ > 0.0) {
                // the Eclipse TD/RM do not explicitly specify the relation of the gas
                // density and the temperature, but equation (69.49) (for Eclipse 2011.1)
                // implies that the temperature dependence of the gas phase is rho(T, p) =
                // rho(tref_, p)/tref_*T ...
                for (int i = 0; i < n; ++i) {
                    double alpha = T[i]/tref_;
                    output_b[i] *= alpha;
                    output_dbdp[i] *= alpha;
                    output_dbdr[i] *= alpha;
                }
            }
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
        int getPvtRegionIndex_(const int* pvtRegionIdx, int cellIdx) const
        {
            if (!pvtRegionIdx)
                return 0;
            return pvtRegionIdx[cellIdx];
        }

        // the PVT propertied for the isothermal case
        std::shared_ptr<const PvtInterface> isothermalPvt_;

        // The PVT properties needed for temperature dependence of the viscosity. We need
        // to store one value per PVT region.
        const TableContainer* gasvisctTables_;
        std::string gasvisctColumnName_;
        int gasCompIdx_;

        // The PVT properties needed for temperature dependence of the density.
        double tref_;
    };

}

#endif

