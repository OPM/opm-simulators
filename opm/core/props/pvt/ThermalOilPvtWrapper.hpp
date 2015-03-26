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

#ifndef OPM_THERMAL_OIL_PVT_WRAPPER_HPP
#define OPM_THERMAL_OIL_PVT_WRAPPER_HPP

#include <opm/core/props/pvt/PvtInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <vector>

namespace Opm
{
    /// Class which wraps another (i.e., isothermal) PVT object into one which adds
    /// temperature dependence of oil
    class ThermalOilPvtWrapper : public PvtInterface
    {
    public:
        ThermalOilPvtWrapper()
        {}


        /// set the tables which specify the temperature dependence of the oil viscosity
        void initFromDeck(std::shared_ptr<const PvtInterface> isothermalPvt,
                          Opm::DeckConstPtr deck,
                          Opm::EclipseStateConstPtr eclipseState)
        {
            isothermalPvt_ = isothermalPvt;

            int numRegions;
            if (deck->hasKeyword("PVTO"))
                numRegions = eclipseState->getPvtoTables().size();
            else if (deck->hasKeyword("PVDO"))
                numRegions = eclipseState->getPvdoTables().size();
            else if (deck->hasKeyword("PVCDO"))
                numRegions = deck->getKeyword("PVCDO")->size();
            else
                OPM_THROW(std::runtime_error, "Oil phase was not initialized using a known way");

            // viscosity
            if (deck->hasKeyword("VISCREF")) {
                oilvisctTables_ = &eclipseState->getOilvisctTables();
                Opm::DeckKeywordConstPtr viscrefKeyword = deck->getKeyword("VISCREF");

                assert(oilvisctTables_->size() == numRegions);
                assert(viscrefKeyword->size() == numRegions);

                viscrefPress_.resize(numRegions);
                viscrefRs_.resize(numRegions);
                muRef_.resize(numRegions);
                for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                    DeckRecordConstPtr viscrefRecord = viscrefKeyword->getRecord(regionIdx);
                    viscrefPress_[regionIdx] = viscrefRecord->getItem("REFERENCE_PRESSURE")->getSIDouble(0);
                    viscrefRs_[regionIdx] = viscrefRecord->getItem("REFERENCE_RS")->getSIDouble(0);

                    // temperature used to calculate the reference viscosity [K]. the
                    // value does not really matter if the underlying PVT object really
                    // is isothermal...
                    double Tref = 273.15 + 20;

                    // compute the reference viscosity using the isothermal PVT object.
                    double tmp1, tmp2;
                    isothermalPvt_->mu(1,
                                       &regionIdx,
                                       &viscrefPress_[regionIdx],
                                       &Tref,
                                       &viscrefRs_[regionIdx],
                                       &muRef_[regionIdx],
                                       &tmp1,
                                       &tmp2);
                }
            }

            // quantities required for density. note that we just always use the values
            // for the first EOS. (since EOS != PVT region.)
            tref_ = 0.0;
            if (deck->hasKeyword("THERMEX1")) {
                oilCompIdx_ = deck->getKeyword("OCOMPIDX")->getRecord(0)->getItem("OIL_COMPONENT_INDEX")->getInt(0) - 1;

                // always use the values of the first EOS
                tref_ = deck->getKeyword("TREF")->getRecord(0)->getItem("TEMPERATURE")->getSIDouble(oilCompIdx_);
                pref_ = deck->getKeyword("PREF")->getRecord(0)->getItem("PRESSURE")->getSIDouble(oilCompIdx_);
                cref_ = deck->getKeyword("CREF")->getRecord(0)->getItem("COMPRESSIBILITY")->getSIDouble(oilCompIdx_);
                thermex1_ = deck->getKeyword("THERMEX1")->getRecord(0)->getItem("EXPANSION_COEFF")->getSIDouble(oilCompIdx_);
            }
        }

        virtual void mu(const int n,
                        const int* pvtRegionIdx,
                        const double* p,
                        const double* T,
                        const double* z,
                        double* output_mu) const
        {
            if (oilvisctTables_)
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

            if (!oilvisctTables_)
                // isothermal case
                return;

            // temperature dependence
            for (int i = 0; i < n; ++i) {
                int regionIdx = getPvtRegionIndex_(pvtRegionIdx, i);

                // calculate the viscosity of the isothermal keyword for the reference
                // pressure given by the VISCREF keyword.
                double muRef = muRef_[regionIdx];

                // compute the viscosity deviation due to temperature
                double muOilvisct = (*oilvisctTables_)[regionIdx].evaluate("Viscosity", T[i]);
                double alpha = muOilvisct/muRef;

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

            if (!oilvisctTables_)
                // isothermal case
                return;

            // temperature dependence
            for (int i = 0; i < n; ++i) {
                int regionIdx = getPvtRegionIndex_(pvtRegionIdx, i);

                // calculate the viscosity of the isothermal keyword for the reference
                // pressure given by the VISCREF keyword.
                double muRef = muRef_[regionIdx];

                // compute the viscosity deviation due to temperature
                double muOilvisct = (*oilvisctTables_)[regionIdx].evaluate("Viscosity", T[i]);
                double alpha = muOilvisct/muRef;

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
            // isothermal case
            isothermalPvt_->B(n, pvtRegionIdx, p, T, z, output_B);

            if (thermex1_ <= 0.0)
                // isothermal case
                return;

            // deal with the temperature dependence of the oil phase. we use equation
            // (3.208) from the Eclipse 2011.1 Reference Manual, but we calculate rho_ref
            // using the isothermal keyword instead of using the value for the
            // components, so the oil compressibility is already dealt with there. Note
            // that we only do the part for the oil component here, the part for
            // dissolved gas is ignored so far.
            double cT1 = thermex1_;
            double TRef = tref_;
            for (int i = 0; i < n; ++i) {
                double alpha = (1 + cT1*(T[i] - TRef));
                output_B[i] *= alpha;
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

            if (thermex1_ <= 0.0)
                // isothermal case
                return;

            // deal with the temperature dependence of the oil phase. we use equation
            // (3.208) from the Eclipse 2011.1 Reference Manual, but we calculate rho_ref
            // using the isothermal keyword instead of using the value for the
            // components, so the oil compressibility is already dealt with there. Note
            // that we only do the part for the oil component here, the part for
            // dissolved gas is ignored so far.
            double cT1 = thermex1_;
            double TRef = tref_;
            for (int i = 0; i < n; ++i) {
                double alpha = (1 + cT1*(T[i] - TRef));
                output_B[i] *= alpha;
                output_dBdp[i] *= alpha;
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

            if (thermex1_ <= 0.0)
                // isothermal case
                return;

            // deal with the temperature dependence of the oil phase. we use equation
            // (3.208) from the Eclipse 2011.1 Reference Manual, but we calculate rho_ref
            // using the isothermal keyword instead of using the value for the
            // components, so the oil compressibility is already dealt with there. Note
            // that we only do the part for the oil component here, the part for
            // dissolved gas is ignored so far.
            double cT1 = thermex1_;
            double TRef = tref_;
            for (int i = 0; i < n; ++i) {
                double alpha = 1.0/(1 + cT1*(T[i] - TRef));
                output_b[i] *= alpha;
                output_dbdp[i] *= alpha;
                output_dbdr[i] *= alpha;
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

            if (thermex1_ <= 0.0)
                // isothermal case
                return;

            // deal with the temperature dependence of the oil phase. we use equation
            // (3.208) from the Eclipse 2011.1 Reference Manual, but we calculate rho_ref
            // using the isothermal keyword instead of using the value for the
            // components, so the oil compressibility is already dealt with there. Note
            // that we only do the part for the oil component here, the part for
            // dissolved gas is ignored so far.
            double cT1 = thermex1_;
            double TRef = tref_;
            for (int i = 0; i < n; ++i) {
                double alpha = 1.0/(1 + cT1*(T[i] - TRef));
                output_b[i] *= alpha;
                output_dbdp[i] *= alpha;
                output_dbdr[i] *= alpha;
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
        std::vector<double> viscrefPress_;
        std::vector<double> viscrefRs_;
        std::vector<double> muRef_;

        const std::vector<Opm::OilvisctTable>* oilvisctTables_;

        // The PVT properties needed for temperature dependence of the density. This is
        // specified as one value per EOS in the manual, but we unconditionally use the
        // expansion coefficient of the first EOS...
        int oilCompIdx_;
        double tref_;
        double pref_;
        double cref_;
        double thermex1_;
    };

}

#endif

