/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/Point2D.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

namespace Opm
{
    double PolymerProperties::cMax() const
    {
        return c_max_;
    }

    double PolymerProperties::mixParam() const
    {
        return mix_param_;
    }

    double PolymerProperties::rockDensity() const
    {
        return rock_density_;
    }

    double PolymerProperties::deadPoreVol() const
    {
        return dead_pore_vol_;
    }

    double PolymerProperties::resFactor() const
    {
        return res_factor_;
    }

    double PolymerProperties::cMaxAds() const
    {
        return c_max_ads_;
    }

    int PolymerProperties::adsIndex() const
    {
        return ads_index_;
    }

    const std::vector<double>&
    PolymerProperties::shearWaterVelocity() const
    {
        return water_vel_vals_;
    }

    const std::vector<double>&
    PolymerProperties::shearViscosityReductionFactor() const
    {
        return shear_vrf_vals_;
    }

    double PolymerProperties:: plyshlogRefConc() const
    {
        return plyshlog_ref_conc_;
    }

    bool PolymerProperties::hasPlyshlogRefSalinity() const
    {
        return has_plyshlog_ref_salinity_;
    }

    bool PolymerProperties::hasPlyshlogRefTemp() const
    {
        return has_plyshlog_ref_temp_;
    }

    double PolymerProperties::plyshlogRefSalinity() const
    {
        return plyshlog_ref_salinity_;
    }

    double PolymerProperties::plyshlogRefTemp() const
    {
        return plyshlog_ref_temp_;
    }

    bool PolymerProperties::hasPlyshlog() const
    {
        return has_plyshlog_;
    }

    bool PolymerProperties::hasShrate() const
    {
        return has_shrate_;
    }


    double PolymerProperties::shrate() const
    {
        return shrate_;
    }


    double
    PolymerProperties::shearVrf(const double velocity) const
    {
        return Opm::linearInterpolation(water_vel_vals_, shear_vrf_vals_, velocity);
    }

    double
    PolymerProperties::shearVrfWithDer(const double velocity, double& der) const
    {
        der =  Opm::linearInterpolationDerivative(water_vel_vals_, shear_vrf_vals_, velocity);
        return Opm::linearInterpolation(water_vel_vals_, shear_vrf_vals_, velocity);
    }

    double PolymerProperties::viscMult(double c) const
    {
        return Opm::linearInterpolation(c_vals_visc_, visc_mult_vals_, c);
    }

    double PolymerProperties::viscMultWithDer(double c, double* der) const
    {
        *der = Opm::linearInterpolationDerivative(c_vals_visc_, visc_mult_vals_, c);
        return Opm::linearInterpolation(c_vals_visc_, visc_mult_vals_, c);
    }

    void PolymerProperties::simpleAdsorption(double c, double& c_ads) const
    {
        double dummy;
        simpleAdsorptionBoth(c, c_ads, dummy, false);
    }

    void PolymerProperties::simpleAdsorptionWithDer(double c, double& c_ads,
                                                    double& dc_ads_dc) const
    {
        simpleAdsorptionBoth(c, c_ads, dc_ads_dc, true);
    }

    void PolymerProperties::simpleAdsorptionBoth(double c, double& c_ads,
                                                 double& dc_ads_dc, bool if_with_der) const
    {
        c_ads = Opm::linearInterpolation(c_vals_ads_, ads_vals_, c);;
        if (if_with_der) {
            dc_ads_dc = Opm::linearInterpolationDerivative(c_vals_ads_, ads_vals_, c);
        } else {
            dc_ads_dc = 0.;
        }
    }

    void PolymerProperties::adsorption(double c, double cmax, double& c_ads) const
    {
        double dummy;
        adsorptionBoth(c, cmax, c_ads, dummy, false);
    }

    void PolymerProperties::adsorptionWithDer(double c, double cmax,
                                              double& c_ads, double& dc_ads_dc) const
    {
        adsorptionBoth(c, cmax, c_ads, dc_ads_dc, true);
    }

    void PolymerProperties::adsorptionBoth(double c, double cmax,
                                           double& c_ads, double& dc_ads_dc,
                                           bool if_with_der) const
    {
        if (ads_index_ == Desorption) {
            simpleAdsorptionBoth(c, c_ads, dc_ads_dc, if_with_der);
        } else if (ads_index_ == NoDesorption) {
            simpleAdsorptionBoth(std::max(c, cmax), c_ads, dc_ads_dc, if_with_der);
        } else {
            OPM_THROW(std::runtime_error, "Invalid Adsoption index");
        }
    }


    void PolymerProperties::effectiveVisc(const double c, const double visc, double& mu_w_eff) const {
        effectiveInvVisc(c, visc, mu_w_eff);
        mu_w_eff = 1./mu_w_eff;
    }

    void PolymerProperties::effectiveViscWithDer(const double c, const double visc, double& mu_w_eff, double dmu_w_eff_dc) const {
        effectiveInvViscWithDer(c, visc, mu_w_eff, dmu_w_eff_dc);
        mu_w_eff = 1./mu_w_eff;
        dmu_w_eff_dc = -dmu_w_eff_dc*mu_w_eff*mu_w_eff;
    }

    void PolymerProperties::effectiveInvVisc(const double c, const double visc, double& inv_mu_w_eff) const
    {
        double dummy;
        effectiveInvViscBoth(c, visc, inv_mu_w_eff, dummy, false);
    }

    void PolymerProperties::effectiveInvViscWithDer(const double c, const double visc,
                                                 double& inv_mu_w_eff,
                                                 double& dinv_mu_w_eff_dc) const {
        effectiveInvViscBoth(c, visc, inv_mu_w_eff, dinv_mu_w_eff_dc, true);
    }

    void PolymerProperties::effectiveInvViscBoth(const double c, const double visc,
                                                 double& inv_mu_w_eff,
                                                 double& dinv_mu_w_eff_dc,
                                                 bool if_with_der) const {
        const double cbar = c/c_max_;
        const double mu_w = visc;
        double mu_m;
        const double omega = mix_param_;
        double dmu_m_dc;
        if (if_with_der) {
            mu_m = viscMultWithDer(c, &dmu_m_dc)*mu_w;
            dmu_m_dc *= mu_w;
        } else {
            mu_m = viscMult(c)*mu_w;
        }
        double mu_p = viscMult(c_max_)*mu_w;
        double inv_mu_m_omega = std::pow(mu_m, -omega);
        double inv_mu_w_e   = inv_mu_m_omega*std::pow(mu_w, omega - 1.);
        double inv_mu_p_eff = inv_mu_m_omega*std::pow(mu_p, omega - 1.);
        inv_mu_w_eff = (1.0 - cbar)*inv_mu_w_e + cbar*inv_mu_p_eff;
        if (if_with_der) {
            double dinv_mu_w_e_dc = -omega*dmu_m_dc*std::pow(mu_m, -omega - 1)*std::pow(mu_w, omega - 1);
            double dinv_mu_p_eff_dc = -omega*dmu_m_dc*std::pow(mu_m, -omega - 1)*std::pow(mu_p, omega - 1);
            dinv_mu_w_eff_dc = (1 - cbar)*dinv_mu_w_e_dc + cbar*dinv_mu_p_eff_dc +
                1/c_max_*(inv_mu_p_eff - inv_mu_w_e);
        }
    }

    void PolymerProperties::effectiveInvPolyVisc(const double c,
                                                 const double visc,
                                                 double& inv_mu_p_eff) const
    {
        double dummy;
        effectiveInvPolyViscBoth(c, visc, inv_mu_p_eff, dummy, false);

    }

    void PolymerProperties::effectiveInvPolyViscWithDer(const double c,
                                                        const double visc,
                                                        double& inv_mu_p_eff,
                                                        double& d_inv_mu_p_eff_dc) const
    {
        effectiveInvPolyViscBoth(c, visc, inv_mu_p_eff, d_inv_mu_p_eff_dc, true);

    }

    void PolymerProperties::effectiveInvPolyViscBoth(const double c,
                                                     const double visc,
                                                     double& inv_mu_p_eff,
                                                     double& dinv_mu_p_eff_dc,
                                                     const bool if_with_der) const
    {
        const double omega = mix_param_;
        const double mu_w = visc;

        double mu_m = 0.0;
        double dmu_m_dc = 0.0;

        if (if_with_der) {
            mu_m = viscMultWithDer(c, &dmu_m_dc)*mu_w;
            dmu_m_dc *= mu_w;
        } else {
            mu_m = viscMult(c)*mu_w;
        }

        const double inv_mu_m_omega = std::pow(mu_m, -omega);
        const double mu_p = viscMult(c_max_) * mu_w;
        inv_mu_p_eff = inv_mu_m_omega * std::pow(mu_p, omega - 1.);

        if (if_with_der) {
            dinv_mu_p_eff_dc = -omega * dmu_m_dc * std::pow(mu_m, -omega - 1) * std::pow(mu_p, omega - 1);
        }
    }

    void PolymerProperties::effectiveRelperm(const double c,
                                             const double cmax,
                                             const double* relperm,
                                             double& eff_relperm_wat) const {
        double dummy;
        effectiveRelpermBoth(c, cmax, relperm, 0, eff_relperm_wat,
                             dummy, dummy, false);
    }

    void PolymerProperties::effectiveRelpermWithDer (const double c,
                                                     const double cmax,
                                                     const double* relperm,
                                                     const double* drelperm_ds,
                                                     double& eff_relperm_wat,
                                                     double& deff_relperm_wat_ds,
                                                     double& deff_relperm_wat_dc) const {
        effectiveRelpermBoth(c, cmax, relperm,
                             drelperm_ds, eff_relperm_wat,
                             deff_relperm_wat_ds, deff_relperm_wat_dc,
                             true);
    }

    void PolymerProperties::effectiveRelpermBoth(const double c,
                                                 const double cmax,
                                                 const double* relperm,
                                                 const double* drelperm_ds,
                                                 double& eff_relperm_wat,
                                                 double& deff_relperm_wat_ds,
                                                 double& deff_relperm_wat_dc,
                                                 bool if_with_der) const {
        double c_ads;
        double dc_ads_dc;
        adsorptionBoth(c, cmax, c_ads, dc_ads_dc, if_with_der);
        double rk = 1 + (res_factor_ - 1)*c_ads/c_max_ads_;
        eff_relperm_wat = relperm[0]/rk;
        if (if_with_der) {
            deff_relperm_wat_ds = (drelperm_ds[0]-drelperm_ds[2])/rk; //derivative with respect to sw
            //\frac{\partial k_{rw_eff}}{\parital c} = -\frac{krw}{rk^2}\frac{(RRF-1)}{c^a_{max}}\frac{\partial c^a}{\partial c}.
            deff_relperm_wat_dc = -(res_factor_ - 1)*dc_ads_dc*relperm[0]/(rk*rk*c_max_ads_);
        } else {
            deff_relperm_wat_ds = -1.0;
            deff_relperm_wat_dc = -1.0;
        }
    }

    void PolymerProperties::effectiveMobilities(const double c,
                                                const double cmax,
                                                const double* visc,
                                                const double* relperm,
                                                double* mob) const
    {
        double dummy;
        double dummy_pointer[4];
        effectiveMobilitiesBoth(c, cmax, visc, relperm,
                                dummy_pointer, mob, dummy_pointer, dummy, false);
    }


    void PolymerProperties::effectiveMobilitiesWithDer(const double c,
                                                       const double cmax,
                                                       const double* visc,
                                                       const double* relperm,
                                                       const double* drelpermds,
                                                       double* mob,
                                                       double* dmobds,
                                                       double& dmobwatdc) const
    {
        effectiveMobilitiesBoth(c, cmax, visc,
                                relperm, drelpermds, mob, dmobds,
                                dmobwatdc, true);
    }

    void PolymerProperties::effectiveMobilitiesBoth(const double c,
                                                    const double cmax,
                                                    const double* visc,
                                                    const double* relperm,
                                                    const double* drelperm_ds,
                                                    double* mob,
                                                    double* dmob_ds,
                                                    double& dmobwat_dc,
                                                    bool if_with_der) const
    {
        double inv_mu_w_eff;
        double dinv_mu_w_eff_dc;
        effectiveInvViscBoth(c, visc[0], inv_mu_w_eff, dinv_mu_w_eff_dc, if_with_der);
        double eff_relperm_wat;
        double deff_relperm_wat_ds;
        double deff_relperm_wat_dc;

        effectiveRelpermBoth(c, cmax, relperm,
                             drelperm_ds, eff_relperm_wat,
                             deff_relperm_wat_ds, deff_relperm_wat_dc,
                             if_with_der);

        // The "function" eff_relperm_wat is defined as a function of only sw (so that its
        // partial derivative with respect to so is zero).

        mob[0] = eff_relperm_wat*inv_mu_w_eff;
        mob[1] = relperm[1]/visc[1];

        if (if_with_der) {
            dmobwat_dc = eff_relperm_wat*dinv_mu_w_eff_dc
                + deff_relperm_wat_dc*inv_mu_w_eff;
            dmob_ds[0*2 + 0] = deff_relperm_wat_ds*inv_mu_w_eff;
	    // one have to deside which variables to derive
	    // here the full derivative is written out
            dmob_ds[0*2 + 1] = 0.0*(drelperm_ds[0*2 + 1] - drelperm_ds[1*2 + 1])/visc[1];
            dmob_ds[1*2 + 0] = -0.0*deff_relperm_wat_ds*inv_mu_w_eff;
            dmob_ds[1*2 + 1] = (drelperm_ds[1*2 + 1] - drelperm_ds[0*2 + 1])/visc[1];
        }
    }

    void PolymerProperties::effectiveTotalMobility(const double c,
                                                   const double cmax,
                                                   const double* visc,
                                                   const double* relperm,
                                                   double& totmob) const
    {
        double dummy1[4];
        double dummy2[2];
        effectiveTotalMobilityBoth(c, cmax, visc, relperm, dummy1,
                                   totmob, dummy2, false);
    }

    void PolymerProperties::effectiveTotalMobilityWithDer(const double c,
                                                          const double cmax,
                                                          const double* visc,
                                                          const double* relperm,
                                                          const double* drelperm_ds,
                                                          double& totmob,
                                                          double* dtotmob_dsdc) const
    {
        effectiveTotalMobilityBoth(c, cmax, visc, relperm, drelperm_ds,
                                   totmob, dtotmob_dsdc, true);
    }

    void PolymerProperties::effectiveTotalMobilityBoth(const double c,
                                                       const double cmax,
                                                       const double* visc,
                                                       const double* relperm,
                                                       const double* drelperm_ds,
                                                       double& totmob,
                                                       double* dtotmob_dsdc,
                                                       bool if_with_der) const
    {
        double mob[2];
        double dmob_ds[4];
        double dmobwat_dc;
        effectiveMobilitiesBoth(c, cmax, visc, relperm, drelperm_ds,
                                mob, dmob_ds, dmobwat_dc, if_with_der);
        totmob = mob[0] + mob[1];
        if (if_with_der) {
            dtotmob_dsdc[0] = dmob_ds[0*2 + 0] -  dmob_ds[1*2 + 0]
                + dmob_ds[0*2 + 1] - dmob_ds[1*2 + 1]; //derivative with respect to sw
            dtotmob_dsdc[1] = dmobwat_dc; //derivative with respect to c
        }
    }

    void PolymerProperties::computeMc(const double& c, double& mc) const
    {
        double dummy;
        computeMcBoth(c, mc, dummy, false);
    }

    void PolymerProperties::computeMcWithDer(const double& c, double& mc,
                                             double& dmc_dc) const
    {
        computeMcBoth(c, mc, dmc_dc, true);
    }

    void PolymerProperties::computeMcBoth(const double& c, double& mc,
                                          double& dmc_dc, bool if_with_der) const
    {
        double cbar = c/c_max_;
        double omega = mix_param_;
        double r = std::pow(viscMult(c_max_), 1 - omega); // viscMult(c_max_)=mu_p/mu_w
        mc = c/(cbar + (1 - cbar)*r);
        if (if_with_der) {
            dmc_dc = r/std::pow(cbar + (1 - cbar)*r, 2);
        } else {
            dmc_dc = 0.;
        }
    }

    bool PolymerProperties::computeShearMultLog(std::vector<double>& water_vel, std::vector<double>& visc_mult, std::vector<double>& shear_mult) const
    {

        double refConcentration = plyshlogRefConc();
        double refViscMult = viscMult(refConcentration);

        std::vector<double> shear_water_vel = shearWaterVelocity();
        std::vector<double> shear_vrf = shearViscosityReductionFactor();

        std::vector<double> logShearWaterVel;
        std::vector<double> logShearVRF;

        logShearWaterVel.resize(shear_water_vel.size());
        logShearVRF.resize(shear_water_vel.size());

        // converting the table using the reference condition
        for (size_t i = 0; i < shear_vrf.size(); ++i) {
            shear_vrf[i] = (refViscMult * shear_vrf[i] - 1.) / (refViscMult - 1);
            logShearWaterVel[i] = std::log(shear_water_vel[i]);
        }

        shear_mult.resize(water_vel.size());

        // the mimum velocity to apply the shear-thinning
        const double minShearVel = shear_water_vel[0];
        const double maxShearVel = shear_water_vel.back();
        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

        for (size_t i = 0; i < water_vel.size(); ++i) {

            if (visc_mult[i] - 1. < epsilon || std::abs(water_vel[i]) < minShearVel) {
                 shear_mult[i] = 1.0;
                 continue;
            }

            for (size_t j = 0; j < shear_vrf.size(); ++j) {
                logShearVRF[j] = (1 + (visc_mult[i] - 1.0) * shear_vrf[j]) / visc_mult[i];
                logShearVRF[j] = std::log(logShearVRF[j]);
            }

            // const double logWaterVelO = std::log(water_vel[i]);
            const double logWaterVelO = std::log(std::abs(water_vel[i]));

            size_t iIntersection; // finding the intersection on the iIntersectionth table segment
            bool foundSegment = false;

            for (iIntersection = 0; iIntersection < shear_vrf.size() - 1; ++iIntersection) {

                double temp1 = logShearVRF[iIntersection] + logShearWaterVel[iIntersection] - logWaterVelO;
                double temp2 = logShearVRF[iIntersection + 1] + logShearWaterVel[iIntersection + 1] - logWaterVelO;

                // ignore the cases the temp1 or temp2 is zero first for simplicity.
                // several more complicated cases remain to be implemented.
                if( temp1 * temp2 < 0.){
                    foundSegment = true;
                    break;
                }
            }

            if (foundSegment == true) {
                detail::Point2D lineSegment[2];
                lineSegment[0] = detail::Point2D{logShearWaterVel[iIntersection], logShearVRF[iIntersection]};
                lineSegment[1] = detail::Point2D{logShearWaterVel[iIntersection + 1], logShearVRF[iIntersection + 1]};

                detail::Point2D line[2];
                line[0] = detail::Point2D{0, logWaterVelO};
                line[1] = detail::Point2D{logWaterVelO, 0};

                detail::Point2D intersectionPoint;

                bool foundIntersection = detail::Point2D::findIntersection(lineSegment, line, intersectionPoint);

                if (foundIntersection) {
                    shear_mult[i] = std::exp(intersectionPoint.getY());
                } else {
                    std::cerr << " failed in finding the solution for shear-thinning multiplier " << std::endl;
                    return false; // failed in finding the solution.
                }
            } else {
                // check if the failure in finding the shear multiplier is due to too big water velocity.
                if ((logWaterVelO - logShearVRF.back()) < logShearWaterVel.back()) {
                    std::cout << " the veclocity is " << water_vel[i] << std::endl;
                    std::cout << " max shear velocity is " << maxShearVel << std::endl;
                    std::cerr << " something wrong happend in finding segment" << std::endl;
                    return false;
                } else {
                    shear_mult[i] = std::exp(logShearVRF.back());
                }
            }

        }

        return true;
    }
}
